#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <cassert>
#include <sys/stat.h>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

static const double PI = 3.14159265358979323846;

// Ограничиваем число v интервалом [a, b]
static double clamp(double v, double a, double b) {
    return std::max(a, std::min(b, v));
}

// Сглаженная функция перехода от 0 к 1, была создана, чтобы не было резких переходов, все было плавно
static double smoothstep01(double x) {
    x = clamp(x, 0.0, 1.0);
    return x * x * (3.0 - 2.0 * x);
}

// Функция для того, чтобы были красивые названия у файликов (четыречзначные)
static string pad4(int value) {
    string s = to_string(value);
    while(s.size() < 4) s = "0" + s;
    return s;
}

// Класс расчётной точки
class CalcNode
{
friend class CalcMesh;

protected:
    // Текущие координаты
    double x;
    double y;
    double z;

    // Начальные координаты относительно центра вращения
    double x0;
    double y0;
    double z0;

    // Полярные данные в плоскости XY
    double r0;
    double theta0;
    double w;

    // Поля в точке
    double smth;

    // Скорость
    double vx;
    double vy;
    double vz;

public:
    CalcNode() : x(0.0), y(0.0), z(0.0), x0(0.0), y0(0.0), z0(0.0),
                 r0(0.0), theta0(0.0), w(0.0), smth(0.0), vx(0.0), vy(0.0), vz(0.0)
    {
    }
};

// Класс тетраэдра, элемента сетки
class Element
{
friend class CalcMesh;

protected:
    // Индексы узлов, образующих этот элемент сетки
    unsigned long nodesIds[4];
};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 3D-сетка из расчётных точек
    vector<CalcNode> nodes;
    vector<Element> elements;

    // Центр вращения и максимальный радиус узлов относительно центра
    double cx;
    double cy;
    double cz;
    double rmax;

    // Параметры движения
    double omegaSpin; // базовая угловая скорость (насколько быстро объект или среда вращается вокруг оси)
    double torsionAmp; // амплитуда крутильной деформации
    double torsionFreq; // частота этой деформации

    // Параметры скалярной функции
    double waveFreqMat; // как быстро меняется внутренняя волна со временем
    double waveFreqLab; // как быстро меняется внешняя волна со временем
    int mMat; // сколько повторений у внутреннего узора по кругу;
    int mLab; // сколько повторений у внешнего узора по кругу

public:
    // Конструктор сетки из заданного stl-файла
    CalcMesh(
        const vector<double>& nodesCoord,
        const vector<size_t>& nodeTags,
        const vector<size_t>& tetrsNodeTags,
        double _cx, double _cy, double _cz,
        double _omegaSpin, double _torsionAmp, double _torsionFreq,
        double _waveFreqMat, double _waveFreqLab, int _mMat, int _mLab
    ) :
        cx(_cx), cy(_cy), cz(_cz), rmax(0.0),
        omegaSpin(_omegaSpin), torsionAmp(_torsionAmp), torsionFreq(_torsionFreq),
        waveFreqMat(_waveFreqMat), waveFreqLab(_waveFreqLab), mMat(_mMat), mLab(_mLab)
    {
        // Заполнеям узлы
        nodes.resize(nodesCoord.size() / 3);
        for(size_t i = 0; i < nodes.size(); i++) {
            CalcNode n;
            n.x = nodesCoord[i * 3 + 0];
            n.y = nodesCoord[i * 3 + 1];
            n.z = nodesCoord[i * 3 + 2];

            n.x0 = n.x - cx;
            n.y0 = n.y - cy;
            n.z0 = n.z - cz;

            // Мы считаем полярные координаты в плоскости XY, еще максимальный радиус по всей модели
            n.r0 = sqrt(n.x0 * n.x0 + n.y0 * n.y0);
            n.theta0 = atan2(n.y0, n.x0);
            rmax = max(rmax, n.r0);

            nodes[i] = n;
        }

        // Здесь для каждой точки вычисляется вес w в зависимости от расстояния до центра
        // В центре вес равен 0, а на краю равен 1
        // Мы используем вес, чтобы разные части объекта двигались по-разному
        double r_in = 0.40 * rmax;
        double r_out = 1.00 * rmax;
        for(size_t i = 0; i < nodes.size(); i++) {
            double xr = (nodes[i].r0 - r_in) / (r_out - r_in);
            nodes[i].w = smoothstep01(xr);
        }

        // Делаем отображение "gmsh-тег -> индекс точки"
        // Узлы необязталеьно идут подряд, так что дополнительно это проверяем таким образом
        size_t maxTag = 0;
        for(size_t i = 0; i < nodeTags.size(); i++) {
            maxTag = max(maxTag, nodeTags[i]);
        }
        vector<long> tag2id(maxTag + 1, -1);
        for(size_t i = 0; i < nodeTags.size(); i++) {
            tag2id[nodeTags[i]] = i;
        }

        // Элементы
        assert(tetrsNodeTags.size() % 4 == 0);
        elements.resize(tetrsNodeTags.size() / 4);
        for(size_t i = 0; i < elements.size(); i++) {
            for(size_t j = 0; j < 4; j++) {
                size_t gmshTag = tetrsNodeTags[i * 4 + j];
                elements[i].nodesIds[j] = static_cast<unsigned long>(tag2id[gmshTag]);
            }
        }
    }

    // Считаем сетку на момент времени t
    void doTimeStep(double t) {
        for(size_t i = 0; i < nodes.size(); i++) {
            CalcNode &n = nodes[i];

            // Небольшая крутильная деформация возле ступицы
            double torsion = torsionAmp * sin(2.0 * PI * torsionFreq * t) * (1.0 - n.w);

            // Полный угол, типо общий угол поворота = основное вращение + локальная деформация
            double angle = omegaSpin * t + torsion;

            // Мгновенная угловая скорость (посчитали чисто производную угла по времени)
            double omegaEff = omegaSpin + torsionAmp * (2.0 * PI * torsionFreq) * cos(2.0 * PI * torsionFreq * t) * (1.0 - n.w);

            // Вращаем вокруг оси Z
            double ca = cos(angle);
            double sa = sin(angle);
            double xr = ca * n.x0 - sa * n.y0;
            double yr = sa * n.x0 + ca * n.y0;
            double zr = n.z0;

            n.x = cx + xr;
            n.y = cy + yr;
            n.z = cz + zr;

            // Поле скоростей
            n.vx = -omegaEff * yr;
            n.vy =  omegaEff * xr;
            n.vz = 0.0;

            // Скалярная функция
            double thetaMat = n.theta0;
            double thetaLab = thetaMat + angle;

            double waveMat = sin(mMat * thetaMat - 2.0 * PI * waveFreqMat * t);
            double waveLab = sin(mLab * thetaLab + 2.0 * PI * waveFreqLab * t);

            n.smth = (1.0 - n.w) * waveMat + n.w * waveLab;
        }
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    string snapshot(size_t snapNumber, const string& outDir) {
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");

        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        for(size_t i = 0; i < nodes.size(); i++) {
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);

            double _vel[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
            vel->InsertNextTuple(_vel);

            smth->InsertNextValue(nodes[i].smth);
        }

        unstructuredGrid->SetPoints(dumpPoints);
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->SetVectors(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);
        unstructuredGrid->GetPointData()->SetScalars(smth);

        for(size_t i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId(0, elements[i].nodesIds[0]);
            tetra->GetPointIds()->SetId(1, elements[i].nodesIds[1]);
            tetra->GetPointIds()->SetId(2, elements[i].nodesIds[2]);
            tetra->GetPointIds()->SetId(3, elements[i].nodesIds[3]);
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }
        
        // Создаём снапшот в файле с заданным именем
        string fileName = "gear-step-" + pad4(snapNumber) + ".vtu";
        string filePath = outDir + "/" + fileName;

        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(filePath.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();

        return fileName;
    }
};

int main(int argc, char** argv) {
    mkdir("../gear_output", 0777);

    // Параметры по времени
    const int numSteps = 200;
    const double dt = 0.02;

    // Параметры движения шестерёнки
    const double spinRevsPerSec = 0.5;
    const double omegaSpin = 2.0 * PI * spinRevsPerSec;
    const double torsionAmp = 0.35;
    const double torsionFreq = 1.0;

    // Параметры скалярной функции
    const double waveFreqMat = 1.0; // как быстро меняется внутренняя волна со временем
    const double waveFreqLab = 0.7; // как быстро меняется внешняя волна со временем
    const int mMat = 6; // сколько повторений у внутреннего узора по кругу
    const int mLab = 12; // сколько повторений у внешнего узора по кругу

    const size_t GMSH_TETR_CODE = 4;

    gmsh::initialize();
    gmsh::model::add("gear");

    try {
        gmsh::merge("../Gear.stl");
    } catch(...) {
        cout << "Could not load STL: " << "Gear.stl" << endl;
        gmsh::finalize();
        return -1;
    }

    // Восстановление геометрии по STL
    double angle = 40.0;
    bool forceParametrizablePatches = false;
    bool includeBoundary = true;
    double curveAngle = 180.0;
    gmsh::model::mesh::classifySurfaces(angle * PI / 180.0, includeBoundary, forceParametrizablePatches, curveAngle * PI / 180.0);
    gmsh::model::mesh::createGeometry();

    // Объём по поверхностям
    vector<pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);
    vector<int> surfTags;
    for(size_t i = 0; i < surfaces.size(); i++) {
        surfTags.push_back(surfaces[i].second);
    }
    int sl = gmsh::model::geo::addSurfaceLoop(surfTags);
    gmsh::model::geo::addVolume({sl});
    gmsh::model::geo::synchronize();

    // Подбираем шаг сетки
    double xmin, ymin, zmin, xmax, ymax, zmax;
    gmsh::model::getBoundingBox(-1, -1, xmin, ymin, zmin, xmax, ymax, zmax);
    double diag = sqrt((xmax - xmin) * (xmax - xmin) +
                       (ymax - ymin) * (ymax - ymin) +
                       (zmax - zmin) * (zmax - zmin));
    double meshSize = diag / 45.0;
    if(meshSize <= 0.0) meshSize = 1.0;

    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", to_string(meshSize));
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    gmsh::model::mesh::generate(3);

    // Извлекаем узлы
    vector<size_t> nodeTags;
    vector<double> nodesCoord;
    vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // Извлекаем тетраэдры
    vector<size_t>* tetrsNodesTags = nullptr;
    vector<int> elementTypes;
    vector<vector<size_t>> elementTags;
    vector<vector<size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(size_t i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] == GMSH_TETR_CODE) {
            tetrsNodesTags = &elementNodeTags[i];
            break;
        }
    }

    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " << nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    // Центр вращения
    double cx = 0.5 * (xmin + xmax);
    double cy = 0.5 * (ymin + ymax);
    double cz = 0.5 * (zmin + zmax);

    CalcMesh mesh(nodesCoord, nodeTags, *tetrsNodesTags, cx, cy, cz,
                  omegaSpin, torsionAmp, torsionFreq,
                  waveFreqMat, waveFreqLab, mMat, mLab);

    gmsh::finalize();

    // PVD для временной серии
    ofstream pvd("../gear_output/gear.pvd");

    pvd << "<?xml version=\"1.0\"?>\n"; // стандарт XML
    pvd << "<VTKFile type=\"Collection\">\n"; // говорим ParaView, что это коллекция файлов
    pvd << "<Collection>\n"; // начало списка кадров

    // создаём кадры анимации
    for(int step = 0; step <= numSteps; step++) {
        double t = step * dt;
        mesh.doTimeStep(t); // пересчёт сетки на момент времени t
        string name = mesh.snapshot(step, "../gear_output"); // запись кадра .vtu
        pvd << "<DataSet timestep=\"" << t << "\" file=\"" << name << "\"/>\n"; // добавляем кадр в список временной серии
    }

    pvd << "</Collection>\n</VTKFile>\n"; // конец списка кадров
    cout << "Откройте gear.pvd в ParaView и посмотрите анимацию)" << endl;
    return 0;
}