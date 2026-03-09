#include <set>
#include <cmath>
#include <sstream> // Нужно для сборки длинной строки формулы по частям
#include <iomanip> // Нужно для точности при переводе double -> string
#include <string>
#include <vector>
#include <array>
#include <gmsh.h>

int main(int argc, char **argv)
{
    gmsh::initialize();
    gmsh::model::add("gear_mesh");

    try {
        gmsh::merge("../Gear.stl");
    } catch (...) {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return 0;
    }

    const double angle = 30;
    const bool forceParametrizablePatches = false;
    const bool includeBoundary = true;
    const double curveAngle = 180;

    gmsh::model::mesh::classifySurfaces(
        angle * M_PI / 180.0,
        includeBoundary,
        forceParametrizablePatches,
        curveAngle * M_PI / 180.0
    );
    gmsh::model::mesh::createGeometry();

    std::vector<std::pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);

    std::vector<int> surfaceTags;
    for(const auto &surf : surfaces) surfaceTags.push_back(surf.second);

    int loop = gmsh::model::geo::addSurfaceLoop(surfaceTags);
    gmsh::model::geo::addVolume({loop});
    gmsh::model::geo::synchronize();

    // Получаем ограничивающий прямоугольник модели
    double xmin, ymin, zmin, xmax, ymax, zmax;
    gmsh::model::getBoundingBox(-1, -1, xmin, ymin, zmin, xmax, ymax, zmax);

    // Центр модели и характерные размеры
    const double cx = 0.5 * (xmin + xmax);
    const double cy = 0.5 * (ymin + ymax);
    const double cz = 0.5 * (zmin + zmax);
    const double R = 0.5 * std::max(xmax - xmin, ymax - ymin);
    const double H = std::max(zmax - zmin, 1e-9);

    // Здесь задаем параметры зон, где хотим изменить плотность сетки
    const double rTeeth = 0.88 * R; // зона зубьев (внешнее кольцо)
    const double wTeeth = 0.09 * R; // ширина сгущения у зубьев
    const double rHub = 0.42 * R; // зона ступицы (внутреннее кольцо)
    const double wHub = 0.16 * R; // ширина сгущения у ступицы
 
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);

    int field = gmsh::model::mesh::field::add("MathEval");

    // Перевод число в строку с опредленной точностью для формулы
    auto num = [](double v) {
        std::ostringstream os;
        os << std::setprecision(12) << v;
        return os.str();
    };

    // Функциия сдвига (удобно вычитать одну координату из другой)
    auto shift = [&](const char *var, double center) {
        return std::string("(") + var + "-(" + num(center) + "))";
    };

    // Считаем квадрат выражения s*s.
    auto sq = [](const std::string &s) {
        return "(" + s + ")*(" + s + ")";
    };

    std::ostringstream expr; // Для того, чтобы в одну строку написать всю формулу для поля
    const std::string dx = shift("x", cx); // x относительно центра
    const std::string dy = shift("y", cy); // y относительно центра
    const std::string dz = shift("z", cz); // z относительно центра
    const std::string rad = "sqrt(" + sq(dx) + " + " + sq(dy) + ")"; // радиус в XY
    const std::string teethBand = "(" + rad + "-(" + num(rTeeth) + "))/" + num(wTeeth); // насколько точка близка к зоне зубьев
    const std::string hubBand = "(" + rad + "-(" + num(rHub) + "))/" + num(wHub); // то же, что и выше, но для Hub
    const std::string wave = "sin(14*" + dx + "/" + num(R) + ")*cos(14*" + dy + "/" + num(R) + ")"; // волновая функция для красивой сетки
    const std::string zNorm = dz + "/" + num(0.5 * H); // нормированная координата по толщине

    // Формула поля F(x,y,z)
    expr
        << "0.3 + 0.7*(" << rad << "/" << num(R) << ")"
        << " - 0.2*exp(-( " << teethBand << "*" << teethBand << " ))"
        << " - 0.15*exp(-( " << hubBand << "*" << hubBand << " ))"
        << " + 0.17*((" << wave << ")*(" << wave << "))"
        << " + 0.04*((" << zNorm << ")*(" << zNorm << "))";

    gmsh::model::mesh::field::setString(field, "F", expr.str());
    gmsh::model::mesh::field::setAsBackgroundMesh(field);
    gmsh::model::mesh::generate(3);
    
    std::vector<std::pair<int, int>> surfEntities;
    gmsh::model::getEntities(surfEntities, 2);

    // Пастельные цвета для красоты
    const std::vector<std::array<int, 3>> pastel = {
        {202, 176, 255},
        {186, 146, 255},
        {216, 191, 255},
        {230, 190, 255},
        {176, 224, 230}
    };

    for(size_t i = 0; i < surfEntities.size(); ++i) {
        const auto &c = pastel[i % pastel.size()];
        gmsh::model::setColor({surfEntities[i]}, c[0], c[1], c[2]);
    }

    std::vector<std::pair<int, int>> volEntities;
    gmsh::model::getEntities(volEntities, 3);

    for(const auto &v : volEntities) {
        gmsh::model::setColor({v}, 186, 146, 255);
    }

    gmsh::write("gear_mesh.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}