#include <set>
#include <gmsh.h>

int main(int argc, char **argv) {

  gmsh::initialize();
  gmsh::model::add("Tokamak");

  double lc = 1.5;

/* ТОР */

  double R1 = 3.0; //внутренний радиус тора
  double R2 = R1 + 0.1; // внешний радуис тора

/* центр тора */
  double y1 = 14.0;
  double x1 = 0.0;
  double z1 = 0.0;

/* начало координат и центр тора */
  gmsh::model::occ::addPoint(0, 0, 0, lc, 1);
  gmsh::model::occ::addPoint(x1, y1, z1, lc, 2); 

/* задаем тор */
  std::vector<double> angles = {0, M_PI / 2, M_PI, 3 * M_PI / 2};
  int ang_size = static_cast<int>(angles.size());
  int innerStartTag = 3;
  int outerStartTag = innerStartTag + angles.size();

/* внутрення окружность */
  for (std::size_t i = 0; i < angles.size(); ++i) {
    double ang = angles[i];
    gmsh::model::occ::addPoint(x1, y1 + R1 * cos(ang), z1 + R1 * sin(ang), lc, innerStartTag + static_cast<int>(i));
  }

/* дуги внутренней окружности, у них id = 100-103 */
  for (int i = 0; i < ang_size; ++i) {
    int start_tag = innerStartTag + i;
    int end_tag = innerStartTag + (i + 1) % ang_size;
    int id = 100 + i;
    gmsh::model::occ::addCircleArc(start_tag, 2, end_tag, id);
  }

/* внешняя окружность */
  for (std::size_t i = 0; i < angles.size(); ++i) {
    double ang = angles[i];
    gmsh::model::occ::addPoint(x1, y1 + R2 * cos(ang), z1 + R2 * sin(ang), lc, outerStartTag + static_cast<int>(i));
  }

/* дуги внешней окружности, у них id = 200-203 */
  for (int i = 0; i < ang_size; ++i) {
    int start_tag = outerStartTag + i;
    int end_tag = outerStartTag + (i + 1) % ang_size;
    int id = 200 + i;
    gmsh::model::occ::addCircleArc(start_tag, 2, end_tag, id);
  }

/* добавляем контур поверхностей */
  gmsh::model::occ::addCurveLoop({200, 201, 202, 203}, 300); // внешний
  gmsh::model::occ::addCurveLoop({100, 101, 102, 103}, 301); // внутренний
  gmsh::model::occ::addPlaneSurface({300, -301}, 400);

/* крутим вокург оси z на угол 2pi, сначала делала на geo, но там не получилось, а в occ уже можно на 2pi вращать */
  std::vector<std::pair<int,int>> out1;
  gmsh::model::occ::revolve({{2, 400}}, 0, 0, 0, 0, 0, 1, 2 * M_PI, out1);
  int torusTag = out1[1].second;
  gmsh::model::occ::removeAllDuplicates();

/* ИНДУКТОР */

  double R3 = 2.5; //внутренний радиус тора
  double R4 = R1 + 2.0; // внешний радуис тора
  
  for (std::size_t i = 0; i < angles.size(); ++i) {
    double ang = angles[i];
    gmsh::model::occ::addPoint(R3 * cos(ang), R3 * sin(ang), 0, lc, 1000 + static_cast<int>(i));
  }

  for (int i = 0; i < ang_size; ++i) {
    int start_tag = 1000 + i;
    int end_tag = 1000 + (i + 1) % ang_size;
    int id = 1100 + i;
    gmsh::model::occ::addCircleArc(start_tag, 1, end_tag, id);
  }

  for (std::size_t i = 0; i < angles.size(); ++i) {
    double ang = angles[i];
    gmsh::model::occ::addPoint(R4 * cos(ang), R4 * sin(ang), 0, lc, 2000 + static_cast<int>(i));
  }

  for (int i = 0; i < ang_size; ++i) {
    int start_tag = 2000 + i;
    int end_tag = 2000 + (i + 1) % ang_size;
    int id = 2100 + i;
    gmsh::model::occ::addCircleArc(start_tag, 1, end_tag, id);
  }

  gmsh::model::occ::addCurveLoop({2100, 2101, 2102, 2103}, 3000); // внешний
  gmsh::model::occ::addCurveLoop({1100, 1101, 1102, 1103}, 3001); // внутренний
  gmsh::model::occ::addPlaneSurface({3000, -3001}, 3003);
  gmsh::model::occ::translate({{2, 3003}}, 0, 0, -2);
  std::vector<std::pair<int, int>> out2;
  gmsh::model::occ::extrude({{2, 3003}}, 0, 0, 4, out2);
  int coilTag  = out2[1].second;

/* ПОЛОИДАЛЬНЫЕ КАТУШКИ */

  double a = 1.5;
  double b = 1.0;

  int xp = 0;
  int yp = 23;
  int zp = 7;

  int p1 = gmsh::model::occ::addPoint(xp, yp, zp, lc);
  int p2 = gmsh::model::occ::addPoint(xp, yp + b, zp, lc);
  int p3 = gmsh::model::occ::addPoint(xp, yp + b, zp + a, lc);
  int p4 = gmsh::model::occ::addPoint(xp, yp, zp + a, lc);

  int l1 = gmsh::model::occ::addLine(p1, p2);
  int l2 = gmsh::model::occ::addLine(p2, p3);
  int l3 = gmsh::model::occ::addLine(p3, p4);
  int l4 = gmsh::model::occ::addLine(p4, p1);

  int loop = gmsh::model::occ::addCurveLoop({l1, l2, l3, l4});
  int surf = gmsh::model::occ::addPlaneSurface({loop});
  std::vector<std::pair<int,int>> out3;
  gmsh::model::occ::revolve({{2, surf}}, 0, 0, 0, 0, 0, 1, 2 * M_PI, out3);

  int poloidVolume = -1;
  for (auto &e : out3) {
    if (e.first == 3) {
      poloidVolume = e.second;
    }
  }

  std::vector<std::pair<int,int>> polCopy;
  gmsh::model::occ::copy({{3, poloidVolume}}, polCopy);
  gmsh::model::occ::mirror(polCopy, 0, 0, 1, 0);

  int poloidVolumeMirror = -1;
  for (auto &e : polCopy) {
    if (e.first == 3) poloidVolumeMirror = e.second;
  }

  gmsh::model::occ::removeAllDuplicates();

/* ТОРОИДАЛЬНЫЕ КАТУШКИ */

  double c = 2.0;
  double d = 1.0;

  int xt = 14;
  int yt = -c/2;
  int zt = 6;

  int p5 = gmsh::model::occ::addPoint(xt, yt, zt, lc);
  int p6 = gmsh::model::occ::addPoint(xt, yt + c, zt, lc);
  int p7 = gmsh::model::occ::addPoint(xt, yt + c, zt + d, lc);
  int p8 = gmsh::model::occ::addPoint(xt, yt, zt + d, lc);

  int l5 = gmsh::model::occ::addLine(p5, p6);
  int l6 = gmsh::model::occ::addLine(p6, p7);
  int l7 = gmsh::model::occ::addLine(p7, p8);
  int l8 = gmsh::model::occ::addLine(p8, p5);

  int loop_t = gmsh::model::occ::addCurveLoop({l5, l6, l7, l8});
  int surf_t = gmsh::model::occ::addPlaneSurface({loop_t});
  std::vector<std::pair<int,int>> out_t1;
  gmsh::model::occ::revolve({{2, surf_t}}, 14, 0, 0, 0, 1, 0, 2 * M_PI, out_t1);

  int torCoilVol = -1;
  for (auto &e : out_t1) {
    if (e.first == 3) { 
      torCoilVol = e.second; 
      break; 
    }
  }

  const int N = 8;
  std::vector<std::pair<int,int>> torCoils;
  torCoils.push_back({3, torCoilVol});

  for (int k = 1; k < N; ++k) {
    std::vector<std::pair<int,int>> cp;
    gmsh::model::occ::copy({{3, torCoilVol}}, cp);
    double ang = 2.0 * M_PI * k / N;
    gmsh::model::occ::rotate(cp, 0, 0, 0, 0, 0, 1, ang);
    for (auto &e : cp) {
      if (e.first == 3) { 
        torCoils.push_back(e);
      }
    }
  }

  gmsh::model::occ::removeAllDuplicates();

/* синхронизация, генерация сетки, вывод*/
  gmsh::model::occ::synchronize();
  gmsh::model::mesh::generate(3);

  gmsh::model::setColor({{3, poloidVolume}}, 0, 0, 255);
  gmsh::model::setColor({{3, poloidVolumeMirror}}, 0, 0, 255);
 
  gmsh::model::setColor({{3, torusTag}}, 187, 195, 255);
  gmsh::model::setColor({{3, coilTag}}, 255, 255, 128);

  for (auto &v : torCoils) {
    gmsh::model::setColor({v}, 255, 0, 0);
  }
 
  std::set<std::string> args(argv, argv + argc);
  if (!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}