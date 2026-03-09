#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();
  gmsh::model::add("tokamak_stage1_simple");

  // Тороидальная оболочка = внешний тор - внутренний тор.
  const double R = 6.0; // Большой радиус тора
  const double r_out = 1.00; // Внешний малый радиус
  const double r_in = 0.85; // Внутренний малый радиус
  const double h = r_out - r_in; // Толщина стенки

  int tor_out = gmsh::model::occ::addTorus(0.0, 0.0, 0.0, R, r_out);
  int tor_in = gmsh::model::occ::addTorus(0.0, 0.0, 0.0, R, r_in);

  std::vector<std::pair<int, int>> out;
  std::vector<std::vector<std::pair<int, int>>> outMap;
  gmsh::model::occ::cut({{3, tor_out}}, {{3, tor_in}}, out, outMap);
  gmsh::model::occ::synchronize();

  // Берём получившийся объём камеры и красим его в лавандовый.
  int vol_tag = -1;

  for(const auto &e : out) {
    if(e.first == 3) {
      vol_tag = e.second;
      break;
    }
  }

  if(vol_tag != -1) gmsh::model::setColor({{3, vol_tag}}, 230, 230, 250);

  // По заданию нам нужно 3-4 тетраэдра на толщину стенки, поэтому при h=0.15 нужно примерно поулчить 0.04-0.05.
  // Я буду брать значения между этими двумя заданными границами
  gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 0.04);
  gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.05);

  gmsh::model::mesh::generate(3);

  gmsh::write("tokamak_stage1_simple.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();
  return 0;
}
