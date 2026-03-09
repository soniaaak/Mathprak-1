import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem.petsc import LinearProblem
import ufl

# T - финальный момент моделирования,
# num_steps - число шагов по времени,
# dt - шаг по времени,
# alpha - коэффициент диффузии (теплопроводности). Чем он выше, тем быстрее все "размажется", другого слова не придумала
T = 12.0
num_steps = 700
dt = T / num_steps
alpha = 0.025

# Делаем сетку квадрата и задаем тип функций на ней
nx = ny = 170
domain = mesh.create_unit_square(MPI.COMM_WORLD, nx, ny)
V = fem.functionspace(domain, ("CG", 1))

# На границе температура равна 0.
u_D = fem.Function(V)
u_D.interpolate(lambda x: np.zeros(x.shape[1]))

# Ищем границу и точки сетки на границе
fdim = domain.topology.dim - 1 # Размерность границы (для квадрата это линии)
facets = mesh.locate_entities_boundary(
    domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool)
) # Все части внешней границы.
dofs = fem.locate_dofs_topological(V, fdim, facets) # Точки сетки на границе
bc = fem.dirichletbc(u_D, dofs) # В этих точках ставим температуру u_D (в нашей задаче это 0)

u_n = fem.Function(V) # Начальная температура

def initial_shape(x):
    # Координаты точек сетки
    X = x[0]
    Y = x[1]
    
    # Удобная математическая функция для задания локальных пятен
    def bump(x0, y0, sx, sy, a):
        return a * np.exp(-(((X - x0) / sx) ** 2 + ((Y - y0) / sy) ** 2))

    # Собираем "кошачью мордочку" из пятен
    face = bump(0.50, 0.52, 0.17, 0.17, 7.5)
    ears = bump(0.38, 0.70, 0.05, 0.11, 5.0) + bump(0.62, 0.70, 0.05, 0.11, 5.0)
    cheeks = bump(0.42, 0.47, 0.07, 0.05, 2.5) + bump(0.58, 0.47, 0.07, 0.05, 2.5)
    chin = bump(0.50, 0.38, 0.08, 0.05, 2.0)

    eyes = bump(0.43, 0.55, 0.03, 0.03, 4.0) + bump(0.57, 0.55, 0.03, 0.03, 4.0)
    nose = bump(0.50, 0.49, 0.025, 0.02, 3.5)
    mouth = bump(0.47, 0.45, 0.03, 0.02, 2.5) + bump(0.53, 0.45, 0.03, 0.02, 2.5)

    glow = bump(0.50, 0.52, 0.30, 0.22, 0.7)
    value = face + ears + cheeks + chin + glow - eyes - nose - mouth
    
    return np.maximum(value, 0.0)

# Записываем начальную температуру в нашу сетку.
u_n.interpolate(initial_shape)
u_n.name = "temperature"

# Готовим формулу для одного шага по времени
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(domain, 0.0)

F = u * v * ufl.dx + dt * alpha * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx - (u_n + dt * f) * v * ufl.dx
a = ufl.lhs(F)
L = ufl.rhs(F)

# a и L — это левая и правая части вариационной формы, из которых собираются матрица и правая часть линейной системы
# bcs — граничные условия
# petsc_options задают тип линейного решателя: тут испольузем прямое LU-разложение
problem = LinearProblem(
    a, L, bcs=[bc],
    petsc_options_prefix="heat_",
    petsc_options={"ksp_type": "preonly", "pc_type": "lu"}
)

# Файл, куда сохраняем результат для просмотра (например, в ParaView).
xdmf = io.XDMFFile(domain.comm, "heat_cat_face.xdmf", "w")
xdmf.write_mesh(domain)
xdmf.write_function(u_n, 0.0)

# Решаем шаг, сохраняем результат, обновляем прошлое значение.
t = 0.0
for n in range(num_steps):
    t += dt
    uh = problem.solve()
    uh.name = "temperature"
    xdmf.write_function(uh, t)
    u_n.x.array[:] = uh.x.array # Новое решение становится старым для следующего шага

xdmf.close()
