from gurobipy import *
import numpy as np

# Constant creation
f = [400, 414, 326]
a = [18, 25, 20]
C = [[22, 33, 24],
     [33, 23, 30],
     [20, 25, 27], ]
D = [206 + 40, 274 + 40, 220 + 40]
dl = [206, 274, 220]
du = [40, 40, 40]
k = 0  # Iterative counting

# Create model
MP = Model()  # Master-problem
SP = Model()  # Sub-problem(KKT)
SDSP = Model()  # Sub-problem (strong duality)
# Construction of Master-problem
# addVars
y = MP.addVars(len(f), lb=0, ub=1, obj=f, vtype=GRB.INTEGER, name='y')
z = MP.addVars(len(a), lb=0, obj=a, vtype=GRB.CONTINUOUS, name='z')
g = MP.addVars(3, lb=0, ub=1.0, name='g')
η = MP.addVar(obj=1.0, name='η')

# addConstrs
Column1 = MP.addConstrs((z[i] <= 800 * y[i] for i in range(3)), name='column1')
Column4 = MP.addConstr(quicksum(z[i] for i in range(3)) >= 772, name='z')
Column5 = MP.addConstr(quicksum(g[i] for i in range(2)) <= 1.2, name='column5')
Column6 = MP.addConstr(quicksum(g[i] for i in range(3)) <= 1.8, name='column6')

# MP.write("MP.lp")  # model print and visual inspection model,can open it with Notepad++
MP.optimize()  # Solve Model
LB = MP.objval  # get optimum value of model
# 添加变量
Mx = np.zeros((3, 3))
Mλ = np.zeros((3))
Mπ = np.zeros((3))
for i in range(3):
    for j in range(3):
        Mx[i][j] = min(D[j], z[i].x)
        Mλ[i] = max(C[i][0], C[i][1], C[i][2])
        Mπ[i] = max(C[0][i], C[1][i], C[2][i])
# 子问题求解kkt
x = SP.addVars(3, 3, lb=0, obj=np.array(C) * -1, vtype=GRB.CONTINUOUS, name='x')
g = SP.addVars(3, lb=0, ub=1.0, name='g')
d = [206 + 40 * g[0], 274 + 40 * g[1], 220 + 40 * g[2]]
α = SP.addVars(3, 3, vtype=GRB.BINARY, name='α')
β = SP.addVars(3, vtype=GRB.BINARY, name='β')
γ = SP.addVars(3, vtype=GRB.BINARY, name='γ')
λ = SP.addVars(3, vtype=GRB.CONTINUOUS, name='λ')
π = SP.addVars(3, vtype=GRB.CONTINUOUS, name='π')
A = [252, 0, 520]
S1 = SP.addConstrs(((quicksum(x[i, j] for j in range(3))) <= z[i].x for i in range(3)), name='SPcolumn1') #0:3
S2 = SP.addConstrs(((quicksum(x[i, j] for i in range(3))) >= d[j] for j in range(3)), name='SPcolumn2') #3:6
S3 = SP.addConstrs(((λ[j] - π[i]) <= C[i][j] for i in range(3) for j in range(3)), name='SPcolumn3') #6:15
S4 = SP.addConstrs((Mx[i][j] * α[i, j] >= x[i, j] for i in range(3) for j in range(3)), name='SPcolumn4') #15:24
S5 = SP.addConstrs(
    ((C[i][j] - λ[j] + π[i]) <= (C[i][j] + Mπ[i]) * (1 - α[i, j]) for i in range(3) for j in range(3)),
    name='SPcolumn5') #24:32
S6 = SP.addConstrs((λ[j] <= Mλ[j] * β[j] for j in range(3)), name='SPcolumn6') #32:35
S7 = SP.addConstrs(((quicksum(x[i, j] for i in range(3)) - d[j]) <= 40 * (1 - β[j]) for j in range(3)),
                   name='SPcolumn7') #35:44
S8 = SP.addConstrs((π[i] <= Mπ[i] * γ[i] for i in range(3)), name='SPcolumn8')
S9 = SP.addConstrs(((z[i].x - quicksum(x[i, j] for j in range(3))) <= (1 - γ[i]) * z[i].x for i in range(3)),
                   name='SPcolumn9')
SP.addConstr(quicksum(g[i] for i in range(2)) <= 1.2, name='SP10')
SP.addConstr(quicksum(g[i] for i in range(3)) <= 1.8, name='SP11')
SP.write("SP.lp")
SP.optimize()
d = [dl[i] + du[i] * g[i].x for i in range(3)]
Q = SP.objval
UB = LB - η.x - Q
while UB - LB > 10e-4:
    xx = MP.addVars(3, 3, lb=0, vtype=GRB.CONTINUOUS, name='x')
    Column2 = MP.addConstrs(((quicksum(xx[i, j] for j in range(3))) <= z[i] for i in range(3)), name='column2')
    Column3 = MP.addConstrs(((quicksum(xx[i, j] for i in range(3))) >= d[j] for j in range(3)), name='column3')
    Column7 = MP.addConstr(quicksum(C[i][j] * xx[i, j] for i in range(3) for j in range(3)) <= η)
    MP.optimize()
    LB = MP.objval
    SP.remove(SP.getConstrs()[0:6])
    SP.remove(SP.getConstrs()[15:23])
    SP.remove(SP.getConstrs()[36:39])
    SP.remove(SP.getConstrs()[42:45])
    S1 = SP.addConstrs(((quicksum(x[i, j] for j in range(3))) <= z[i].x for i in range(3)), name='SPcolumn1')
    S2 = SP.addConstrs(((quicksum(x[i, j] for i in range(3))) >= d[j] for j in range(3)), name='SPcolumn2')
    S7 = SP.addConstrs(((quicksum(x[i, j] for i in range(3)) - d[j]) <= 40 * (1 - β[j]) for j in range(3)),
                       name='SPcolumn7')
    S9 = SP.addConstrs(((z[i].x - quicksum(x[i, j] for j in range(3))) <= (1 - γ[i]) * z[i].x for i in range(3)),
                       name='SPcolumn9')
    for i in range(3):
        for j in range(3):
            Mx[i][j] = min(D[j], z[i].x)
    S4 = SP.addConstrs((Mx[i][j] * α[i, j] >= x[i, j] for i in range(3) for j in range(3)), name='SPcolumn4')
    SP.write('SP.lp')
    SP.optimize()
    UB = LB - η.x - SP.objval
    k = k + 1
print(LB)
