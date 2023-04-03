# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 17:21:33 2023

@author: wyx
"""
from ExampleMatrix import *
import numpy as np

MP = Model('MP')

# construct main problem
c = np.array([400,414,326])
a = np.array([18,25,20])
b = np.array([22,33,24,33,23,30,20,25,27])
b_mp = b.reshape(3,3)
dl = [206, 274, 220]
du = [40, 40, 40]
bigM = 10**5
LB = -GRB.INFINITY
UB = GRB.INFINITY
epsilon = 1e-5
k = 1

x = MP.addVars(3,3,vtype=GRB.CONTINUOUS)
y = MP.addVars(3,obj=c,vtype=GRB.BINARY)
z = MP.addVars(3,obj=a,vtype=GRB.CONTINUOUS)
d = MP.addVars(3, lb=0, name='d')
eta = MP.addVar(obj=1,vtype=GRB.CONTINUOUS)
g = MP.addVars(3,ub=1,vtype=GRB.CONTINUOUS)

# 建立主问题约束
MP_Cons_1 = MP.addConstrs((z[i] <= 800*y[i] for i in range(3)))
MP_Cons_2 = MP.addConstr((quicksum(z[i] for i in range(3)) >= 772))
MP.optimize()
MP_obj = MP.ObjVal
S_y = np.array([y[i].x for i in range(len(y))])
S_z = np.array([z[i].x for i in range(len(z))])
LB = max(MP_obj, LB)

SP = Model('SP')
x = SP.addMVar((9,),vtype=GRB.CONTINUOUS,name='x')
pi = SP.addMVar(G.shape[0],vtype=GRB.CONTINUOUS,name='pi')
g_sub = SP.addMVar((3,),ub=1,vtype=GRB.CONTINUOUS,name='g_sub')
v = SP.addMVar((G.shape[0],),vtype=GRB.BINARY,name='v')
w = SP.addMVar((G.shape[1],),vtype=GRB.BINARY,name='w')

G1 = SP.addConstr(G@x>=h-M@g_sub-E@np.concatenate([S_y,S_z]),name="G1")
SP.addConstr(G.T@pi<=b,name='pi')
SP.addConstr(pi<=bigM*v,name='v')
G2 = SP.addConstr(G@x-h+E@np.concatenate([S_y,S_z])+M@g_sub<=bigM*(1-v),name='G2')
SP.addConstr(x<=bigM*w,name='w1')
SP.addConstr(b-G.T@pi<=bigM*(1-w),name='w2')
SP.addConstr(g_sub[:2].sum()<=1.2,name='g1')
SP.addConstr(g_sub.sum()<=1.8,name='g2')
SP.setObjective(b@x,GRB.MAXIMIZE)
SP.optimize()
S_g = np.array([g_sub[i].x for i in range(len(g))])
SP_obj = SP.ObjVal    
UB = min(UB, c@S_y+a@S_z+SP_obj)
while abs(UB-LB)>=epsilon:
    # todo 
    # x_new = MP.addMVar((9,),vtype=GRB.CONTINUOUS,name='x_new')
    # MP.addConstr(eta>=b.T@x_new)
    # MP.addConstr(E[:,:3]@y+E[:,3:]@z+G@x_new>=h-M@S_g)
    # MP.optimize()
    # LB = max(LB, MP.objVal)
    
    # add x^{k+1}
    x_new = MP.addVars(3,3, vtype=GRB.CONTINUOUS)
    # eta>=bTx^{k+1}
    MP.addConstr(eta >= quicksum(x_new[i,j]*b_mp[i][j] for i in range(3) for j in range(3)))
    # Ey+Gx^{k+1}>=h-Mu_{k+1}
    MP.addConstrs((quicksum(x_new[i,j] for j in range(3)) <= z[i] for i in range(3)))
    MP.addConstrs((quicksum(x_new[i, j] for i in range(3)) >= 40*g_sub[j].x+bias[j] for j in range(3)))
    MP.optimize()
    
    S_y = np.array([y[i].x for i in range(len(y))])
    S_z = np.array([z[i].x for i in range(len(z))])
    LB = max(LB, MP.objval)
    # update the SP constrs according to the MP solution
    SP.remove(G1)
    SP.remove(G2)
    SP.addConstr(G@x>=h-M@g_sub-E@np.concatenate([S_y,S_z]),name="G1")
    SP.addConstr(G@x-h+E@np.concatenate([S_y,S_z])+M@g_sub<=bigM*(1-v),name='G2')
    SP.update()
    SP.optimize()
    # obtain the optimal y^{k+1}
    S_y = np.array([y[i].x for i in range(len(y))])
    S_z = np.array([z[i].x for i in range(len(z))])
    SP_obj = SP.ObjVal
    UB = min(UB, c@S_y+a@S_z+eta.x)
    k += 1
    # go back to the MP
    print("经过{}次迭代".format(k))
    print("上界为：{}".format(UB))
    print("下界为：{}".format(LB))

    
    







