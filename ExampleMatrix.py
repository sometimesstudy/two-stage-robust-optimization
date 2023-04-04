# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 17:21:51 2023

@author: wyx
"""

from gurobipy import *
import numpy as np

model = Model('findMatrix')

bias = [206, 274, 220]

y_i = model.addVars(3,)
z_i = model.addVars(3,)
x_ij = model.addVars(3,3)
g_j = model.addVars(3,)

model.addConstrs((quicksum(x_ij[i,j] for j in range(3))*(-1)>=z_i[i]*(-1) for i in range(3)))
model.addConstrs((quicksum(x_ij[i,j] for i in range(3))>=40*g_j[j]+bias[j] for j in range(3)))
model.addConstrs((-z_i[i] >= -800*y_i[i] for i in range(3)))
model.addConstr((quicksum(z_i[i] for i in range(3)) >= 772))

A = model.getA().toarray()  # Constraint matrix
b = np.array(model.RHS)  # Right-hand side as NumPy array
sense = np.array(model.sense)  # Sense array as NumPy array
Age = A[sense == '>', :]  # Submatrix corresponding to greater-or-equal constraints
bge = b[sense == '>']  # RHS for GE constraints

Ay, by, G, E, M, D, h = [], [], [], [], [], [], []
s1 = len(y_i) + len(z_i) # the lenghth of first stage variables
# find the first stage constrains Ay>=d
for row in np.arange(Age.shape[0]):
    if np.all(Age[row][s1:] == 0): # extract the fisrst stage matrix
        Ay.append(Age[row][0:s1])
        by.append(bge[row])
    else: # extract the first-second stage nexus variables  
        E.append(Age[row][:s1])
        G.append(Age[row][s1:s1+len(x_ij)])
        M.append(Age[row][s1+len(x_ij):])
        h.append(bge[row])
Ay, by, G, E, M, h = np.array(Ay), np.array(by), np.array(G), np.array(E), np.array(M), np.array(h)
