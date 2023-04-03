# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 17:21:51 2023

@author: wyx
"""

from gurobipy import *
import numpy as np

model = Model('findMatrix')

bias = [206, 274, 220]

x_ij = model.addVars(3,3)
y_i = model.addVars(3,)
z_i = model.addVars(3,)
g_j = model.addVars(3,)

model.addConstrs((quicksum(x_ij[i,j] for j in range(3))*(-1)>=z_i[i]*(-1) for i in range(3)))
model.addConstrs((quicksum(x_ij[i,j] for i in range(3))>=40*g_j[j]+bias[j] for j in range(3)))

A = model.getA().toarray()  # Constraint matrix
b = np.array(model.RHS)  # Right-hand side as NumPy array
sense = np.array(model.sense)  # Sense array as NumPy array
Age = A[sense == '>', :]  # Submatrix corresponding to greater-or-equal constraints
bge = b[sense == '>']  # RHS for GE constraints

G = Age[:,:9]
E = Age[:,9:15]
M = Age[:,15:18]
h = bge