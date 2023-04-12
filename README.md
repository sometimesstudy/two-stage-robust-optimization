# Two-Stage-Robust-Optimization
Two stage robust optimization using a column-and-constraint generation (C&amp;CG) method.
This is an example code for the classical C&CG algorithm for two stage robust optimization, which is programmed by original Python Gurobi solver.

please refer to 

"Zeng B, Zhao L. Solving two-stage robust optimization problems using a column-and-constraint generation method[J]. Operations Research Letters, 2013, 41(5): 457-461."

for more details.

All the constrains were transformed into Matrix formulations.

C&CG iteration process

| interations   | LB  |UB|
|  ---- | ----  |----  |
| 1  | 14296 | 35238 |
| 2  | 33680 | 33680 |

Benders decomposition iteration process

| interations   | LB  |UB|
|  ---- | ----  |----  |
| 1  | 14296 | 35238 |
| 2  | 14665 | 35238 |
| 3  | 14860 | 35238 |
| 4  | 15227 | 35238 |
| 5  | 30532 | 34556 |
| 6  | 30956 | 34556 |
| 7  | 30958 | 34556 |
| 8  | 31383 | 34556 |
| 9  | 33127 | 33680 |
| 10  | 33543 | 33680 |
| 11  | 33680 | 33680 |
