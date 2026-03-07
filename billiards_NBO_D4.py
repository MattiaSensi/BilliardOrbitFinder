# If you use this code, kindly cite our paper:
# https://doi.org/10.1007/s00023-026-01670-7

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# n corresponds to the D_n symmetry class; for convexity, one needs
# alpha<1/(1+n^2). Refer to Example 4.1, (7) and Appendix A

n = 4
alp = 0.05    # alpha

# parametrization of the D_n table; the table is parametrized for x \in [0,1)
def Gamma1(x, n, alp):
    return (1 + alp * np.cos(2 * n * np.pi * x)) * np.cos(2 * np.pi * x)

def Gamma2(x, n, alp):
    return (1 + alp * np.cos(2 * n * np.pi * x)) * np.sin(2 * np.pi * x)

# gradient flow
def HH(A, B, n, alp):
    # the same parametrization given in lines 16 and 17
    def Gamma1(x, n, alp):
        return (1 + alp * np.cos(2 * n * np.pi * x)) * np.cos(2 * np.pi * x)
    def Gamma2(x, n, alp):
        return (1 + alp * np.cos(2 * n * np.pi * x)) * np.sin(2 * np.pi * x)

    # derivative of Gamma1, Gamma2 with respect to x (a multiplicative 2*pi
    # contribution was removed to reduce computational time)
    def G1(x, n, alp):
        return -n * alp * np.sin(2 * n * np.pi * x) * np.cos(2 * np.pi * x) - (1 + alp * np.cos(2 * n * np.pi * x)) * np.sin(2 * np.pi * x)
    def G2(x, n, alp):
        return -n * alp * np.sin(2 * n * np.pi * x) * np.sin(2 * np.pi * x) + (1 + alp * np.cos(2 * n * np.pi * x)) * np.cos(2 * np.pi * x)

    # this function combines F^- and F^+ as in the paper, eq. (23) and (24)
    denom = np.sqrt((Gamma1(A, n, alp) - Gamma1(B, n, alp))**2 + (Gamma2(A, n, alp) - Gamma2(B, n, alp))**2)
    Acca = (1 / denom) * (G1(B, n, alp) * (Gamma1(B, n, alp) - Gamma1(A, n, alp)) +
                          G2(B, n, alp) * (Gamma2(B, n, alp) - Gamma2(A, n, alp)))
    return Acca

# p to define a (p,q) orbit, starting from an (p,q) state.
# q is determined by the initial conditions.
p = 12

# x0 contains the initial conditions for the gradient flow; we are working
# with the values x which will then be evaluated in Gamma1 and Gamma2,
# so these will be considered by the code to be mod 1
x0 = np.zeros(p)

x0[0] = 1/2 + 0.05
x0[5] = 7/8

# here we impose D_4 symmetry
x0[2]  = x0[5] + 1/4
x0[11] = x0[5] + 1/2
x0[8]  = x0[5] + 3/4
x0[9]  = x0[0] + 1/4
x0[6]  = x0[0] + 1/2
x0[3]  = x0[0] + 3/4
x0[10] = 1/4 - x0[6]
x0[1]  = 7/4 - x0[9]
x0[7]  = 3/4 - x0[3]
x0[4]  = 5/4 - x0[0]

# standard forward integration of the dynamical system
tspan = (0, 100)

def floww(t, X):
    dydt = np.zeros(p)
    X = X.copy()

    X[5] = 7/8

    # We enforce the same symmetry as in the initial conditions
    X[2]  = X[5] + 1/4
    X[11] = X[5] + 1/2
    X[8]  = X[5] + 3/4
    X[9]  = X[0] + 1/4
    X[6]  = X[0] + 1/2
    X[3]  = X[0] + 3/4
    X[10] = 1/4 - X[6]
    X[1]  = 7/4 - X[9]
    X[7]  = 3/4 - X[3]
    X[4]  = 5/4 - X[0]

    dydt[0] = HH(X[1], X[0], n, alp) + HH(X[-1], X[0], n, alp)
    for j in range(1, p-1):
        dydt[j] = HH(X[j-1], X[j], n, alp) + HH(X[j+1], X[j], n, alp)
    dydt[-1] = HH(X[-2], X[-1], n, alp) + HH(X[0], X[-1], n, alp)

    return dydt

# integrate using solve_ivp 
sol = solve_ivp(floww, tspan, x0, max_step=0.01, atol=1e-9, rtol=1e-11)
X = sol.y.T
x1 = X[-1, :]

# discretization to plot the D_n table
howmany = 2000

# lines 103-116 plot the initial condition in cyan
plt.figure()
plt.axis('equal')
plt.axis('off')
ttime = np.linspace(0, 1, howmany)

plt.plot(Gamma1(ttime, n, alp), Gamma2(ttime, n, alp), color='k', linewidth=1)
for jj in range(len(x0)-1):
    plt.plot([Gamma1(x0[jj], n, alp), Gamma1(x0[jj+1], n, alp)],
             [Gamma2(x0[jj], n, alp), Gamma2(x0[jj+1], n, alp)],
             color='c', linewidth=1.5)
plt.plot([Gamma1(x0[-1], n, alp), Gamma1(x0[0], n, alp)],
         [Gamma2(x0[-1], n, alp), Gamma2(x0[0], n, alp)],
         color='c', linewidth=1.5)
plt.gcf().set_facecolor('w')

# lines 120-132 plot the orbit corresponding to the orbit,
# an approximation of the desired NBO, in blue
plt.figure()
plt.axis('equal')
plt.axis('off')
ttime = np.linspace(0, 1, howmany)
plt.plot(Gamma1(ttime, n, alp), Gamma2(ttime, n, alp), color='k', linewidth=1)
for jj in range(len(x1)-1):  # plots final orbit
    plt.plot([Gamma1(x1[jj], n, alp), Gamma1(x1[jj+1], n, alp)],
             [Gamma2(x1[jj], n, alp), Gamma2(x1[jj+1], n, alp)],
             color='b', linewidth=1.5)
plt.plot([Gamma1(x1[-1], n, alp), Gamma1(x1[0], n, alp)],
         [Gamma2(x1[-1], n, alp), Gamma2(x1[0], n, alp)],
         color='b', linewidth=1.5)
plt.gcf().set_facecolor('w')

plt.show()

