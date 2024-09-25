# BilliardOrbitFinder
Matlab code to find Birkhoff and non-Birkhoff orbits in a convex billiard using gradient flow. For questions about the code, contact me at the email address you find here: https://mattiasensi.github.io/

These codes are provided as supplementary material of the preprint "Non-Birkhoff periodic orbits in symmetric billiards", https://arxiv.org/abs/XXXX.YYYYY, see Appendix B therein.

If you use these codes for your research, we kindly ask you to cite our preprint as well.

The gist of the algorithm is as follows:

Consider a billiard table of the form, for x in [0,1),

γ(x)=(γ_1(x),γ_2(x))

and the corresponding derivative

γ'(x)=(γ'_1(x),γ'_2(x)).

Notice that these tables may depend on one or more parameters, e.g. the lengths of the axes of an elliptic billiard. This is the case for all the examples we include here.

We then provide an initial condition for the dynamics in the form of a sequence of N values x_i, i=1,2,...,N in [0,1).

These initial conditions can be forced to possess some symmetry (e.g., for N=4, we may provide x_1, x_2 and then give the relations x_3=-x_2, x_4=-x_1), and such symmetry can be enforced throughout the gradient dynamics, by including them in the ODE solver (floww) in the Matlab code.

The ODE solver takes as input the initial conditions provided in the code, and integrates the system forward in time so that the values of x_i approximately reach an equilibrium in the prescribed symmetry space. By construction, such equilibrium is a proper orbit for the billiard, i.e. the angle of incidence and the angle of refraction on each point of the orbit are the same.

The examples included in this folder are the following:

For elliptic billiards, we include

1. billiards_NBO_manu.m, in which the initial conditions are provided as in formula (XX) of the manuscript

2. billiards_NBO_timerev.m, in which the initial conditions satisfy a time reversal symmetry, and therefore form a candy-shape

For a D4-symmetric Limacon-type billiards

3. billiards_NBO_D4.m, which was used to produce figure 1c) in the manuscript
