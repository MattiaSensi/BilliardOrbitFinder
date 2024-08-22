# BilliardOrbitFinder
Matlab code to find Birkhoff and non-Birkhoff orbits in convex billiard via gradient flow.

These codes are provided as supplementary material of the preprint "Non-Birkhoff choreographies in symmetric billiards", https://arxiv.org/abs/XXXX.YYYYY, see Appendix B therein.

If you use these codes for your research, we kindly ask you to cite our preprint as well.

The gist of the algorithm is as follows:

Given a table of the form, for θ in [0,2π),

Γ(θ)=(Γ_X(θ),Γ_Y(θ))

and its derivative

Γ'(θ)=(Γ'_X(θ),Γ'_Y(θ)),

we substitute, for convenience, θ=2πt, to obtain a parametrization of the form, for t in [0,1),

Γ(t)=(Γ_X(2πt),Γ_Y(2πt))

and the corresponding derivative

Γ'(t)=(Γ'_X(2πt),Γ'_Y(2πt)).

Notice that these tables usually depend on one or more parameters, e.g. a and b are the usual semi-axes for an elliptic billiard.

We then provide an initial condition for the dynamics in the form of a sequence of N values t_i, i=1,2,...,N in [0,1).

These initial conditions can be provided through some symmetry (e.g., for N=4, we may provide t_1, t_2 and then give the relations t_3=-t_2, t_4=-t_1), and such symmetry can be enforced throughout the dynamics, by including them in the ODE solver (floww) in the Matlab code.

The examples included in this folder are the following:

For elliptic billiards, we include

1. billiards_NBO_manu.m, in which the initial conditions are provided as in formula (XX) of the manuscript

2. billiards_NBO_timerev.m, in which the initial conditions satisfy a time reversal symmetry, and form a candy-shaped orbit

3. billiards_NBO_updown.m, in which the initial conditions satisfy an up-down symmetry, and form a candy-shaped orbit

For D_n symmetric billiards, we include

4. billiard_NBO_D4.m, which was used to produce figure 2b) in the manuscript
