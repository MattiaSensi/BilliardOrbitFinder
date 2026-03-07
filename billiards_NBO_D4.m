%If you use this code, kindly cite our paper:
%https://doi.org/10.1007/s00023-026-01670-7

%n corresponds to the D_n symmetry class; for convexity, one needs
%alpha<1/(1+n^2). Refer to Example 4.1, (7) and Appendix A
n=4;
alp=0.05;    %alpha

%parametrization of the D_n table; the table is parametrized for x \in [0,1)
Gamma1 = @(x,n,alp) (1+alp.*cos(2*n*pi*x)).*cos(2*pi*x);
Gamma2 = @(x,n,alp) (1+alp.*cos(2*n*pi*x)).*sin(2*pi*x);

%p to define a (p,q) orbit, starting from an (p,q) state.
%q is determined by the initial conditions.
p=12;

%x0 contains the initial conditions for the gradient flow; we are working
%with the values x which will then be evaluated in Gamma1 and Gamma2,
%so these will be considered by the code to be mod 1
x0=zeros(p,1);

x0(1)=1/2+0.05;
x0(6)=7/8;

%here we impose D_4 symmetry
x0(3)=x0(6)+1/4; x0(12)=x0(6)+1/2; x0(9)=x0(6)+3/4;
x0(10)=x0(1)+1/4; x0(7)=x0(1)+1/2; x0(4)=x0(1)+3/4;
x0(11)=1/4-x0(7); x0(2)=7/4-x0(10); x0(8)=3/4-x0(4); x0(5)=5/4-x0(1);

%standard forward integration of the dynamical system with ode45
tspan = [0 100];
opts = odeset('MaxStep',.01,'AbsTol',1e-09,'RelTol',1e-11);
[t,X]=ode45(@(t,X) floww(t, X, p, n, alp), tspan, x0, opts);

%in x1 we store the final values of the integration, which will be close
%to the equilibrium (i.e. the NBO) is tspan is large enough
x1=X(end,:);

%discretization to plot the D_n table
howmany=2000;

%lines 43-57 plot the initial condition in cyan
figure

hold on
ttime=linspace(0,1,howmany);

plot(Gamma1(ttime,n,alp),Gamma2(ttime,n,alp),'Color','k','LineWidth',1)
for jj=1:length(x0)-1
    line([Gamma1(x0(jj),n,alp) Gamma1(x0(jj+1),n,alp)],[Gamma2(x0(jj),n,alp) Gamma2(x0(jj+1),n,alp)],'Color','c','LineWidth',1.5);
end
line([Gamma1(x0(end),n,alp) Gamma1(x0(1),n,alp)],[Gamma2(x0(end),n,alp) Gamma2(x0(1),n,alp)],'Color','c','LineWidth',1.5)
axis equal
axis off
set(gcf,'color','w');

hold off

%lines 61-75 plot the orbit corresponding to the orbit,
%an approximation of the desired NBO, in blue
figure

hold on
ttime=linspace(0,1,howmany);
hold on
plot(Gamma1(ttime,n,alp),Gamma2(ttime,n,alp),'Color','k','LineWidth',1)
for jj=1:length(x1)-1  %plots final orbit
    line([Gamma1(x1(jj),n,alp) Gamma1(x1(jj+1),n,alp)],[Gamma2(x1(jj),n,alp) Gamma2(x1(jj+1),n,alp)],'Color','b','LineWidth',1.5);
end
line([Gamma1(x1(end),n,alp) Gamma1(x1(1),n,alp)],[Gamma2(x1(end),n,alp) Gamma2(x1(1),n,alp)],'Color','b','LineWidth',1.5)
axis equal
axis off
set(gcf,'color','w');

hold off

%gradient flow
function [ dydt ] = floww (t, X, N, n, alp)

dydt = zeros(N,1);

X(6)=7/8;

%We enforce the same symmetry as in the initial conditions
X(3)=X(6)+1/4; X(12)=X(6)+1/2; X(9)=X(6)+3/4;
X(10)=X(1)+1/4; X(7)=X(1)+1/2; X(4)=X(1)+3/4;
X(11)=1/4-X(7); X(2)=7/4-X(10); X(8)=3/4-X(4); X(5)=5/4-X(1);

dydt(1)= HH(X(2),X(1),n,alp)+HH(X(N),X(1),n,alp);
for j=2:N-1
    dydt(j) = HH(X(j-1),X(j),n,alp)+HH(X(j+1),X(j),n,alp);
end
dydt(N)= HH(X(N-1),X(N),n,alp)+HH(X(1),X(N),n,alp);

end

function Acca = HH(A,B,n,alp)
%the same parametrization given in lines 10 and 11
Gamma1 = @(x,n,alp) (1+alp.*cos(2*n*pi*x)).*cos(2*pi*x);
Gamma2 = @(x,n,alp) (1+alp.*cos(2*n*pi*x)).*sin(2*pi*x);

%derivative of Gamma1, Gamma2 with respect to x (a multiplicative 2*pi
%contribution was removed to reduce computational time)
G1 = @(x,n,alp) -n*alp*sin(2*n*pi*x)*cos(2*pi*x)-(1+alp*cos(2*n*pi*x))*sin(2*pi*x);
G2 = @(x,n,alp) -n*alp*sin(2*n*pi*x)*sin(2*pi*x)+(1+alp*cos(2*n*pi*x))*cos(2*pi*x);

%this function combines F^- and F^+ as in the paper, eq. (23) and (24)
Acca=(1/(sqrt((Gamma1(A,n,alp)-Gamma1(B,n,alp))^2+(Gamma2(A,n,alp)-Gamma2(B,n,alp))^2)))*(G1(B,n,alp)*(Gamma1(B,n,alp)-Gamma1(A,n,alp)) +G2(B,n,alp)*(Gamma2(B,n,alp)-Gamma2(A,n,alp)));

end
