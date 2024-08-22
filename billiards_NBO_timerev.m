%If you use this code, kindly cite our preprint:
%https://arxiv.org/abs/AAAA.BBBBB

%a and b are the semi-axes of the ellipse
a=5;
b=4.5;

%parametrization of the ellipse; the table is parametrized for x \in [0,1)
Gamma1 = @(x,a,b) a*cos(2*pi*x);
Gamma2 = @(x,a,b) b*sin(2*pi*x);

%p and N to define a (Np,Nq) orbit, starting from an (Np,Nq) state.
%q is determined by the initial conditions.
p=2;
N=4;

%x0 contains the initial conditions for the gradient flow; we are working
%with the values x which will then be evaluated in Gamma1 and Gamma2,
%so these will be considered by the code to be mod 1
x0=zeros(N*p,1);

%here we impose a time reversal symmetry with a candy-shaped orbit
for ii=1:N*p/2
    x0(ii)=(-1)^(ii+1)*(2*ii-1)/(2*N*p);
end
for ii=N*p/2+1:N*p
    x0(ii)=-x0(N*p-ii+1);
end

%standard forward integration of the dynamical system with ode45
tspan = [0 1000];
[t,X]=ode45(@(t,X) floww(t, X, a, b, N*p), tspan, x0);

%in x1 we store the final values of the integration, which will be close
%to the equilibrium (i.e. the NBO) is tspan is large enough
x1=X(end,:);

%discretization to plot the ellipse
howmany=2000;

%lines 42-56 plot the orbit corresponding to the given initial condition in cyan
figure

hold on
ttime=linspace(0,1,howmany);

plot(Gamma1(ttime,a,b),Gamma2(ttime,a,b),'Color','k','LineWidth',1)
for jj=1:length(x0)-1
    line([Gamma1(x0(jj),a,b) Gamma1(x0(jj+1),a,b)],[Gamma2(x0(jj),a,b) Gamma2(x0(jj+1),a,b)],'Color','c','LineWidth',1.5);
end
line([Gamma1(x0(end),a,b) Gamma1(x0(1),a,b)],[Gamma2(x0(end),a,b) Gamma2(x0(1),a,b)],'Color','c','LineWidth',1.5)
axis equal
axis off
set(gcf,'color','w');

hold off

%lines 60-74 plot the orbit corresponding to the final orbit,
%an approximation of the desired NBO, in blue
figure

hold on
ttime=linspace(0,1,howmany);

plot(Gamma1(ttime,a,b),Gamma2(ttime,a,b),'Color','k','LineWidth',1)
for jj=1:length(x1)-1
    line([Gamma1(x1(jj),a,b) Gamma1(x1(jj+1),a,b)],[Gamma2(x1(jj),a,b) Gamma2(x1(jj+1),a,b)],'Color','b','LineWidth',1.5);
end
line([Gamma1(x1(end),a,b) Gamma1(x1(1),a,b)],[Gamma2(x1(end),a,b) Gamma2(x1(1),a,b)],'Color','b','LineWidth',1.5)
axis equal
axis off
set(gcf,'color','w');

hold off

%gradient flow
function [ dydt ] = floww (t, X, a, b, N)

dydt = zeros(N,1);

%We enforce the same symmetry as in the initial conditions
for ii=N/2+1:N
    X(ii)=-X(N-ii+1);
end

dydt(1)= HH(X(2),X(1),a,b)+HH(X(N),X(1),a,b);
for j=2:N-1
    dydt(j) = HH(X(j-1),X(j),a,b)+HH(X(j+1),X(j),a,b);
end
dydt(N)= HH(X(N-1),X(N),a,b)+HH(X(1),X(N),a,b);

end

function Acca = HH(A,B,a,b)
%the same parametrization given in lines 9 and 10
Gamma1 = @(x,a,b) a*cos(2*pi*x);
Gamma2 = @(x,a,b) b*sin(2*pi*x);

%derivative of Gamma1, Gamma2 with respect to x
G1 = @(x,a,b) 2*pi*(-a*sin(2*pi*x));
G2 = @(x,a,b) 2*pi*(b*cos(2*pi*x));

%function H as in the paper, equation (XX)
Acca=(1/(sqrt((Gamma1(A,a,b)-Gamma1(B,a,b))^2+(Gamma2(A,a,b)-Gamma2(B,a,b))^2)))*(G1(B,a,b)*(Gamma1(B,a,b)-Gamma1(A,a,b)) +G2(B,a,b)*(Gamma2(B,a,b)-Gamma2(A,a,b)));
end