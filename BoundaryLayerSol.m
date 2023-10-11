% Assignment of Constant Parameters
U = 10;
dx = 0.0005;
dy = 0.0005;
Re = 10000;
delta = 5/sqrt(Re);
L=1;
H = 2*delta;
nu = U*L/Re;
nx =L/dx;
ny=H/dy;

% Velocities
u = zeros(nx+1,ny+1);
v = zeros(nx+1,ny+1);

% Boundary Conditions
u(1,:) = 1;
u(:,end) = 1;
u(:,1)=0;
v(1,:) = 0;
u(1,1)=0;

% Iterative Calculation of Velocities
for i=1:nx
 for j = 2:ny
 u(i+1,j)=u(i,j)+(dx/(dy*dy*Re))*(u(i,j+1)-2*u(i,j)+u(i,j-1))-(0.5*dx/dy)*(v(i,j)/u(i,j))*(u(i,j+1)-u(i,j-1));
 v(i+1,j) = v(i+1,j-1) -(dy/(2*dx))*(u(i+1,j)-u(i,j)+u(i+1,j-1)-u(i,j-1));
 end
end


% Rows at Given points
x1 = 0.0005;
x2 = 0.5;
x_point1=(x1/dx)+1;
x_point2=(x2/dx)+1;

% u at Given points
u_point1=u(x_point1,:);
u_point2=u(x_point2,:);

% v at Given points
v_point1=v(x_point1,:);
v_point2=v(x_point2,:);

% Grid Spacing
x=linspace(0,L,ny+1);
y=linspace(0,H,ny+1);

% At x=0.0005
etaB_x1=y.*sqrt(U/(nu*x1));
[eta_x1,f_x1]=blasiusSol(etaB_x1);
y_etaBx1=(eta_x1./sqrt(U/(nu*x1)));
u_Bx1=f_x1(:,2);
v_Bx1=(0.5.*sqrt(nu/(x1*U)).*(eta_x1.*f_x1(:,2)-f_x1(:,1)));

% At x=0.5
etaB_x2=y.*sqrt(U/(nu*x2));
[eta_x2,f_x2]=blasiusSol(etaB_x2);
y_etaBx2=(eta_x2./sqrt(U/(nu*x2)));
u_Bx2=f_x2(:,2);
v_Bx2=(0.5.*sqrt(nu/(x2*U)).*(eta_x2.*f_x2(:,2)-f_x2(:,1)));

% Transpose of Velocities
u_T=u.';
v_T=v.';

figure(1)
contourf(u_T,30);
colormap("hot");
colorbar();
title('U');
xlabel('x/L');
ylabel('y/delta');

figure(2)
contourf(v_T,15);
colormap("turbo");
colorbar()
title('V')
xlabel('x/L');
ylabel('y/delta');

figure(3)
plot(y,u_point1,'-gx','DisplayName','Computational');
hold on
title('u/U at X = 0.0005');
ylabel('u/U');
xlabel('y/delta');
plot(y_etaBx1,u_Bx1,'-bx','DisplayName','Blasius Solution');
hold off
legend('Computational','Blasius Solution');

figure(4)
plot(y,v_point1,'-gx','DisplayName','Computational');
hold on
title('v/U at X = 0.0005');
ylabel('v/U');
xlabel('y/delta');
plot(y_etaBx1,v_Bx1,'-bx','DisplayName','Blasius Solution');
hold off
legend('Computational','Blasius Solution');

figure(5)
plot(y,u_point2,'-gx','DisplayName','Computational');
hold on
title('u/U at X = 0.5');
ylabel('u/U');
xlabel('y/delta');
plot(y_etaBx2,u_Bx2,'-bx','DisplayName','Blasius Solution');
hold off
legend('Computational','Blasius Solution');

figure(6)
plot(y,v_point2,'-gx','DisplayName','Computational');
hold on
title('v/U at X = 0.5');
ylabel('v/U');
xlabel('y/delta');
plot(y_etaBx2,v_Bx2,'-bx','DisplayName','Blasius Solution');
hold off
legend('Computational','Blasius Solution');

function[eta,f]=blasiusSol(eta_B)
[eta,f]=ode45(@D,eta_B,[0 0 0.332]);
end

function de_L = D(eta,f)
de_L=zeros(3,1);
de_L(1)=f(2);
de_L(2)=f(3);
de_L(3)=-(1/2)*f(1)*f(3);
end