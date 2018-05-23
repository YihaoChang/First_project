close all;
clear all;
clc

%% The parameters of the simulation
% part1:
rho = 1;
Re = 40;
D = 0.2;
% Re = rho*V*D/mu = V*D/nu = Q/A*D/nu  , let V =1 and Re= D/nu
nu = D/Re;
% nu = mu/ rho
mu = rho * nu;
% Part2:
c = 1e+07;
epsilon = 1/(c * mu);
% part3:
ball_center=[0 ,0];

%% Create the mesh
domain = [-13.4 * D, 16.5 * D, -8.35 * D, 8.35 * D];
subdomain = [-D, D, -D, D];

% non-uniform mesh 250*160
nx1 = 85;
nx2 = 60;% mesh the subdomain
nx3 = 105;

ny1 = 50;
ny2 = 60;% mesh the subdomain
ny3 = 50;

nx = nx1 + nx2 + nx3;
ny = ny1 + ny2 + ny3;

x1 = linspace(domain(1), subdomain(1), nx1 + 1);
x2 = linspace(subdomain(1), subdomain(2), nx2 + 1);
x3 = linspace(subdomain(2), domain(2), nx3 + 1);
x = unique([x1, x2, x3])';
y1 = linspace(domain(3), subdomain(3), ny1 + 1);
y2 = linspace(subdomain(3), subdomain(4), ny2+1);
y3 = linspace(subdomain(4), domain(4), ny3+1);
y = unique([y1, y2, y3])';
xc = (x(2 : end) + x(1 : end - 1)) / 2;
yc = (y(2 : end) + y(1 : end - 1)) / 2;

[Xu, Yu]= ndgrid(x(2:end-1),yc);
[Xv, Yv]= ndgrid(xc,y(2:end-1));
[Xc, Yc] = ndgrid(xc, yc);

hx = diff(x);
hy = diff(y);
hxc = diff(xc);
hyc = diff(yc);
subhx = (subdomain(2)-subdomain(1)) / nx2;
subhy = (subdomain(4)-subdomain(3)) / ny2;

%% Time discretize parameter
CFL = 0.4; % original 0.37
ft = 50;
dt = max(subhx,subhy) * CFL;
total_step_number = ceil(ft/dt);
nt = total_step_number;
%% Discretize the interface
resolution = 101;
nb = linspace(0,2*pi,resolution);
bx = ball_center(1) + cos(nb) * (1/2) * D;
by = ball_center(2) + sin(nb) * (1/2) * D;
cylinder = [bx',by'];

eta1=inpolygon(Xu,Yu,bx,by);
eta2=inpolygon(Xv,Yv,bx,by);

%% Construct the matrix by using finite difference with staggered mesh
% for Lap velocity and time step
Lap_u = get_Laplcain_U(nx, ny, hx, hyc);
Lap_v = get_Laplcain_V(nx, ny, hxc, hy);

Iu = speye((nx-1)*ny);
Iv = speye(nx*(ny-1));

A1 = (1/dt) * Iu - nu * Lap_u;
A2 = (1/dt) * Iv - nu * Lap_v;
A = blkdiag(A1,A2);
% for  gradient -p
[Bx, By] = getGp(nx,ny,hxc,hyc);
B = [Bx;By];

% for div velocity
[Cx, Cy] = getDu(nx,ny,hx,hy);
C = [Cx,Cy];

%% Initial value
u = zeros(nx-1, ny);
v = zeros(nx, ny-1);
p = zeros(nx, ny);

Virtual_force_1 = zeros(nx-1, ny);
Virtual_force_2 = zeros(nx, ny-1); 

drag = zeros(nt+1,1);
lift = zeros(nt+1,1);

%% Start simulation
fprintf('progress bar\n--20%%--40%%--60%%--80%%--100%%\n')
timestep = 0;
for nt = 1:nt
    timestep = timestep + 1;
    t = timestep * dt;
    % g:=0  
    g = zeros(nx,ny);
    % r_u = r_v :=0
    r_u = zeros(nx-1,ny);
    r_v = zeros(nx,ny-1); 
     %% Boundary condition
    uW = yc*0 + 1;
    uE = u(end,:);
    
    uall = [uW';u;uE];
    
    uN = uall(:,end);
    uS = uall(:,1);
    
    uall = [(8/3)*uS-2*uall(:,1)+(1/3)*uall(:,2), uall, (8/3)*uN-2*uall(:,end)+(1/3)*uall(:,end-1)];
    
    vN = xc*0;
    vS = xc*0;
    vall = [vS v vN];
    
    vW = y*0;
    vE = vall(end,:);
    vall = [(8/3)*vW'-2*vall(1,:)+(1/3)*vall(2,:);vall;(8/3)*vE-2*vall(end,:)+(1/3)*vall(end-1,:)];
    
    % G = g-bc
    G = g;
    % G(:,1)   = G(:,1) + vS./hy(1);
    % G(:,end) = G(:,end) - vN./hy(end);
    G(1,:)   = G(1,:) + uW'./hx(1);
    %G(end,:) = G(end,:) - uE./hx(end);
    
    % R = r -bc
    R_u = r_u;
    R_v = r_v;
    % R_u(:,1)   = R_u(:,1) + 2*uS(2:end-1)./hy(1)^2;
    % R_u(:,end) = R_u(:,end) + 2*uN(2:end-1)./hy(end).^2;
      R_u(1,:)   = R_u(1,:) + uW'./hx(1).^2;
    % R_u(end,:) = R_u(end,:) + uE./hx(end).^2;
    
    % R_v(:,1)   = R_v(:,1) + vS./hy(1).^2;
    % R_v(:,end) = R_v(:,end) + vN./hy(end).^2;
    % R_v(1,:)   = R_v(1,:) + 2*vW(2:end-1)'./hx(1).^2;
    % R_v(end,:) = R_v(end,:) + 2*vE(2:end-1)./hx(end).^2;
    %% Nonlinear
    [nonlinear_u, nonlinear_v] = treat_nonlinesr(u,v,uall,vall,nx,ny,hx,hy,hxc,hyc);
    
    %% Solve the problem
% operater: A = d/dt(u) - nu * Lap(u), B = [-d/dx;-d/dy], C = [d/dx,d/dy]
% Recall: G = g - bc, R =r-bc
% p = p + (dt/epsilon) * (G - C * u)
% R - nonlinear + F = Au + grad(p) = Au - Bp = Au - B(p+(dt/epsilon)*(G-C*u))
%                   = Au - Bp - (dt/epsilon)*B*G + (dt/epsilon)*B*C*u
% Au + (dt/epsilon)*B*C*u = R - nonlinear + B*p + (dt/epsilon)*B*G +F
%-----------------------------------------------------------------
    % Stage1: calculate the virtual force
    lhs = A + (dt / epsilon) * (B * C);
    rhs = (1 / dt) * [u(:); v(:)] + [R_u(:);R_v(:)]-[nonlinear_u(:); nonlinear_v(:)]...
      + B * p(:) + (dt/epsilon)*B*G(:) + [Virtual_force_1(:);Virtual_force_2(:)] ;
    usol = lhs \ rhs;
    

    u_1 = reshape(usol(1:(nx-1)*ny), nx-1, ny);
    v_1 = reshape(usol((nx-1)*ny+1:end), nx, ny-1);
    
    pre_Virtual_force_1 = eta1.*(-u_1)/dt;
    pre_Virtual_force_2 = eta2.*(-v_1)/dt;
    
    % Stage2: calculate the velocity field, pressure and save the virtual force
    lhs = A + (dt / epsilon) * (B * C);
    rhs = (1 / dt) * [u(:); v(:)] + [R_u(:);R_v(:)]-[nonlinear_u(:); nonlinear_v(:)]...
      + B * p(:) + (dt/epsilon)*B*G(:) + [pre_Virtual_force_1(:);pre_Virtual_force_2(:)] ;
    usol = lhs \ rhs;
    

    u_2 = reshape(usol(1:(nx-1)*ny), nx-1, ny);
    v_2 = reshape(usol((nx-1)*ny+1:end), nx, ny-1);
    
    cor_Virtual_force_1 = eta1.*(-u_2)/dt;
    cor_Virtual_force_2 = eta2.*(-v_2)/dt;  
   
    AA_Virtual_force_1 = pre_Virtual_force_1 + cor_Virtual_force_1;
    AA_Virtual_force_2 = pre_Virtual_force_2 + cor_Virtual_force_2;
    
    u = u_2 + dt * cor_Virtual_force_1;
    v = v_2 + dt * cor_Virtual_force_2;
    
    p = p(:) + (dt / epsilon) * (G(:) - C * usol);
    p = reshape(p, nx, ny);
    
    %% Restructure the velocity with boundary 
    U = zeros(nx+1,ny);
    U(1,:) = 1;
    U(2:end-1,:) = u;
    U(end,:) = U(end-1,:);
   
    V = zeros(nx,ny+1);
    V(:,2:end-1) = v;
    
%% Drag  and lift
drag(nt+1) = -(sum(sum(AA_Virtual_force_1)) * min(hx)^2) / (D/2);
lift(nt+1) = -(sum(sum(AA_Virtual_force_2)) * min(hy)^2) / (D/2);   
    %% Visualize

 [X,Y] = ndgrid(x(2:end-1),y(2:end-1));
  if ~mod(nt, 20)
% Compute the vortucity
 omega = diff(V(:,2:end-1))./ hxc(:,ones(ny-1,1)) - (diff(U(2:end-1,:)')./ hyc(:,ones(nx-1,1)))';

 fig1=figure(1);

 positive_level = setdiff(0:15,0);
 nagative_level = setdiff(-15:0,0);

 contour(X,Y,omega,nagative_level); hold on
 contour(X,Y,omega,positive_level);
% fill(cylinder(:,1),cylinder(:,2),'k','linewidth',1.0);  % particle
 fill(cylinder(:,1),cylinder(:,2),'k','linewidth',1.0);  hold off;


 set(gca,'fontsize',16);
% title(['t=', num2str(t)]);
 title(sprintf('t = %6.4f',nt*dt));
 axis equal
 axis([-0.5,3,-1,1]);
 box on
 drawnow
         %saveas(fig1,sprintf('%05d.eps',k),'psc2');
         saveas(fig1,sprintf('%05d.png',nt),'png');
 
 
 
 %
    if floor(25*nt/total_step_number)>floor(25*(nt-1)/total_step_number)
        fprintf('*')
    end
  end 
end
% fprintf(' ... done\n')
plot_drag_lift(dt,total_step_number,drag,lift)
























