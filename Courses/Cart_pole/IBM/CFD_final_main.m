%% Main code
clear; clc;
% Define parameters
dt = 0.01;
Re = 150;
dx = 1/32;
dy = dx;
Lx = 12;
Ly = 6;
D = 1;
bc = 1;
Nx = Lx/dx;
Ny = Ly/dy;
omega_u = 1;
omega_p = 1.5;

% ifluid
x_center = 4;
y_center = Ly/2;
r = D/2;
x = linspace(0, Lx, Nx+2);
y = linspace(0, Ly, Ny+2);
[X, Y] = meshgrid(x, y);
dist_squared = (X-x_center).^2+(Y-y_center).^2;
ifluid = zeros(Ny+2, Nx+2);
ifluid(dist_squared >= r^2) = 1;

% Initialization
un_m1 = zeros(Ny+2,Nx+2);
Un_m1 = zeros(Ny+2,Nx+2);
vn_m1 = zeros(Ny+2,Nx+2);
Vn_m1 = zeros(Ny+2,Nx+2);
if bc == 1
    [un_m1,Un_m1,vn_m1,Vn_m1] = BC(un_m1,Un_m1,vn_m1,Vn_m1,dy,Ly,Nx,Ny,ifluid);
else 
    [un_m1,Un_m1,vn_m1,Vn_m1] = BC_wall(un_m1,Un_m1,vn_m1,Vn_m1,dy,Ly,Nx,Ny,ifluid);
end
un = un_m1;
Un = Un_m1;
vn = vn_m1;
Vn = Vn_m1;
p_old = zeros(Ny+2,Nx+2);
p = zeros(Ny+2,Nx+2);

u_sample = zeros(Ny+2,Nx+2,200);
v_sample = zeros(Ny+2,Nx+2,200);
p_sample = zeros(Ny+2,Nx+2,200);

% Main loop
clc;
res = 1;
n = 0;
while n<10000 % 100s
    t_total=tic;
    n = n+1;
    % Step 1
    NL_u =  AB2(un_m1,un,Un_m1,Un,Vn_m1,Vn,Nx,Ny,dx,dy,dt);
    NL_v =  AB2(vn_m1,vn,Un_m1,Un,Vn_m1,Vn,Nx,Ny,dx,dy,dt);
    dfus_u = diffusion(un,Nx,Ny,dt,Re,dx,dy,ifluid);
    dfus_v = diffusion(vn,Nx,Ny,dt,Re,dx,dy,ifluid);
    tu=tic;
    [u_star,U_star,iter_u] = line_sor(NL_u,dfus_u,un,dt,Re,dx,dy,Nx,Ny,1,omega_u,Ly,ifluid);
    time_u_star = toc(tu);
    tv=tic;
    [v_star,V_star,iter_v] = line_sor(NL_v,dfus_v,vn,dt,Re,dx,dy,Nx,Ny,2,omega_u,Ly,ifluid);
    time_v_star = toc(tv);
    if bc == 1
        [u_star,U_star,v_star,V_star] = BC(u_star,U_star,v_star,V_star,dy,Ly,Nx,Ny,ifluid);
    else
        [u_star,U_star,v_star,V_star] = BC_wall(u_star,U_star,v_star,V_star,dy,Ly,Nx,Ny,ifluid);
    end
    
    % Step 2
    tp=tic;
    [p,iter_p] = Possion(U_star,V_star,dx,dy,dt,Nx,Ny,p_old,omega_p,ifluid);
    time_p = toc(tp);

    % Step 3
    [un_p1,Un_p1] = fractional_step(u_star,U_star,p,dt,dx,dy,Nx,Ny,1,ifluid);
    [vn_p1,Vn_p1] = fractional_step(v_star,V_star,p,dt,dx,dy,Nx,Ny,2,ifluid);
    if bc == 1
        [un_p1,Un_p1,vn_p1,Vn_p1] = BC(un_p1,Un_p1,vn_p1,Vn_p1,dy,Ly,Nx,Ny,ifluid);
    else
        [un_p1,Un_p1,vn_p1,Vn_p1] = BC_wall(un_p1,Un_p1,vn_p1,Vn_p1,dy,Ly,Nx,Ny,ifluid);
    end
    umax = max(un_p1, [], 'all');
    vmax = max(vn_p1, [], 'all');
    cfl_u = umax*dt/dx^2;
    cfl_v = vmax*dt/dy^2;
    iunstable = 0;
    if cfl_u > 1 || cfl_v > 1
        iunstable = 1;
    end
    
    % Step 4
    tres=tic;
    [max_residual,sum_residual] = residual(un_p1,un,vn_p1,vn,p,p_old);
    res = sum_residual;
    time_res = toc(tres);

    % Update u,v
    un_m1 = un;
    vn_m1 = vn;
    Un_m1 = Un;
    Vn_m1 = Vn;

    un = un_p1;
    vn = vn_p1;
    Un = Un_p1;
    Vn = Vn_p1;
    p_old = p;
    time_total = toc(t_total);

    if rem(n, 1) == 0
        fprintf('n=%d, T=%f, res=%f, iter_u=%d, iter_v=%d, iter_p=%d, iunstable=%d, cfl_u=%f, cfl_v=%f, max_residual=%f, sum_residual=%f\n', ...
         n, n*dt, res, iter_u, iter_v, iter_p, iunstable, cfl_u, cfl_v, max_residual, sum_residual);
        fprintf('T_u=%.2f%%, T_v=%.2f%%, T_p=%.2f%%, T_res=%.2f%%\n', ...
         time_u_star/time_total * 100, time_v_star/time_total  * 100, time_p/time_total  * 100, time_res/time_total  * 100);
    end

    rsampling = 50;
    if rem(n, rsampling) == 0
        u_sample(:,:,n/rsampling) = un;
        v_sample(:,:,n/rsampling) = vn;
        p_sample(:,:,n/rsampling) = p;
    end
    
    if mod(n,2000) == 0
        filename = sprintf('Data_%d.mat', n/2000);
        save(filename);
    end

end

%% Step 1 u_star
function [u_star,U_star,k] = line_sor(NL,dfus,un,dt,Re,dx,dy,Nx,Ny,uv,omega,Ly,ifluid)
    src = NL+dfus;
    R = dt/(2*Re);
    rtot = 1;
    uk = un;
    k = 0;
    cw = ones(Ny+2,Nx+2)*(-R/dx^2);
    ce = ones(Ny+2,Nx+2)*(-R/dx^2);
    cn = ones(Ny+2,Nx+2)*(-R/dy^2);
    cs = ones(Ny+2,Nx+2)*(-R/dy^2);
    cp = ones(Ny+2,Nx+2)*(1+2*R/dx^2+2*R/dy^2);
    [cp,ce,cw,cn,cs,src] = coeff_u(cp,ce,cw,cn,cs,src,ifluid,Nx,Ny);
    while rtot > 1e-6
        k = k+1;
        temp = uk;
        for i = 2:Ny+1
            B = src(i,2:Nx+1)-cs(i,2:Nx+1).*temp(i-1,2:Nx+1)-cn(i,2:Nx+1).*temp(i+1,2:Nx+1);
            temp(i,2:Nx+1) = tdmac(cw(i,3:Nx+1),cp(i,2:Nx+1),ce(i,2:Nx),B,Nx);
        end
        temp(2:Ny+1,2:Nx+1)=(1-omega).*uk(2:Ny+1,2:Nx+1)+omega.*temp(2:Ny+1,2:Nx+1);
        if uv == 1
            [temp,~,~,~] = BC(temp,temp,temp,temp,dy,Ly,Nx,Ny,ifluid);
        else
            [~,~,temp,~] = BC(temp,temp,temp,temp,dy,Ly,Nx,Ny,ifluid);
        end
        uk_1 = temp;
        udiff = uk_1(2:end-1,2:end-1) - uk(2:end-1,2:end-1);
        rtot = rms(udiff,'all');
        uk = uk_1;
    end
    u_star = uk_1;
    U_star = zeros(Ny+2,Nx+2);
    if uv == 1
        U_star(2:Ny+1,2:Nx+2) = (u_star(2:Ny+1,2:Nx+2)+u_star(2:Ny+1,1:Nx+1))/2;
    else
        U_star(2:Ny+2,2:Nx+1) = (u_star(2:Ny+2,2:Nx+1)+u_star(1:Ny+1,2:Nx+1))/2;
    end

end

function NL = AB2(un_m1,un,Un_m1,Un,Vn_m1,Vn,Nx,Ny,dx,dy,dt)
NL = zeros(Ny+2,Nx+2);
    for i = 2:Ny+1
        NL(i,2:Nx+1) = -1/2*dt*( ...
            3*((Un(i,3:Nx+2).*(un(i,3:Nx+2)+un(i,2:Nx+1))/2 - Un(i,2:Nx+1).*(un(i,2:Nx+1)+un(i,1:Nx))/2)/dx ...
            + (Vn(i+1,2:Nx+1).*(un(i+1,2:Nx+1)+un(i,2:Nx+1))/2 - Vn(i,2:Nx+1).*(un(i,2:Nx+1)+un(i-1,2:Nx+1))/2)/dy) ...
            -((Un_m1(i,3:Nx+2).*(un_m1(i,3:Nx+2)+un_m1(i,2:Nx+1))/2 - Un_m1(i,2:Nx+1).*(un_m1(i,2:Nx+1)+un_m1(i,1:Nx))/2)/dx ...
            + (Vn_m1(i+1,2:Nx+1).*(un_m1(i+1,2:Nx+1)+un_m1(i,2:Nx+1))/2 - Vn_m1(i,2:Nx+1).*(un_m1(i,2:Nx+1)+un_m1(i-1,2:Nx+1))/2)/dy) ...
            );
    end
end

function dfus = diffusion(un,Nx,Ny,dt,Re,dx,dy,ifluid)
    dfus = zeros(Ny+2,Nx+2);
    src = zeros(Ny+2,Nx+2);
    R = dt/2/Re;
    cw = ones(Ny+2,Nx+2)*R/dx^2;
    ce = ones(Ny+2,Nx+2)*R/dx^2;
    cn = ones(Ny+2,Nx+2)*R/dy^2;
    cs = ones(Ny+2,Nx+2)*R/dy^2;
    cp = ones(Ny+2,Nx+2)*(1+R*(-2/dx^2-2/dy^2));
    [cp,ce,cw,cn,cs,src]=coeff_u(cp,ce,cw,cn,cs,src,ifluid,Nx,Ny);
    dfus(2:Ny+1,2:Nx+1) = cp(2:Ny+1,2:Nx+1).*un(2:Ny+1,2:Nx+1)+ce(2:Ny+1,2:Nx+1).*un(2:Ny+1,3:Nx+2) ...
        +cw(2:Ny+1,2:Nx+1).*un(2:Ny+1,1:Nx)+cn(2:Ny+1,2:Nx+1).*un(3:Ny+2,2:Nx+1) ...
        +cs(2:Ny+1,2:Nx+1).*un(1:Ny,2:Nx+1);
end

%% Step 2 p
function [p,k] = Possion(U_star,V_star,dx,dy,dt,Nx,Ny,p_old,omega,ifluid)
    pk = p_old;
    rtot = 1;
    U_star(2:Ny-1,end)=U_star(2:Ny-1,end)+mean(U_star(2:Ny-1,2)-U_star(2:Ny-1,end));
    src = zeros(Ny+2,Nx+2);
    src(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).* ...
        (((U_star(2:Ny+1,3:Nx+2)-U_star(2:Ny+1,2:Nx+1))/dx+(V_star(3:Ny+2,2:Nx+1)-V_star(2:Ny+1,2:Nx+1))/dy)/dt);
    k = 0;
    cw = ones(Ny+2,Nx+2)*(1/dx^2);
    ce = ones(Ny+2,Nx+2)*(1/dx^2);
    cn = ones(Ny+2,Nx+2)*(1/dy^2);
    cs = ones(Ny+2,Nx+2)*(1/dy^2);
    cp = ones(Ny+2,Nx+2)*(-2/dx^2-2/dy^2);
    [cp,ce,cw,cn,cs,src] = coeff_p(cp,ce,cw,cn,cs,src,ifluid,Nx,Ny);
    while rtot > 1e-6
        k = k+1;
        temp = pk;
        for i = 2:Ny+1
            B = src(i,2:Nx+1)-cs(i,2:Nx+1).*temp(i-1,2:Nx+1)-cn(i,2:Nx+1).*temp(i+1,2:Nx+1);
            temp(i,2:Nx+1) = tdmac(cw(i,3:Nx+1),cp(i,2:Nx+1),ce(i,2:Nx),B,Nx);
        end
        temp(2:Ny+1,2:Nx+1)=(1-omega).*pk(2:Ny+1,2:Nx+1)+omega.*temp(2:Ny+1,2:Nx+1);
       
        temp(1,2:Nx+1) = temp(2,2:Nx+1);
        temp(end,2:Nx+1) = temp(end-1,2:Nx+1);
        temp(2:Ny+1,1) = temp(2:Ny+1,2);
        temp(2:Ny+1,end) = temp(2:Ny+1,end-1);
        pk_1 = temp;
        
        pdiff = pk_1(2:end-1,2:end-1) - pk(2:end-1,2:end-1);
        rtot = rms(pdiff,'all');
        pk = pk_1;
        if rem(k, 100) == 0
            fprintf('k=%d, residual=%f\n', ...
            k, rtot);
        end
    end
    p = pk_1;

end

%% step 3 un_1
function [un_1,Un_1] = fractional_step(u_star,U_star,p,dt,dx,dy,Nx,Ny,uv,ifluid)
    un_1 = zeros(Ny+2,Nx+2);
    Un_1 = zeros(Ny+2,Nx+2);
    pcc = zeros(Ny+2,Nx+2);
    pfc = zeros(Ny+2,Nx+2);
    src = zeros(Ny+2,Nx+2);
    if uv == 1
        cw = -ones(Ny+2,Nx+2);
        ce = ones(Ny+2,Nx+2);
        cn = zeros(Ny+2,Nx+2);
        cs = zeros(Ny+2,Nx+2);
        cp = zeros(Ny+2,Nx+2);
        [cp,ce,cw,cn,cs,src] = coeff_p(cp,ce,cw,cn,cs,src,ifluid,Nx,Ny);
        pcc(2:Ny+1,2:Nx+1) = (ce(2:Ny+1,2:Nx+1).*p(2:Ny+1,3:Nx+2)+cp(2:Ny+1,2:Nx+1).*p(2:Ny+1,2:Nx+1)+cw(2:Ny+1,2:Nx+1).*p(2:Ny+1,1:Nx))/2/dx;

        cw = -ones(Ny+2,Nx+2);
        ce = zeros(Ny+2,Nx+2);
        cn = zeros(Ny+2,Nx+2);
        cs = zeros(Ny+2,Nx+2);
        cp = ones(Ny+2,Nx+2);
        [cp,ce,cw,cn,cs,src] = coeff_p(cp,ce,cw,cn,cs,src,ifluid,Nx,Ny);
        pfc(2:Ny+1,2:Nx+2) = (cp(2:Ny+1,2:Nx+2).*p(2:Ny+1,2:Nx+2)+cw(2:Ny+1,2:Nx+2).*p(2:Ny+1,1:Nx+1))/dx;

        un_1(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(u_star(2:Ny+1,2:Nx+1)-dt*pcc(2:Ny+1,2:Nx+1));
        Un_1(2:Ny+1,2:Nx+2) = U_star(2:Ny+1,2:Nx+2)-dt*pfc(2:Ny+1,2:Nx+2);
    else
        cw = zeros(Ny+2,Nx+2);
        ce = zeros(Ny+2,Nx+2);
        cn = ones(Ny+2,Nx+2);
        cs = -ones(Ny+2,Nx+2);
        cp = zeros(Ny+2,Nx+2);
        [cp,ce,cw,cn,cs,src] = coeff_p(cp,ce,cw,cn,cs,src,ifluid,Nx,Ny);
        pcc(2:Ny+1,2:Nx+1) = (cn(2:Ny+1,2:Nx+1).*p(3:Ny+2,2:Nx+1)+cp(2:Ny+1,2:Nx+1).*p(2:Ny+1,2:Nx+1)+cs(2:Ny+1,2:Nx+1).*p(1:Ny,2:Nx+1))/2/dy;

        cw = zeros(Ny+2,Nx+2);
        ce = zeros(Ny+2,Nx+2);
        cn = zeros(Ny+2,Nx+2);
        cs = -ones(Ny+2,Nx+2);
        cp = ones(Ny+2,Nx+2);
        [cp,ce,cw,cn,cs,src] = coeff_p(cp,ce,cw,cn,cs,src,ifluid,Nx,Ny);
        pfc(2:Ny+2,2:Nx+1) = (cp(2:Ny+2,2:Nx+1).*p(2:Ny+2,2:Nx+1)+cs(2:Ny+2,2:Nx+1).*p(1:Ny+1,2:Nx+1))/dy;

        un_1(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(u_star(2:Ny+1,2:Nx+1)-dt*pcc(2:Ny+1,2:Nx+1));
        Un_1(2:Ny+2,2:Nx+1) = U_star(2:Ny+2,2:Nx+1)-dt*pfc(2:Ny+2,2:Nx+1);
    end    
end

%% Step 4 Residual
function [max_residual,sum_residual] = residual(un_p1,un,vn_p1,vn,p,p_old)
    dff_u = un_p1(2:end-1,2:end-1)-un(2:end-1,2:end-1);
    dff_v = vn_p1(2:end-1,2:end-1)-vn(2:end-1,2:end-1);
    dff_p = p(2:end-1,2:end-1)-p_old(2:end-1,2:end-1);
    max_residual = max(abs(dff_u) + abs(dff_v) + abs(dff_p), [], 'all');
    sum_residual = rms(dff_u + dff_v+dff_p,'all');
end

%% Boundary condition
function [un,Un,vn,Vn] = BC(un_old,Un_old,vn_old,Vn_old,dy,Ly,Nx,Ny,ifluid)
  % left
    Un_old(:,2) = 1;
    un_old(:,1) = 2-un_old(:,2);
    vn_old(:,1) = -vn_old(:,2);
  % right
    Un_old(:,end) = Un_old(:,end-1);
    un_old(:,end) = un_old(:,end-1);
    vn_old(:,end) = vn_old(:,end-1);
  % upper
    un_old(end,:) = un_old(end-1,:);
    vn_old(end,:) = -vn_old(end-1,:);
    Vn_old(end,:) = 0;
  % bottom
    un_old(1,:) = un_old(2,:);
    vn_old(1,:) = -vn_old(2,:);
    Vn_old(2,:) = 0;
  % inner
    Un_old(2:Ny+1,2:Nx+2) = ifluid(2:Ny+1,2:Nx+2).*Un_old(2:Ny+1,2:Nx+2).*ifluid(2:Ny+1,1:Nx+1);
    Vn_old(2:Ny+2,2:Nx+1) = ifluid(2:Ny+2,2:Nx+1).*Vn_old(2:Ny+2,2:Nx+1).*ifluid(1:Ny+1,2:Nx+1);
  % update
    un = un_old;
    vn = vn_old;
    Un = Un_old;
    Vn = Vn_old;
end

function [un,Un,vn,Vn] = BC_wall(un_old,Un_old,vn_old,Vn_old,dy,Ly,Nx,Ny,ifluid)
  % left
    Un_old(:,2) = 1;
    un_old(:,1) = 2-un_old(:,2);
    vn_old(:,1) = -vn_old(:,2);
  % right
    Un_old(:,end) = un_old(:,end-1);
    un_old(:,end) = un_old(:,end-1);
    vn_old(:,end) = vn_old(:,end-1);
  % upper
    un_old(end,:) = -un_old(end-1,:);
    vn_old(end,:) = -vn_old(end-1,:);
    Vn_old(end,:) = 0;
  % bottom
    un_old(1,:) = -un_old(2,:);
    vn_old(1,:) = -vn_old(2,:);
    Vn_old(2,:) = 0;
  % inner
    Un_old(2:Ny+1,2:Nx+2) = ifluid(2:Ny+1,2:Nx+2).*Un_old(2:Ny+1,2:Nx+2).*ifluid(2:Ny+1,1:Nx+1);
    Vn_old(2:Ny+2,2:Nx+1) = ifluid(2:Ny+2,2:Nx+1).*Vn_old(2:Ny+2,2:Nx+1).*ifluid(1:Ny+1,2:Nx+1);
  % update
    un = un_old;
    vn = vn_old;
    Un = Un_old;
    Vn = Vn_old;
end

function [cp,ce,cw,cn,cs,src] = coeff_u(cp,ce,cw,cn,cs,src,ifluid,Nx,Ny)
    cp(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(cp(2:Ny+1,2:Nx+1)- ...
        (1-ifluid(2:Ny+1,1:Nx)).*cw(2:Ny+1,2:Nx+1)- ...
        (1-ifluid(2:Ny+1,3:Nx+2)).*ce(2:Ny+1,2:Nx+1)- ...
        (1-ifluid(3:Ny+2,2:Nx+1)).*cn(2:Ny+1,2:Nx+1)- ...
        (1-ifluid(1:Ny,2:Nx+1)).*cs(2:Ny+1,2:Nx+1))+ ...
        (1-ifluid(2:Ny+1,2:Nx+1));
    cw(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(ifluid(2:Ny+1,1:Nx).*cw(2:Ny+1,2:Nx+1));
    ce(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(ifluid(2:Ny+1,3:Nx+2).*ce(2:Ny+1,2:Nx+1));
    cw(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(ifluid(3:Ny+2,2:Nx+1).*cn(2:Ny+1,2:Nx+1));
    cs(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(ifluid(1:Ny,2:Nx+1).*cs(2:Ny+1,2:Nx+1));
    src(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*src(2:Ny+1,2:Nx+1);
end

function [cp,ce,cw,cn,cs,src] = coeff_p(cp,ce,cw,cn,cs,src,ifluid,Nx,Ny)
    cp(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(cp(2:Ny+1,2:Nx+1)+ ...
        (1-ifluid(2:Ny+1,1:Nx)).*cw(2:Ny+1,2:Nx+1)+ ...
        (1-ifluid(2:Ny+1,3:Nx+2)).*ce(2:Ny+1,2:Nx+1)+ ...
        (1-ifluid(3:Ny+2,2:Nx+1)).*cn(2:Ny+1,2:Nx+1)+ ...
        (1-ifluid(1:Ny,2:Nx+1)).*cs(2:Ny+1,2:Nx+1)) + ...
        (1-ifluid(2:Ny+1,2:Nx+1));
    cw(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(ifluid(2:Ny+1,1:Nx).*cw(2:Ny+1,2:Nx+1));
    ce(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(ifluid(2:Ny+1,3:Nx+2).*ce(2:Ny+1,2:Nx+1));
    cn(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(ifluid(3:Ny+2,2:Nx+1).*cn(2:Ny+1,2:Nx+1));
    cs(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*(ifluid(1:Ny,2:Nx+1).*cs(2:Ny+1,2:Nx+1));
    src(2:Ny+1,2:Nx+1) = ifluid(2:Ny+1,2:Nx+1).*src(2:Ny+1,2:Nx+1);
end

function [x] = tdmac(a,b,c,d,n)
    for i=2:n
        b(i) = b(i) - c(i-1)/b(i-1)*a(i-1);
        d(i) = d(i) - d(i-1)/b(i-1)*a(i-1);
    end
    x(n) = d(n)/b(n);
    for i=n-1:-1:1
        x(i) = (d(i) - c(i)*x(i+1))/b(i);
    end
end
