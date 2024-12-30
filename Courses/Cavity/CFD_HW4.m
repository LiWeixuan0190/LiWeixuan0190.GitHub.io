clear; clc;
%% Define parameters
N = 16;
dx = 1/N;
dt = 0.001;
Re = 100;
omega_u = 1;
omega_p = 1;

%% Initialization
un_m1 = zeros(N+2,N+2);
Un_m1 = zeros(N+2,N+2);
vn_m1 = zeros(N+2,N+2);
Vn_m1 = zeros(N+2,N+2);
[un_m1,Un_m1,vn_m1,Vn_m1] = BC(un_m1,Un_m1,vn_m1,Vn_m1);

un = un_m1;
Un = Un_m1;
vn = vn_m1;
Vn = Vn_m1;
p_old = zeros(N+2,N+2);
p = zeros(N+2,N+2);

%% Main loop
clc;
res = 1;
n = 0;
while res>1e-6
    t_total=tic;
    n = n+1;
    % Step 1
    NL_u =  AB2(un_m1,un,Un_m1,Un,Vn_m1,Vn,N,dx,dt);
    NL_v =  AB2(vn_m1,vn,Un_m1,Un,Vn_m1,Vn,N,dx,dt);
    dfus_u = diffusion(un,N,dt,Re,dx);
    dfus_v = diffusion(vn,N,dt,Re,dx);
    tu=tic;
    [u_star,U_star,iter_u] = line_sor(NL_u,dfus_u,un,dt,Re,dx,N,1,omega_u);
    time_u_star = toc(tu);
    tv=tic;
    [v_star,V_star,iter_v] = line_sor(NL_v,dfus_v,vn,dt,Re,dx,N,2,omega_u);
    time_v_star = toc(tv);
    [u_star,U_star,v_star,V_star] = BC(u_star,U_star,v_star,V_star);
    
    % Step 2
    tp=tic;
    [p,iter_p] = Possion(U_star,V_star,dx,dt,N,p_old,omega_p);
    time_p = toc(tp);

    % Step 3
    [un_p1,Un_p1] = fractional_step(u_star,U_star,p,dt,dx,N,1);
    [vn_p1,Vn_p1] = fractional_step(v_star,V_star,p,dt,dx,N,2);
    [un_p1,Un_p1,vn_p1,Vn_p1] = BC(un_p1,Un_p1,vn_p1,Vn_p1);
    umax = max(un_p1, [], 'all');
    vmax = max(vn_p1, [], 'all');
    cfl_u = umax*dt/dx^2;
    cfl_v = vmax*dt/dx^2;
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

end

%% Write tecplot file
[X, Y] = meshgrid(linspace(0, 1, N), linspace(0, 1, N));

fid = fopen('filename.dat', 'wt');

fprintf(fid, 'TITLE = "Velocity Field"\n');
fprintf(fid, 'VARIABLES = "X" "Y" "U" "V"\n');
fprintf(fid, 'ZONE T="Zone 1", I=%d, J=%d, F=POINT\n', N, N);

for j = 2:N+1
    for i = 2:N+1
        fprintf(fid, '%f %f %f %f\n', X(i-1,j-1), Y(i-1,j-1), un(i,j), vn(i,j));
    end
end
fclose(fid);

%% Step 1 u_star
function [u_star,U_star,k] = line_sor(NL,dfus,un,dt,Re,dx,N,uv,omega)
    b = NL+dfus;
    R = dt/(2*Re*dx^2);
    rtot = 1;
    uk = un;
    k = 0;
    while rtot > 1e-6
        k = k+1;
        temp = uk;
        mainDiagx = (1+4*R)*ones(N,1);
        upperDiag = -R*ones(N-1, 1);
        lowerDiag = -R*ones(N-1, 1);
        triDiagMatrix = diag(mainDiagx) + diag(upperDiag, 1) + diag(lowerDiag, -1);
        for i = 2:N+1
            B = b(i,2:N+1)+R*(temp(i-1,2:N+1)+temp(i+1,2:N+1));
            temp(i,2:N+1) = triDiagMatrix\B';
        end
        temp(2:N+1,2:N+1)=(1-omega).*uk(2:N+1,2:N+1)+omega.*temp(2:N+1,2:N+1);
        if uv == 1
            [temp,~,~,~] = BC(temp,temp,temp,temp);
        else
            [~,~,temp,~] = BC(temp,temp,temp,temp);
        end
        uk_1 = temp;
        rtot = sum(sum(abs(uk_1-uk)));
        uk = uk_1;
    end
    u_star = uk_1;

    if uv == 1
        U_star(2:N+1,2:N+2) = (u_star(2:N+1,2:N+2)+u_star(2:N+1,1:N+1))/2;
    else
        U_star(2:N+2,2:N+1) = (u_star(2:N+2,2:N+1)+u_star(1:N+1,2:N+1))/2;
    end

end

function [un,Un,vn,Vn] = BC(un_old,Un_old,vn_old,Vn_old)
    Un_old(:,2) = 0;
    Un_old(:,end) = 0;
    Un = Un_old;

    Vn_old(2,:) = 0;
    Vn_old(end,:) = 0;
    Vn = Vn_old;

    un_old(end,:) = 2-un_old(end-1,:);
    un_old(1,:) = -un_old(2,:);
    un_old(:,1) = -un_old(:,2);
    un_old(:,end) = -un_old(:,end-1);
    un = un_old;

    vn_old(end,:) = -vn_old(end-1,:);
    vn_old(1,:) = -vn_old(2,:);
    vn_old(:,1) = -vn_old(:,2);
    vn_old(:,end) = -vn_old(:,end-1);
    vn = vn_old;

end

function NL = AB2(un_m1,un,Un_m1,Un,Vn_m1,Vn,N,dx,dt)
dy = dx;
NL = zeros(N+2,N+2);
    for i = 2:N+1
        NL(i,2:N+1) = -1/2*dt*( ...
            3*((Un(i,3:N+2).*(un(i,3:N+2)+un(i,2:N+1))/2 - Un(i,2:N+1).*(un(i,2:N+1)+un(i,1:N))/2)/dx ...
            + (Vn(i+1,2:N+1).*(un(i+1,2:N+1)+un(i,2:N+1))/2 - Vn(i,2:N+1).*(un(i,2:N+1)+un(i-1,2:N+1))/2)/dy ...
            ) ...
            -((Un_m1(i,3:N+2).*(un_m1(i,3:N+2)+un_m1(i,2:N+1))/2 - Un_m1(i,2:N+1).*(un_m1(i,2:N+1)+un_m1(i,1:N))/2)/dx ...
            + (Vn_m1(i+1,2:N+1).*(un_m1(i+1,2:N+1)+un_m1(i,2:N+1))/2 - Vn_m1(i,2:N+1).*(un_m1(i,2:N+1)+un_m1(i-1,2:N+1))/2)/dy) ...
            );
    end
end

function dfus = diffusion(un,N,dt,Re,dx)
    dy = dx;
    dfus = zeros(N+2,N+2);
    for i = 2:N+1
        dfus(i,2:N+1) = un(i,2:N+1)+dt/2/Re*((un(i,3:N+2)-2*un(i,2:N+1)+un(i,1:N))/dx^2 ...
            + (un(i+1,2:N+1)-2*un(i,2:N+1)+un(i-1,2:N+1))/dy^2);
    end
end

%% Step 2 p
function [p,k] = Possion(U_star,V_star,dx,dt,N,p_old,omega)
    dy = dx;
    pk = p_old;
    rtot = 1;
    b = zeros(N+2,N+2);
    b(2:N+1,2:N+1) = ((U_star(2:N+1,3:N+2)-U_star(2:N+1,2:N+1))/dx + (V_star(3:N+2,2:N+1)-V_star(2:N+1,2:N+1))/dy)/dt;
    k = 0;
    while rtot > 1e-6
        k = k+1;
        temp = pk;
        for i = 2:N+1
            temp(i,2:N+1) = -(b(i,2:N+1)*dx^2-(temp(i-1,2:N+1)+temp(i+1,2:N+1)+temp(i,1:N)+temp(i,3:N+2)))/4;
        end
        temp(2:N+1,2:N+1)=(1-omega).*pk(2:N+1,2:N+1)+omega.*temp(2:N+1,2:N+1);
        temp(1,2:N+1) = temp(2,2:N+1);
        temp(N+2,2:N+1) = temp(N+1,2:N+1);
        temp(2:N+1,1) = temp(2:N+1,2);
        temp(2:N+1,N+2) = temp(2:N+1,N+1);
        pk_1 = temp;
        res=zeros(N+2,N+2);
        for i = 2:N+1
            res(i,2:N+1)=pk_1(i-1,2:N+1)+pk_1(i+1,2:N+1)+pk_1(i,1:N)+pk_1(i,3:N+2)-4*pk_1(i,2:N+1)-b(i,2:N+1)*dx^2;
        end
        rtot = sum(sum(abs(res(2:N+1,2:N+1))));
        pk = pk_1;
    end
    p = pk_1;

end

%% step 3 un_1
function [un_1,Un_1] = fractional_step(u_star,U_star,p,dt,dx,N,uv)
    un_1 = zeros(N+2,N+2);
    Un_1 = zeros(N+2,N+2);
    pcc = zeros(N+2,N+2);
    pfc = zeros(N+2,N+2);
    if uv == 1
        pcc(2:N+1,2:N+1) = (p(2:N+1,3:N+2)-p(2:N+1,1:N))/2/dx;
        pfc(2:N+1,2:N+2) = (p(2:N+1,2:N+2)-p(2:N+1,1:N+1))/dx;
        un_1(2:N+1,2:N+1) = u_star(2:N+1,2:N+1)-dt*pcc(2:N+1,2:N+1);
        Un_1(2:N+1,2:N+2) = U_star(2:N+1,2:N+2)-dt*pfc(2:N+1,2:N+2);
    else
        pcc(2:N+1,2:N+1) = (p(3:N+2,2:N+1)-p(1:N,2:N+1))/2/dx;
        pfc(2:N+2,2:N+1) = (p(2:N+2,2:N+1)-p(1:N+1,2:N+1))/dx;
        un_1(2:N+1,2:N+1) = u_star(2:N+1,2:N+1)-dt*pcc(2:N+1,2:N+1);
        Un_1(2:N+2,2:N+1) = U_star(2:N+2,2:N+1)-dt*pfc(2:N+2,2:N+1);
    end    
end

%% Step 4 Residual
function [max_residual,sum_residual] = residual(un_p1,un,vn_p1,vn,p,p_old)
    dff_u = abs(un_p1-un);
    dff_v = abs(vn_p1-vn);
    dff_p = abs(p-p_old);
    max_residual = max(dff_u + dff_v+dff_p, [], 'all');
    sum_residual = sum(sum(dff_u +dff_v+dff_p));
    
end
