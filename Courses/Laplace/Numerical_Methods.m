%% Final Project
clear;clc;

%% Plot residual
omega = 1.0;
N=32;
[c_inner32,u32_init]=initial(N);
[u32, rtot32] = Gauss_Seidel(c_inner32,u32_init,N,omega);
N=96;
[c_inner96,u96_init]=initial(N);
[u96, rtot96] = Gauss_Seidel(c_inner96,u96_init,N,omega);
N=160;
[c_inner160,u160_init]=initial(N);
[u160, rtot160] = Gauss_Seidel(c_inner160,u160_init,N,omega);
N=224;
[c_inner224,u224_init]=initial(N);
[u224, rtot224] = Gauss_Seidel(c_inner224,u224_init,N,omega);
figure;
semilogy(rtot32, 'Color', [0.8, 0, 0], 'linewidth', 1);
hold on                                                     
semilogy(rtot96, 'Color', [0, 0, 0.8], 'linewidth', 1);
hold on
semilogy(rtot160, 'Color', [0, 0.8, 0], 'linewidth', 1);
hold on
semilogy(rtot224, 'Color', [0.7, 0, 0.7], 'linewidth', 1);
legend('Gauss-Seidel scheme (32\times32 grid)','Gauss-Seidel scheme (96\times96 grid)','Gauss-Seidel scheme (160\times160 grid)','Gauss-Seidel scheme (224\times224 grid)')
title('Total Residual versus interative index of', 'Gauss-Seidel scheme with different resolutions')
xlabel('Iterative Step')
ylabel('Total Residual')

%% Contour
% 32 grids 
L=1.0;
N=32;
x = linspace(0, L, N+1);
[X, Y] = meshgrid(x, x);
figure('Position',[100, 100, 600, 500]);
contourf(X,Y,u32,'ShowText','on',"LabelFormat","%0.1f")
colorbar;
xlabel('x')
ylabel('y')
title('Contour of the converged solution (32\times32 grid)')

% 96 grids
L=1.0;
N=96;
x = linspace(0, L, N+1);
[X, Y] = meshgrid(x, x);
figure('Position',[100, 100, 600, 500]);
contourf(X,Y,u96,10,'ShowText',true,"LabelFormat","%0.1f")
colorbar;
xlabel('x')
ylabel('y')
title('Contour of the converged solution (96\times96 grid)')

% 160 grids
L=1.0;
N=160;
x = linspace(0, L, N+1);
[X, Y] = meshgrid(x, x);
figure('Position',[100, 100, 600, 500]);
contourf(X,Y,u160,10,'ShowText',true,"LabelFormat","%0.1f")
colorbar;
xlabel('x')
ylabel('y')
title('Contour of the converged solution (160\times160 grid)')

% 224 grids
L=1.0;
N=224;
x = linspace(0, L, N+1);
[X, Y] = meshgrid(x, x);
figure('Position',[100, 100, 600, 500]);
contourf(X,Y,u224,10,'ShowText',true,"LabelFormat","%0.1f")
colorbar;
xlabel('x')
ylabel('y')
title('Contour of the converged solution (224\times224 grid)')

%% check correctness
[res32, err32, u_exact, u32]=find_error(32);
[res96, err96, u_exact, u96]=find_error(96);
[res160, err160, u_exact, u160]=find_error(160);
[res224, err224, u_exact, u224]=find_error(224);

%% Overrelaxation
N=32;
[c_inner32,u32_init]=initial(N);
omega = 1.0;
[u32_w10, rtot32_w10] = Gauss_Seidel(c_inner32,u32_init,N,omega);
omega = 1.1;
[u32_w11, rtot32_w11] = Gauss_Seidel(c_inner32,u32_init,N,omega);
omega = 1.3;
[u32_w13, rtot32_w13] = Gauss_Seidel(c_inner32,u32_init,N,omega);
omega = 1.5;
[u32_w15, rtot32_w15] = Gauss_Seidel(c_inner32,u32_init,N,omega);

figure;
semilogy(rtot32_w10,'Color',[0.8, 0, 0], 'linewidth', 1);
hold on
semilogy(rtot32_w11, 'Color',[0, 0.8, 0], 'linewidth', 1);
hold on
semilogy(rtot32_w13, 'Color',[0, 0, 0.8], 'linewidth', 1);
hold on
semilogy(rtot32_w15, 'Color',[0.7, 0, 0.7], 'linewidth', 1);
legend('\omega = 1.0', '\omega = 1.1','\omega = 1.3', '\omega = 1.5')
title({'Convergence procedure of point SOR Gauss-Seidel Scheme' ...
    ' with different overrelaxation parameters (32\times32 grids)'});
xlabel('Iterative Step')
ylabel('Total Residual')

%% Find optimal
omega = linspace(1.0,1.78,30);
iter=[];
for i = 1:length(omega)
    [u32, rtot32] = Gauss_Seidel(c_inner32,u32_init,N,omega(i));
    iter=[iter,length(rtot32)];
end
figure;
plot(omega,iter,'b-','linewidth',1)
title({'Iterative step needed for convergence with different relaxation parameters', ...
    'of Gauss-Seidel scheme (32\times32 grids)'});
xlabel('Relaxation Parameter')
ylabel('Iterative Step for Convergence')
[minvalue,i]=min(iter);
omega_opt=omega(i)

%% Calculate running time
tic
N=32;
omega=1;
[c_inner32,u32_init]=initial(N);
[u32, rtot32] = Gauss_Seidel(c_inner32,u32_init,N,omega);
elapsedTime32 = toc;

tic
N=96;
omega=1;
[c_inner96,u96_init]=initial(N);
[u96, rtot96] = Gauss_Seidel(c_inner96,u96_init,N,omega);
elapsedTime96 = toc;

tic
N=160;
omega=1;
[c_inner160,u160_init]=initial(N);
[u160, rtot160] = Gauss_Seidel(c_inner160,u160_init,N,omega);
elapsedTime160 = toc;

tic
N=224;
omega=1;
[c_inner224,u224_init]=initial(N);
[u224, rtot224] = Gauss_Seidel(c_inner224,u224_init,N,omega);
elapsedTime224 = toc;

grid_points=[32^2 96^2 160^2 224^2];
elapsedTime=[elapsedTime32 elapsedTime96 elapsedTime160 elapsedTime224];
logDx = log10(grid_points);
logDy = log10(elapsedTime);
p = polyfit(logDx, logDy, 1);
fitLogy = polyval(p, logDx);
figure(1)
loglog(grid_points, elapsedTime, 'bo'); 
hold on;
loglog(grid_points, 10.^fitLogy, 'b-'); 
hold off;
xlabel('Number of grid points');
ylabel('CPU time (s)');
title('CPU time variation of the solver with different numbers of grid points')
legend('CPU time', sprintf('CPU time fitting trend, slope = %.2f on log-log scale', p(1)))

%% Accuracy
omega = 1.0;
N=96;
[c_inner96,u96_init]=initial(N);
[u96, rtot96] = Gauss_Seidel(c_inner96,u96_init,N,omega);
N=120;
[c_inner120,u120_init]=initial(N);
[u120, rtot120] = Gauss_Seidel(c_inner120,u120_init,N,omega);
N=128;
[c_inner128,u128_init]=initial(N);
[u128, rtot128] = Gauss_Seidel(c_inner128,u128_init,N,omega);
N=140;
[c_inner140,u140_init]=initial(N);
[u140, rtot140] = Gauss_Seidel(c_inner140,u140_init,N,omega);
N=154;
[c_inner154,u154_init]=initial(N);
[u154, rtot154] = Gauss_Seidel(c_inner154,u154_init,N,omega);
N=160;
[c_inner160,u160_init]=initial(N);
[u160, rtot160] = Gauss_Seidel(c_inner160,u160_init,N,omega);
N=480;
[c_inner480,u480_init]=initial(N);
[u480, rtot480] = Gauss_Seidel(c_inner480,u480_init,N,omega);

[XFine, YFine] = meshgrid(linspace(0, 1, 481), linspace(0, 1, 481));

[X96, Y96] = meshgrid(linspace(0, 1, 97), linspace(0, 1, 97));
inter96 = interp2(XFine, YFine, u480, X96, Y96);
errorMatrix = abs(u96 - inter96);
rmsError96 = sqrt(mean(errorMatrix.^2));

[X120, Y120] = meshgrid(linspace(0, 1, 121), linspace(0, 1, 121));
inter120 = interp2(XFine, YFine, u480, X120, Y120);
errorMatrix = abs(u120 - inter120);
rmsError120 = sqrt(mean(errorMatrix.^2));

[X128, Y128] = meshgrid(linspace(0, 1, 129), linspace(0, 1, 129));
inter128 = interp2(XFine, YFine, u480, X128, Y128);
errorMatrix = abs(u128 - inter128);
rmsError128 = sqrt(mean(errorMatrix.^2));

[X140, Y140] = meshgrid(linspace(0, 1, 141), linspace(0, 1, 141));
inter140 = interp2(XFine, YFine, u480, X140, Y140);
errorMatrix = abs(u140 - inter140);
rmsError140 = sqrt(mean(errorMatrix.^2));

[X154, Y154] = meshgrid(linspace(0, 1, 155), linspace(0, 1, 155));
inter154 = interp2(XFine, YFine, u480, X154, Y154);
errorMatrix = abs(u154 - inter154);
rmsError154 = sqrt(mean(errorMatrix.^2));

[X160, Y160] = meshgrid(linspace(0, 1, 161), linspace(0, 1, 161));
inter160 = interp2(XFine, YFine, u480, X160, Y160);
errorMatrix = abs(u160 - inter160);
rmsError160 = sqrt(mean(errorMatrix.^2));

dx = [1/160 1/154 1/140 1/128 1/120 1/96];
er = [rmsError160 rmsError154 rmsError140 rmsError128 rmsError120 rmsError96];

logDx = log10(dx);
logEr = log10(er);
p1 = polyfit(logDx, logEr, 1);
fitLogEr = polyval(p1, logDx);
figure;
loglog(dx, er, 'bo');
hold on;
loglog(dx, 10.^fitLogEr, 'b-', 'linewidth',1.5);
slope = p1(1);  

xlabel('\Deltax');
ylabel('Error');
title('Error data with different grid sizes')
legend('Error', sprintf('Error fitting trend, slope = %.2f on log-log scale', slope))

%% Four-cylinder case
% Compare convergence
omega = 1.0;
N=32;
[c_inner32, u32_init] = initial_4c(N);
[u32_4c, rtot32_4c] = Gauss_Seidel(c_inner32,u32_init,N,omega);
[c_inner32,u32_init]=initial(N);
[u32, rtot32] = Gauss_Seidel(c_inner32,u32_init,N,omega);
figure;
semilogy(rtot32, 'Color', [0, 0, 0.8], 'linewidth', 1);
hold on
semilogy(rtot32_4c, '--', 'Color', [0, 0, 0.8], 'linewidth', 1);
legend('One circle','Four circles')
title('Convergence curves of two configurations', '(32\times32 grids)')
xlabel('Iterative Step')
ylabel('Total Residual')

%% Compare running time - Four-cylinder
omega=1;
tic
N=32;
[c_inner32,u32_init]=initial(N);
[u32, rtot32] = Gauss_Seidel(c_inner32,u32_init,N,omega);
elapsedTime32 = toc;

tic
N=96;
[c_inner96,u96_init]=initial(N);
[u96, rtot96] = Gauss_Seidel(c_inner96,u96_init,N,omega);
elapsedTime96 = toc;

tic
N=160;
[c_inner160,u160_init]=initial(N);
[u160, rtot160] = Gauss_Seidel(c_inner160,u160_init,N,omega);
elapsedTime160 = toc;

tic
N=224;
[c_inner224,u224_init]=initial(N);
[u224, rtot224] = Gauss_Seidel(c_inner224,u224_init,N,omega);
elapsedTime224 = toc;

grid_points=[32^2 96^2 160^2 224^2];
elapsedTime=[elapsedTime32 elapsedTime96 elapsedTime160 elapsedTime224];
figure;
loglog(grid_points, elapsedTime, 'bo-'); 

tic
N=32;
[c_inner32,u32_init]=initial_4c(N);
[u32, rtot32] = Gauss_Seidel(c_inner32,u32_init,N,omega);
elapsedTime32 = toc;

tic
N=96;
[c_inner96,u96_init]=initial_4c(N);
[u96, rtot96] = Gauss_Seidel(c_inner96,u96_init,N,omega);
elapsedTime96 = toc;

tic
N=160;
[c_inner160,u160_init]=initial_4c(N);
[u160, rtot160] = Gauss_Seidel(c_inner160,u160_init,N,omega);
elapsedTime160 = toc;

tic
N=224;
[c_inner224,u224_init]=initial_4c(N);
[u224, rtot224] = Gauss_Seidel(c_inner224,u224_init,N,omega);
elapsedTime224 = toc;

elapsedTime=[elapsedTime32 elapsedTime96 elapsedTime160 elapsedTime224];
hold on
loglog(grid_points, elapsedTime, 'ro-'); 

xlabel('Number of grid points');
ylabel('CPU time (s)');
title({'CPU time variation of the solver with different numbers of grid points', 'at two configurations'})
legend('CPU time-one circle', 'CPU time-four circles')

%% Overrelaxation - Four-cylinder
N=32;
[c_inner32,u32_init]=initial_4c(N);
omega = 1.0;
[u32_w10, rtot32_w10] = Gauss_Seidel(c_inner32,u32_init,N,omega);
omega = 1.1;
[u32_w11, rtot32_w11] = Gauss_Seidel(c_inner32,u32_init,N,omega);
omega = 1.3;
[u32_w13, rtot32_w13] = Gauss_Seidel(c_inner32,u32_init,N,omega);
omega = 1.5;
[u32_w15, rtot32_w15] = Gauss_Seidel(c_inner32,u32_init,N,omega);

figure;
semilogy(rtot32_w10,'Color',[0.8, 0, 0], 'linewidth', 1);
hold on
semilogy(rtot32_w11, 'Color',[0, 0.8, 0], 'linewidth', 1);
hold on
semilogy(rtot32_w13, 'Color',[0, 0, 0.8], 'linewidth', 1);
hold on
semilogy(rtot32_w15, 'Color',[0.7, 0, 0.7], 'linewidth', 1);
legend('\omega = 1.0', '\omega = 1.1','\omega = 1.3', '\omega = 1.5')
title({'Convergence procedure of point SOR Gauss-Seidel Scheme' ...
    ' with different overrelaxation parameters (32\times32 grids)' ...
    '(Four-cylinder case)'});
xlabel('Iterative Step')
ylabel('Total Residual')

%% Find optimal - Four-cylinder
omega = linspace(1.0,1.68,30);
N=32;
[c_inner32,u32_init]=initial_4c(N);
iter=[];
for i = 1:length(omega) 
    [u32, rtot32] = Gauss_Seidel(c_inner32,u32_init,N,omega(i));
    iter=[iter,length(rtot32)];
end
figure;
plot(omega,iter,'b-','linewidth',1)
title({'Iterative step needed for convergence with different relaxation parameters', ...
'of Gauss-Seidel scheme (32\times32 grids)' ...
'(Four-cylinder case)'}, 'fontweight', 'normal');
xlabel('Relaxation Parameter')
ylabel('Iterative Step for Convergence')
[minvalue,i]=min(iter);
omega_opt=omega(i)

%% Calculate running time - Four-cylinder
tic
N=32;
omega=1;
[c_inner32,u32_init]=initial(N);
[u32, rtot32] = Gauss_Seidel(c_inner32,u32_init,N,omega);
elapsedTime32 = toc;

tic
N=96;
omega=1;
[c_inner96,u96_init]=initial(N);
[u96, rtot96] = Gauss_Seidel(c_inner96,u96_init,N,omega);
elapsedTime96 = toc;

tic
N=160;
omega=1;
[c_inner160,u160_init]=initial(N);
[u160, rtot160] = Gauss_Seidel(c_inner160,u160_init,N,omega);
elapsedTime160 = toc;

tic
N=224;
omega=1;
[c_inner224,u224_init]=initial(N);
[u224, rtot224] = Gauss_Seidel(c_inner224,u224_init,N,omega);
elapsedTime224 = toc;

grid_points=[32^2 96^2 160^2 224^2];
elapsedTime=[elapsedTime32 elapsedTime96 elapsedTime160 elapsedTime224];
figure;
loglog(grid_points, elapsedTime, 'bo-'); 

tic
N=32;
[c_inner32,u32_init]=initial_4c(N);
[u32, rtot32] = Gauss_Seidel(c_inner32,u32_init,N,omega);
elapsedTime32 = toc;

tic
N=96;
[c_inner96,u96_init]=initial_4c(N);
[u96, rtot96] = Gauss_Seidel(c_inner96,u96_init,N,omega);
elapsedTime96 = toc;

tic
N=160;
[c_inner160,u160_init]=initial_4c(N);
[u160, rtot160] = Gauss_Seidel(c_inner160,u160_init,N,omega);
elapsedTime160 = toc;

tic
N=224;
[c_inner224,u224_init]=initial_4c(N);
[u224, rtot224] = Gauss_Seidel(c_inner224,u224_init,N,omega);
elapsedTime224 = toc;

elapsedTime=[elapsedTime32 elapsedTime96 elapsedTime160 elapsedTime224];
hold on
loglog(grid_points, elapsedTime, 'ro-'); 

xlabel('Number of grid points');
ylabel('CPU time (s)');
title({'CPU time variation of the solver with different numbers of grid points' ...
    'at two configurations'})
legend('CPU time-one circle', 'CPU time-four circles')

%% Accuracy - Four-cylinder case
omega = 1.0;
N=96;
[c_inner96,u96_init]=initial_4c(N);
[u96, rtot96] = Gauss_Seidel(c_inner96,u96_init,N,omega);
N=154;
[c_inner154,u154_init]=initial_4c(N);
[u154, rtot154] = Gauss_Seidel(c_inner154,u154_init,N,omega);
N=160;
[c_inner160,u160_init]=initial_4c(N);
[u160, rtot160] = Gauss_Seidel(c_inner160,u160_init,N,omega);
N=172;
[c_inner172,u172_init]=initial_4c(N);
[u172, rtot172] = Gauss_Seidel(c_inner172,u172_init,N,omega);
N=180;
[c_inner180,u180_init]=initial_4c(N);
[u180, rtot180] = Gauss_Seidel(c_inner180,u180_init,N,omega);
N=224;
[c_inner224,u224_init]=initial_4c(N);
[u224, rtot224] = Gauss_Seidel(c_inner224,u224_init,N,omega);
N=480;
[c_inner480,u480_init]=initial_4c(N);
[u480, rtot480] = Gauss_Seidel(c_inner480,u480_init,N,omega);

[XFine, YFine] = meshgrid(linspace(0, 1, 481), linspace(0, 1, 481));

[X96, Y96] = meshgrid(linspace(0, 1, 97), linspace(0, 1, 97));
inter96 = interp2(XFine, YFine, u480, X96, Y96);
errorMatrix = abs(u96 - inter96);
rmsError96 = sqrt(mean(errorMatrix.^2));

[X154, Y154] = meshgrid(linspace(0, 1, 155), linspace(0, 1, 155));
inter154 = interp2(XFine, YFine, u480, X154, Y154);
errorMatrix = abs(u154 - inter154);
rmsError154 = sqrt(mean(errorMatrix.^2));

[X160, Y160] = meshgrid(linspace(0, 1, 161), linspace(0, 1, 161));
inter160 = interp2(XFine, YFine, u480, X160, Y160);
errorMatrix = abs(u160 - inter160);
rmsError160 = sqrt(mean(errorMatrix.^2));

[X172, Y172] = meshgrid(linspace(0, 1, 173), linspace(0, 1, 173));
inter172 = interp2(XFine, YFine, u480, X172, Y172);
errorMatrix = abs(u172 - inter172);
rmsError172 = sqrt(mean(errorMatrix.^2));

[X180, Y180] = meshgrid(linspace(0, 1, 181), linspace(0, 1, 181));
inter180 = interp2(XFine, YFine, u480, X180, Y180);
errorMatrix = abs(u180 - inter180);
rmsError180 = sqrt(mean(errorMatrix.^2));

[X224, Y224] = meshgrid(linspace(0, 1, 225), linspace(0, 1, 225));
inter224 = interp2(XFine, YFine, u480, X224, Y224);
errorMatrix = abs(u224 - inter224);
rmsError224 = sqrt(mean(errorMatrix.^2));

dx = [1/224 1/180 1/172 1/160 1/154 1/96];
er = [rmsError224 rmsError180 rmsError172 rmsError160 rmsError154 rmsError96];

logDx = log10(dx);
logEr = log10(er);
p1 = polyfit(logDx, logEr, 1);
fitLogEr = polyval(p1, logDx);
figure;
loglog(dx, er, 'bo');
hold on;
loglog(dx, 10.^fitLogEr, 'b-', 'linewidth',1.5);
slope = p1(1);  

xlabel('\Deltax');
ylabel('Error');
title('Error data with different grid sizes (Four-cylinder case)')
legend('Error', sprintf('Error fitting trend, slope = %.2f on log-log scale', slope))

%% Line SOR
% Convergence
omega = 1.0;
N=32;
[u32sor, rtot32sor] = line_sor(N,omega);
[c_inner32,u32_init]=initial(N);
[u32, rtot32] = Gauss_Seidel(c_inner32,u32_init,N,omega);
figure;
semilogy(rtot32, 'Color', [0, 0, 0.8], 'linewidth', 1);
hold on
semilogy(rtot32sor, '--', 'Color', [0, 0, 0.8], 'linewidth', 1);
legend('Point SOR','Line SOR')
title('Convergence curve of point SOR and line SOR Gauss-Seidel schemes', '(32\times32 grids)')
xlabel('Iterative Step')
ylabel('Total Residual')

%% Contour - Line SOR
L=1.0;
N=32;
x = linspace(0, L, N+1);
[X, Y] = meshgrid(x, x);
figure('Position',[100, 100, 600, 500]);
contourf(X,Y,u32sor,[0.0 0.01 0.04 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.0],'ShowText','off',"LabelFormat","%0.1f")
colorbar;
xlabel('x')
ylabel('y')
title('Contour of the converged solution (32\times32 grid)')

%% Calculate running time - Line SOR
omega=1;
tic
N=32;
[u32, rtot32] = line_sor(N,omega);
elapsedTime32 = toc;

tic
N=96;
[u96, rtot96] = line_sor(N,omega);
elapsedTime96 = toc;

tic
N=160;
[u160, rtot160] = line_sor(N,omega);
elapsedTime160 = toc;

tic
N=224;
[u224, rtot224] = line_sor(N,omega);
elapsedTime224 = toc;

grid_points=[32^2 96^2 160^2 224^2];
elapsedTime=[elapsedTime32 elapsedTime96 elapsedTime160 elapsedTime224];
figure;
loglog(grid_points, elapsedTime, 'ro-'); 

tic
N=32;
[c_inner32,u32_init]=initial_4c(N);
[u32, rtot32] = Gauss_Seidel(c_inner32,u32_init,N,omega);
elapsedTime32 = toc;

tic
N=96;
[c_inner96,u96_init]=initial_4c(N);
[u96, rtot96] = Gauss_Seidel(c_inner96,u96_init,N,omega);
elapsedTime96 = toc;

tic
N=160;
[c_inner160,u160_init]=initial_4c(N);
[u160, rtot160] = Gauss_Seidel(c_inner160,u160_init,N,omega);
elapsedTime160 = toc;

tic
N=224;
[c_inner224,u224_init]=initial_4c(N);
[u224, rtot224] = Gauss_Seidel(c_inner224,u224_init,N,omega);
elapsedTime224 = toc;

elapsedTime=[elapsedTime32 elapsedTime96 elapsedTime160 elapsedTime224];
hold on
loglog(grid_points, elapsedTime, 'bo-'); 

xlabel('Number of grid points');
ylabel('CPU time (s)');
title({'CPU time variation of point SOR and line SOR Gauss-Seidel schemes', 'with different numbers of grid points'})
legend('CPU time-line', 'CPU time-point')

%% Streamlines
L = 1.0;
N = 224;
omega=1;
[c_inner96,u96_init]=initial_4c(N);
[u96, rtot96] = Gauss_Seidel_nmbc(c_inner96,u96_init,N,omega);

x = linspace(0, L, N+1);
y = linspace(0, L, N+1);
[X, Y] = meshgrid(x, y);

[v, u] = gradient(u96, x(2)-x(1), y(2)-y(1));

figure;
streamslice(X, Y, u, v, 3);
title('Streamslice of the Velocity Field');
xlabel('x');
ylabel('y');

%% Functions
function [u,rtota]=line_sor(N,omega)
[c_inner, u_init]=initial_4c(N);
u=u_init;
unew=u;
rtot=1;
rtota=[];
b = ones(N-1, 1);
dx=1/N;
dx2=dx^2;
while rtot>1e-5
rtot=0;
temp=u;

% row
for i = 2:N
    mainDiagx = (1 - c_inner(i, 2:N)) * (-4) + c_inner(i, 2:N);
    b(1:N-1) = -(1 - c_inner(i, 2:N)) .* (temp(i-1, 2:N) + temp(i+1, 2:N)) + c_inner(i, 2:N);
    upperDiag = 1 - c_inner(i, 2:N-1);
    lowerDiag = 1 - c_inner(i, 3:N);

    triDiagMatrix_row = diag(mainDiagx) + diag(upperDiag, 1) + diag(lowerDiag, -1);
    unew(i,2:N)=triDiagMatrix_row\b;
    temp(i,2:N)=(1-omega).*u(i,2:N)+omega.*unew(i,2:N);
    unew(i,2:N)=temp(i,2:N);
end


% column
for j = 2:N
    mainDiagx = (1 - c_inner(2:N, j)) * (-4.0) + c_inner(2:N, j);
    b(1:N-1) = -(1 - c_inner(2:N, j)) .* (temp(2:N, j-1) + temp(2:N, j+1)) + c_inner(2:N, j);
    upperDiag = 1 - c_inner(2:N-1, j);
    lowerDiag = 1 - c_inner(3:N, j);

    triDiagMatrix_row = diag(mainDiagx) + diag(upperDiag, 1) + diag(lowerDiag, -1);
    unew(2:N, j)=triDiagMatrix_row\b;
    temp(2:N, j)=(1-omega).*u(2:N,j)+omega.*unew(2:N,j);
    unew(2:N, j)=temp(2:N, j);
end

u=unew;
for i = 2:N
    for j = 2:N
        r(i,j)=abs( ...
           (1-c_inner(i,j))* ...
           ((u(i-1,j)-2*u(i,j)+u(i+1,j))/dx2+(u(i,j-1)-2*u(i,j)+u(i,j+1))/dx2));
        rtot = rtot + r(i,j);
    end
end
rtota=[rtota,rtot];
end

end

function [residual, error_relative, u_exact, u]=find_error(N)
R=sqrt(2)/2;
r1=0.25;
omega = 1.0;
[c_inner,u_init]=initial(N);
[u_initial]=change_boundary(u_init,N);
[u,rtot] = Gauss_Seidel(c_inner,u_initial,N,omega);

x = linspace(-0.5, 0.5, N+1);
[X, Y] = meshgrid(x, x);
dist_squared = sqrt(X.^2 + Y.^2);
u_exact = (log(dist_squared)-log(R))/(log(r1)-log(R));

u_exact(N/2+1,N/2+1)=1;
index_c_inner = find(c_inner==0);
num_0_c_inner = length(index_c_inner)-4*N;
residual = sum(abs((1-c_inner).*(u_exact-u)))/num_0_c_inner;
error_relative=sqrt(sum((1-c_inner(2:N,2:N)).*(u_exact(2:N,2:N)-u(2:N,2:N)).^2, 'all')/num_0_c_inner);

end

function [c_inner,u]=initial_4c(N)
L = 1.0;
r = 0.125;
cen1 = 0.25;
cen2 = 0.75;
x = linspace(0, L, N+1);
[X, Y] = meshgrid(x, x);

dist_squared1 = (X - cen1).^2 + (Y - cen2).^2;
dist_squared2 = (X - cen2).^2 + (Y - cen2).^2;
dist_squared3 = (X - cen1).^2 + (Y - cen1).^2;
dist_squared4 = (X - cen2).^2 + (Y - cen1).^2;
c_inner = zeros(N+1, N+1);
c_inner(dist_squared1 <= r^2) = 1;
c_inner(dist_squared2 <= r^2) = 1;
c_inner(dist_squared3 <= r^2) = 1;
c_inner(dist_squared4 <= r^2) = 1;

u = zeros(N+1, N+1);
randomMatrix = 2 * rand(N+1, N+1) - 1;
u(2:N,2:N) = (1 - c_inner(2:N,2:N)) .* randomMatrix(2:N,2:N) + c_inner(2:N,2:N) * 1.0;

end

function [u] = change_boundary(u_init,N)
dx=1/N;
R=sqrt(2)/2;
r1=0.25;
u=u_init;
for i = 1:N+1
    xe=-0.5+(i-1)*dx;
    r=sqrt(xe^2+0.5^2);
    u(1,i)=(log(r)-log(R))/(log(r1)-log(R));
    u(N+1,i)=(log(r)-log(R))/(log(r1)-log(R));
    u(i,1)=(log(r)-log(R))/(log(r1)-log(R));
    u(i,N+1)=(log(r)-log(R))/(log(r1)-log(R));
end

end

function [c_inner,u] = initial(N)
L = 1.0;
r = 0.25;
cen = 0.5;
x = linspace(0, L, N+1);
[X, Y] = meshgrid(x, x);

dist_squared = (X - cen).^2 + (Y - cen).^2;
c_inner = ones(N+1, N+1);
c_inner(dist_squared > r^2) = 0;

u = zeros(N+1, N+1);
randomMatrix = 2 * rand(N+1, N+1) - 1;
u(2:N,2:N) = (1 - c_inner(2:N,2:N)) .* randomMatrix(2:N,2:N) + c_inner(2:N,2:N) * 1.0;

end

function [u,rtota] = Gauss_Seidel(c_inner,u_init,N,omega)
u = u_init;
unew=u;
rtot=1;
rtota=[];
while rtot>1e-5
    rtot=0;
    temp=u;
    for i = 2:N
        for j = 2:N
            unew(i,j) = (1 - c_inner(i,j))*(1/4*(temp(i+1,j)+temp(i-1,j)+temp(i,j+1)+temp(i,j-1)))+c_inner(i,j) * 1.0;
            temp(i,j)=(1-omega)*u(i,j)+omega*unew(i,j);
            unew(i,j) = temp(i,j);
        end
    end
    u=unew;
    dx=1/N;
    for i = 2:N
        for j = 2:N
            r(i,j)=abs( ...
               (1-c_inner(i,j))* ...
               ((u(i-1,j)-2*u(i,j)+u(i+1,j))/dx^2+(u(i,j-1)-2*u(i,j)+u(i,j+1))/dx^2));
            rtot = rtot + r(i,j);
        end
    end
    rtota=[rtota,rtot];
end
end

function [u,rtota] = Gauss_Seidel_nmbc(c_inner,u_init,N,omega)
u = u_init;
unew=u;
rtot=1;
rtota=[];
L=1.0;
while rtot>1e-5
    rtot=0;
    temp=u;
    dx=L/N;
   
    for i = 2:N
        for j = 2:N
            unew(i,j) = (1 - c_inner(i,j))*(1/4*(temp(i+1,j)+temp(i-1,j)+temp(i,j+1)+temp(i,j-1)));
            temp(i,j)=(1-omega)*u(i,j)+omega*unew(i,j);
            unew(i,j) = temp(i,j);
        end
    end

    for j = 1:N+1
            u(1,j) = u(2,j)-dx
            u(N+1,j) = u(N,j)+dx;
    end
    for i = 1:N+1
            u(i,1) = u(i,2);
            u(i,N+1) = u(i,N); 
    end
    if c_inner(i,j)==1
        neighbors = [];
        
        neighbors(end+1) = u(i-1,j); 
        neighbors(end+1) = u(i+1,j); 
        neighbors(end+1) = u(i,j+1);
        neighbors(end+1) = u(i,j-1);
        
        u(i,j) = max(neighbors);
    end

    u=unew;
    
    for i = 2:N
        for j = 2:N
            r(i,j)=abs( ...
               (1-c_inner(i,j))* ...
               ((u(i-1,j)-2*u(i,j)+u(i+1,j))/dx^2+(u(i,j-1)-2*u(i,j)+u(i,j+1))/dx^2));
            rtot = rtot + r(i,j);
        end
    end
    rtota=[rtota,rtot];
end
end

