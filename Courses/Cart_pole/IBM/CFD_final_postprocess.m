%% Post-processing code
% Contours
dt_sample = 0.5;
i = 100/dt_sample; % 100s
directory = '';

% Pressure contour
[X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));
figure
contourf(X, Y, p_sample(2:end-1,2:end-1,i),50,'LineStyle','none')
colorbar('southoutside');
load('colorData.mat','BuDRd_18')
colormap(othercolor('BuDRd_18'));
shading flat;
setPivot(0)
title(['Pressure contour at ', sprintf('Re=%d ',Re), '(resolution: ', sprintf('%.d',Nx), '\times', sprintf('%d)',Ny)])
axis equal
xlim([0 Lx])
ylim([0 Ly])
filename = sprintf('p_re%d_%ds.png', Re, i*0.5);
filepath = fullfile(directory, filename);
saveas(gcf, filepath);

% u Velocity contour
[X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));
figure('Position', [100, 100, 600, 400]);
scale_factor = 2;
velocity_magnitude = sqrt(un(2:Ny+1,2:Nx+1).^2 +vn(2:Ny+1,2:Nx+1).^2);
contourf(X, Y, u_sample(2:Ny+1,2:Nx+1,i),50, 'LineColor', 'none');
colorbar('southoutside');
shading flat;
axis equal
title(['u velocity contour at ', sprintf('Re=%d ',Re), '(resolution: ', sprintf('%.d',Nx), '\times', sprintf('%d)',Ny)])
xlim([0 Lx])
ylim([0 Ly])
filename = sprintf('un_re%d_%ds.png', Re, i*0.5);
filepath = fullfile(directory, filename);
saveas(gcf, filepath);

% v Velocity contour
[X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));
figure('Position', [100, 100, 600, 400]);
scale_factor = 2;
velocity_magnitude = sqrt(un(2:Ny+1,2:Nx+1).^2 +vn(2:Ny+1,2:Nx+1).^2);
contourf(X, Y, v_sample(2:Ny+1,2:Nx+1,i),50, 'LineColor', 'none');
colorbar('southoutside');
shading flat;
axis equal
title(['v velocity contour at ', sprintf('Re=%d ',Re), '(resolution: ', sprintf('%.d',Nx), '\times', sprintf('%d)',Ny)])
xlim([0 Lx])
ylim([0 Ly])
filename = sprintf('vn_re%d_%ds.png', Re, i*0.5);
filepath = fullfile(directory, filename);
saveas(gcf, filepath);

% Vorticity contour
[X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));
u = u_sample(2:Ny+1, 2:Nx+1, i);
v = v_sample(2:Ny+1, 2:Nx+1, i);
[w1,w2] = curl(X,Y,u,v);
figure
contourf(X, Y, w1, 80, 'LineStyle', 'none');
load('colorData.mat','BuDRd_18')
colormap(othercolor('BuDRd_18'));
shading flat;
colorbar('southoutside');
setPivot(0)
axis equal
title(['Vorticity contour at ', sprintf('Re=%d ',Re), '(resolution: ', sprintf('%.d',Nx), '\times', sprintf('%d)',Ny)])
filename = sprintf('vor_re%d_%ds.png', Re, i*0.5);
filepath = fullfile(directory, filename);
saveas(gcf, filepath);

close all;

%% t vs Drag Coefficient
[bdsurf,bdn,bds,bde,bdw]=bodycontact(1-ifluid);
for i = 1:n/rsampling
    pre=p_sample(:,:,i);
    dragw=bdw.*pre;
    drage=bde.*pre;
    fdrag=sum(sum(dragw))*dy-sum(sum(drage))*dy;
    cd(i)=(-2*fdrag);
    t(i)=0.5*i;
end
figure;
plot(t,cd)
xlabel('t');
ylabel('C_D');
title('Pressure drag coefficient variation with time')
cdpm=mean(mean(cd));

%% t vs Lift Coefficient
[bdsurf,bdn,bds,bde,bdw]=bodycontact(1-ifluid);
for i = 1:n/rsampling
    pre=p_sample(:,:,i);
    drags=bds.*pre;
    dragn=bdn.*pre;
    fdrag=sum(sum(drags))*dy-sum(sum(dragn))*dy;
    cl(i)=(-2*fdrag);
    tl(i)=0.5*i;
end
figure;
plot(tl,cl)
xlabel('t');
ylabel('C_L');
title('Lift coefficient variation with time')
cdlm = mean(mean(cl));

%% Frequency
[bdsurf,bdn,bds,bde,bdw]=bodycontact(1-ifluid);
p_post = zeros(Ny+2,Nx+2,100);
tp = zeros(100,1);
cdp = zeros(100,1);
p_post(:,:,1:100) = p_sample(:,:,101:1:200); % Extract 50s~100s
for i = 1:100 
    pre=p_post(:,:,i);
    dragw=bdw.*pre;
    drage=bde.*pre;
    fdrag=sum(sum(dragw))*dy-sum(sum(drage))*dy;
    cdp(i)=(-2*fdrag);
    tp(i)=0.5*(i+100);
end
figure;
plot(tp,cdp)
xlabel('tp');
ylabel('C_{Dp}');
title('Pressure drag coefficient variation with time')

L=length(cdp);
Fs = 2;
f = Fs*(0:(L/2))/L;
Y=fft(cdp);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure
plot(f(2:end),P1(2:end)) 
title('Spectrum of pressure drag coefficient')
xlabel("Frequency (f_D)")
ylabel("Spectrum")

%% Functions
function setPivot(varargin)
if nargin==0
    ax=gca;pivot=0;
else
    if isa(varargin{1},'matlab.graphics.axis.Axes')
        ax=varargin{1};
        if nargin>1
            pivot=varargin{2};
        else
            pivot=0;
        end
    else
        ax=gca;pivot=varargin{1};
    end
end
try
    CLimit=get(ax,'CLim');
catch
end
try
    CLimit=get(ax,'ColorLimits');
catch
end
CMap=colormap(ax);
CLen=[pivot-CLimit(1),CLimit(2)-pivot];
if all(CLen>0)
    [CV,CInd]=sort(CLen);
    CRatio=round(CV(1)/CV(2).*300)./300;
    CRatioCell=split(rats(CRatio),'/');
    if length(CRatioCell)>1
        Ratio=[str2double(CRatioCell{1}),str2double(CRatioCell{2})];
        Ratio=Ratio(CInd);
        N=size(CMap,1);
        CList1=CMap(1:floor(N/2),:);
        CList2=CMap((floor(N/2)+1):end,:);
        if mod(N,2)~=0
            CList3=CList2(1,:);CList2(1,:)=[];
            CInd1=kron((1:size(CList1,1))',ones(Ratio(1)*2,1));
            CInd2=kron((1:size(CList2,1))',ones(Ratio(2)*2,1));
            CMap=[CList1(CInd1,:);repmat(CList3,[Ratio(1)+Ratio(2),1]);CList2(CInd2,:)];
        else
            CInd1=kron((1:size(CList1,1))',ones(Ratio(1),1));
            CInd2=kron((1:size(CList2,1))',ones(Ratio(2),1));
            CMap=[CList1(CInd1,:);CList2(CInd2,:)];
        end
        colormap(ax,CMap);
    end
end
end

function [bdsurf,bdn,bds,bde,bdw]=bodycontact(Geobd)
    nx=size(Geobd,2); ny=size(Geobd,1);
    bdn=[Geobd(2:end,:);zeros(1,nx)]-Geobd;
    bds=[zeros(1,nx);Geobd(1:end-1,:)]-Geobd;
    bde=[Geobd(:,2:end),zeros(ny,1)]-Geobd;
    bdw=[zeros(ny,1),Geobd(:,1:end-1)]-Geobd;
    bdn(bdn<0)=0;bds(bds<0)=0;
    bde(bde<0)=0;bdw(bdw<0)=0;
    bdsurf=bdn+bds+bde+bdw; bdsurf(bdsurf>0)=1;
end