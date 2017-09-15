close all
clear all
addpath(genpath(' '))

DIR = pwd;

RESO= 600;
scale = 0;
Nx = 32;
Ny = 56;


Smin = 992;
Smax = 1056;
Tmin = 22;
Tmax = 38;



alpha = 0.15e-6;
R = 0.04; %% 4cm
vot_cha = alpha/(R^2); %% characteristic vorticity scale

V0 = 15e-2;


xlimit = [-0.05 1.02]*R;
ylimit = [-0.05 1.8]*R;
%set(gca,'yaxislocation','right')

%line([0 1 1],[0 0 1.75],'linewidth',2,'color','k')

Time = [0.03,0.055,0.085,0.11,0.14,0.17];
Minute = [5, 10, 15, 20, 25, 30];

for it = 1:length(Time)
time = Time(it);
minute = Minute(it);
file = ['data_u_w_T_S_vort_Nu0.9_t' num2str(time) '_' num2str(minute) 'min.txt'];
data = [ file ];
a = load(data);

x = a(:,1);
y = a(:,2);
u = a(:,3);
v = a(:,4);
t = a(:,5);
s = a(:,6);
vt = a(:,7);

xmin = min(x(:));
xmax = max(x(:));
ymin = min(y(:));
ymax = max(y(:));

x1d = linspace(xmin,xmax,Nx);
y1d = linspace(ymin,ymax-1e-4,Ny);

[X, Y] = meshgrid(x1d,y1d);
U = griddata(x,y,u,X,Y);
V = griddata(x,y,v,X,Y);
T = griddata(x,y,t,X,Y);
S = griddata(x,y,s,X,Y);
VT = griddata(x,y,vt,X,Y);

U = U * alpha / R;
V = V * alpha / R;
X = X * R;
Y = Y * R;

U = U / V0;
V = V / V0;

Xm = -X;
Tm = flip(T,2);
%Sm = flip(S,2);

Sm_Dim = S*(Smax-Smin) + Smin; %% dimensional value of salinity

% figure
% colormap (copper)
% hS = contourf(Xm,Y,Sm_Dim,100,'linecolor','none');
% hbar = colorbar;
% axis equal
% xlim(-flip(xlimit))
% ylim(ylimit)
% set(gca,'yaxislocation','right')
% 
% %line([0 -1 -1],[0 0 1.75],'linewidth',2,'color','k')
% set(gca,'visible','off')
% set(gcf,'color','w')
% pdf =[DIR filesep 'salinity_time' num2str(time)  '_' num2str(RESO)];
% 
% expr = ['export_fig ' pdf ' -png -eps -pdf -transparent -r' num2str(RESO)];
% eval(expr)

figure
colormap (parula)
hold on
hVT = contourf(X,Y,VT*vot_cha,100,'linecolor','none');
quiver(X,Y,U,V,scale,'color','black')
axis equal
hbar = colorbar;
set(hbar,'ylim',[-1.2 0.8])
%  set(hbar,'fontsize',20,'fontname','times new roman','position',[0.8 0.13 0.035 0.77],'ylim',[-11052 7168])  
%  caxis([-11052 7168])
xlim(xlimit)
ylim(ylimit)
set(gca,'yaxislocation','right')

%line([0 1 1],[0 0 1.75],'linewidth',2,'color','k')

% [hx,hy] = format_ticks(gca,{'0','0.5','1'},...
%                               {'0','0.5','1','1.5'},...
%                               [[0 0.5 1]],[0 0.5 1 1.5],[],[],[0.0],[-1.2]);
% set(hx,'fontname','times new roman','fontsize',20)    
% set(hy,'fontname','times new roman','fontsize',20)  
set(gca,'visible','off')
set(gcf,'color','w')
pdf =[DIR filesep 'vector_vorticityPhi_time' num2str(time) '_' num2str(RESO)];

expr = ['export_fig ' pdf ' -png -eps -pdf -transparent -r' num2str(RESO)];
eval(expr)
end

