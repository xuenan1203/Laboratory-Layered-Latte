close all
clear all
addpath(genpath(' '))

load('Fig. 2c_data.mat');

RESO= 600;
scale = 0;
R = 0.04;

% xmin = min(x(:));
% xmax = max(x(:));
% ymin = min(y(:));
% ymax = max(y(:));
% 
% x1d = linspace(xmin,xmax);
% y1d = linspace(ymin,ymax-1e-4);
% 
% [X, Y] = meshgrid(x1d,y1d);
% U = griddata(x,y,u,X,Y);
% V = griddata(x,y,v,X,Y);

V0 = 15e-2;

X = x * 1e-3;
Y = y * 1e-3;
U = u * 1e-3/V0;
V = v * -1e-3/V0;


%xlimit = [-0.005 0.08];
%ylimit = [0 0.08];


figure;
quiver(X,Y,U,V,scale,'color','black')
set(gca,'visible','off')
%xlim(xlimit)
%ylim(ylimit)
%axis equal

eval('export_fig [ADDRESS] -png -eps -pdf -transparent -r600')
