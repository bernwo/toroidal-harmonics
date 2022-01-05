%% CLEAR EVERYTHING
clear
close all
clc
path_to_exports = 'C:/Users/jen_7/Desktop/Toroidal Harmonics/';
%%
close all
a = 1;
r = 1;
Ln = 50;
sigma = linspace(-pi,pi,Ln);
phi = linspace(0,2*pi,Ln);
tau = zeros(1,Ln)+1/r;

[Sigma,Tau,Phi] = meshgrid(sigma,tau,phi);

x = (a.*sinh(Tau).*cos(Phi)) ./ (cosh(Tau)-cos(Sigma)); 
y = (a.*sinh(Tau).*sin(Phi)) ./ (cosh(Tau)-cos(Sigma));
z = (a.*sin(Sigma)) ./ (cosh(Tau)-cos(Sigma)); 

% writematrix([Sigma(:) Tau(:) Phi(:)],'toroidal_coordinates.csv'); %writes .csv
%%
Phi00 = table2array(readtable([path_to_exports 'Toroidal_Harmonics00.csv']));
Phi00 = Phi00(:,1) + Phi00(:,2)*1i;
max00 = max(abs(Phi00));

Phi01 = table2array(readtable([path_to_exports 'Toroidal_Harmonics01.csv']));
Phi01 = Phi01(:,1) + Phi01(:,2)*1i;
max01 = max(abs(Phi01));

Phi11 = table2array(readtable([path_to_exports 'Toroidal_Harmonics11.csv']));
Phi11 = Phi11(:,1) + Phi11(:,2)*1i;
max11 = max(abs(Phi11));

Phi12 = table2array(readtable([path_to_exports 'Toroidal_Harmonics12.csv']));
Phi12 = Phi12(:,1) + Phi12(:,2)*1i;
max12 = max(abs(Phi12));

Phi21 = table2array(readtable([path_to_exports 'Toroidal_Harmonics21.csv']));
Phi21 = Phi21(:,1) + Phi21(:,2)*1i;
max21 = max(abs(Phi21));

Phi22 = table2array(readtable([path_to_exports 'Toroidal_Harmonics22.csv']));
Phi22 = Phi22(:,1) + Phi22(:,2)*1i;
max22 = max(abs(Phi22));
%% PLOTTING
myview = [15 35];
close all
fig1 = figure('WindowState','maximized');
ax1=subplot(2,3,1);
scatter3(x(:),y(:),z(:),20, real(Phi00) , 'filled'); 
caxis([-max00 max00])
axis equal
grid off;
rotate3d on;
colormap(jet);
colorbar();
title('\nu=0,\mu=0')
view(myview)

ax2=subplot(2,3,2);
scatter3(x(:),y(:),z(:),20, real(Phi01) , 'filled'); 
caxis([-max01 max01])
axis equal
grid off;
rotate3d on;
colormap(jet);
colorbar();
title('\nu=0,\mu=1')
view(myview)

ax3=subplot(2,3,3);
scatter3(x(:),y(:),z(:),20, real(Phi11) , 'filled'); 
caxis([-max11 max11])
axis equal
grid off;
rotate3d on;
colormap(jet);
colorbar();
title('\nu=1,\mu=1')
view(myview)

ax4=subplot(2,3,4);
scatter3(x(:),y(:),z(:),20, real(Phi12) , 'filled'); 
caxis([-max12 max12])
axis equal
grid off;
rotate3d on;
colormap(jet);
colorbar();
title('\nu=1,\mu=2')
view(myview)

ax5=subplot(2,3,5);
scatter3(x(:),y(:),z(:),20, real(Phi21) , 'filled'); 
caxis([-max21 max21])
axis equal
grid off;
rotate3d on;
colormap(jet);
colorbar();
title('\nu=2,\mu=1')
view(myview)

ax6=subplot(2,3,6);
scatter3(x(:),y(:),z(:),20, real(Phi22) , 'filled'); 
caxis([-max22 max22])
axis equal
grid off;
rotate3d on;
colormap(jet);
colorbar();
title('\nu=2,\mu=2')
view(myview)

%% GENERATE GIF
pi_range_len = 30;
pi_range = linspace(0,2*pi,pi_range_len);
giffilename = 'toroidal_harmonics.gif';

for iter = 1:(pi_range_len-1)
    h = figure('WindowState','maximized');
    
    ax1=subplot(2,3,1);
    scatter3(x(:),y(:),z(:),20, real(Phi00*exp(-1i*pi_range(iter))) , 'filled'); 
    caxis([-max00 max00])
    axis equal
    grid off;
    rotate3d on;
    colormap(jet);
    colorbar();
    title('\nu=0,\mu=0')
    view(myview)

    ax2=subplot(2,3,2);
    scatter3(x(:),y(:),z(:),20, real(Phi01*exp(-1i*pi_range(iter))) , 'filled'); 
    caxis([-max01 max01])
    axis equal
    grid off;
    rotate3d on;
    colormap(jet);
    colorbar();
    title('\nu=0,\mu=1')
    view(myview)

    ax3=subplot(2,3,3);
    scatter3(x(:),y(:),z(:),20, real(Phi11*exp(-1i*pi_range(iter))) , 'filled'); 
    caxis([-max11 max11])
    axis equal
    grid off;
    rotate3d on;
    colormap(jet);
    colorbar();
    title('\nu=1,\mu=1')
    view(myview)

    ax4=subplot(2,3,4);
    scatter3(x(:),y(:),z(:),20, real(Phi12*exp(-1i*pi_range(iter))) , 'filled'); 
    caxis([-max12 max12])
    axis equal
    grid off;
    rotate3d on;
    colormap(jet);
    colorbar();
    title('\nu=1,\mu=2')
    view(myview)

    ax5=subplot(2,3,5);
    scatter3(x(:),y(:),z(:),20, real(Phi21*exp(-1i*pi_range(iter))) , 'filled'); 
    caxis([-max21 max21])
    axis equal
    grid off;
    rotate3d on;
    colormap(jet);
    colorbar();
    title('\nu=2,\mu=1')
    view(myview)

    ax6=subplot(2,3,6);
    scatter3(x(:),y(:),z(:),20, real(Phi22*exp(-1i*pi_range(iter))) , 'filled'); 
    caxis([-max22 max22])
    axis equal
    grid off;
    rotate3d on;
    colormap(jet);
    colorbar();
    title('\nu=2,\mu=2')
    view(myview)

    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if iter == 1 
      imwrite(imind,cm,giffilename,'gif','DelayTime',0.08,'Loopcount',inf); 
    else 
      imwrite(imind,cm,giffilename,'gif','DelayTime',0.08,'WriteMode','append'); 
    end
    close
    disp(iter)
end