% Multi-phase LBM: Coalescence and Splitting in a Microfluidic Channel
% Written by: Ahmet Burak Yıldırım & Muhammed Uygun, 2023

clear; close all; clc

lx           = 127*1;
ly           = 127*1;
filename     = "channel.gif"; f = figure; set(gcf,'color','w');


G = -0.8;              % Amplitude of the molecular interaction force

omega_droplet = 1;     % Relaxation parameter for droplet
omega_fluid = 1;       % Relaxation parameter for fluid

maxT   = 3000;  
tplot  = 5;      

% D2Q9
W   = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
CX  = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
CY  = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];

[y,x] = meshgrid(1:ly,1:lx);

%%%%%%%%%%%%%%%%%%%%%%%% CHANNEL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WT = 5;         % thickness of the upper and lower walls
entryL = 30;    % length of the entry channels
entryT = 20;    % thickness of the entry channels
jL = 5;         % length of the junction
exitT = 21;     % thickness of the outlet channel
dm = 5;
shape = zeros(lx,ly);
shape(dm:end-dm,1:WT) = 1;          % upper wall
shape(dm:end-dm,end-WT:end) = 1;    % lower wall
shape(dm:entryL,WT+entryT:ly-(WT+entryT)) = 1;
shape(entryL+jL:end-dm,1:(ly-exitT)/2) = 1;
shape(entryL+jL:end-dm,end-(ly-exitT)/2:end) = 1;

% exit channel upper triangle
for index = (entryL+jL):entryL+jL+((ly-(2*WT+exitT))/2)
    shape(entryL+jL:index,WT+1+index-(entryL+jL)) = 0;
end
% exit channel lower triangle
for index = (entryL+jL):entryL+jL+((ly-(2*WT+exitT))/2)
    shape(entryL+jL:index,(ly-WT)-(index-(entryL+jL)+1)) = 0;
end
% seperating wall upper triangle
for index = (entryL):entryL+((ly-2*WT-2*entryT+1)/2)
    shape(entryL:index,(WT+entryT)+index-(entryL)) = 1;
end
% seperating wall lower triangle
for index = (entryL):entryL+((ly-2*WT-2*entryT+1)/2)
    shape(entryL:index,(ly-(WT+entryT))-(index-entryL)) = 1;
end
grid = zeros(lx,ly);
wall = logical(insert(grid,shape));

drho = 1;
delta_rho = -drho*(1-2.0*rand(lx));

center1 = [15, 15];
center2 = [15, ly-15];

for dx = 1:lx
    for dy = 1:ly
        delta_rho(dx,dy) = drho*-1;
        d1 = sqrt((dx-center1(1))^2+(dy-center1(2))^2);
        d2 = sqrt((dx-center2(1))^2+(dy-center2(2))^2);

        if d1<8
            delta_rho(dx,dy) = drho*1;
        end
        if d2<8
            delta_rho(dx,dy) = drho*1;
        end
       
    end
end

for i=1:9
    dropletIn(i,1:lx,1:ly) = W(i).*(1.0 + delta_rho);
    fluidIn(i,1:lx,1:ly) = W(i).*(1.0 - delta_rho);
end

dropletIn = reshape(dropletIn,[9 ly*lx]); dropletIn = dropletIn'; 
fluidIn = reshape(fluidIn,[9 ly*lx]); fluidIn = fluidIn'; 

droplet_wall = dropletIn(wall,[1 6 7 8 9 2 3 4 5]); % Bounce-back boundary condition
fluid_wall = fluidIn(wall,[1 6 7 8 9 2 3 4 5]);     % Bounce-back boundary condition

dropletIn(wall,:) = droplet_wall;
fluidIn(wall,:) = fluid_wall;

dropletIn = dropletIn'; dropletIn = reshape(dropletIn,[9, ly, lx]); 
fluidIn = fluidIn'; fluidIn = reshape(fluidIn,[9, ly ,lx]);

Gomega_droplet = G/omega_droplet;
Gomega_fluid = G/omega_fluid;
for cycle = 1:maxT
    rho_droplet = sum(dropletIn);
    rho_fluid = sum(fluidIn);
    jx_droplet  = reshape ( (CX * reshape(dropletIn,9,lx*ly)), 1,lx,ly);
    jy_droplet  = reshape ( (CY * reshape(dropletIn,9,lx*ly)), 1,lx,ly);
    jx_fluid  = reshape ( (CX * reshape(fluidIn,9,lx*ly)), 1,lx,ly);
    jy_fluid  = reshape ( (CY * reshape(fluidIn,9,lx*ly)), 1,lx,ly);
   
    rhoTot_OMEGA = rho_droplet*omega_droplet + rho_fluid*omega_fluid;
    uTotX = (jx_droplet*omega_droplet+jx_fluid*omega_fluid) ./ rhoTot_OMEGA;
    uTotY = (jy_droplet*omega_droplet+jy_fluid*omega_fluid) ./ rhoTot_OMEGA;
	
    rhoContrib1x = 0.0;
    rhoContrib2x = 0.0;
    rhoContrib1y = 0.0;
    rhoContrib2y = 0.0;

    for i=2:9
        rhoContrib1x = rhoContrib1x + circshift(rho_droplet*W(i), [0,CX(i),CY(i)])*CX(i);
        rhoContrib1y = rhoContrib1y + circshift(rho_droplet*W(i), [0,CX(i),CY(i)])*CY(i);
        
        rhoContrib2x = rhoContrib2x + circshift(rho_fluid*W(i), [0,CX(i),CY(i)])*CX(i);
        rhoContrib2y = rhoContrib2y + circshift(rho_fluid*W(i), [0,CX(i),CY(i)])*CY(i);
    end
    
    %POTENTIAL CONTRIBUTION OF FLUID ON DROPLET
    uTotX1 = uTotX - Gomega_droplet.*rhoContrib2x;
    uTotY1 = uTotY - Gomega_droplet.*rhoContrib2y;

    %POTENTIAL CONTRIBUTION OF DROPLET ON FLUID
    uTotX2 = uTotX - Gomega_fluid.*rhoContrib1x;
    uTotY2 = uTotY - Gomega_fluid.*rhoContrib1y;

    % To introduce velocity
    if cycle < 1300
        uTotX2(1,1,:) = uTotX2(1,1,:) + 3*0.05;
    end
    if cycle > 1300 && cycle < 1800
        uTotX2(1,1,:) = uTotX2(1,1,:);
    end
    if cycle > 1800
        uTotX2(1,1,:) = uTotX2(1,1,:) - 3*0.05;
    end

    uTotX1(wall)=0; uTotY1(wall)=0;
    uTotX2(wall)=0; uTotY2(wall)=0;

   for i=1:9
      cuNS_droplet = 3*(CX(i)*uTotX1+CY(i)*uTotY1);
      cuNS_fluid = 3*(CX(i)*uTotX2+CY(i)*uTotY2);
      
      dropletEq(i,:,:)   = rho_droplet .* W(i) .* ...
                       ( 1 + cuNS_droplet + 0.5*(cuNS_droplet.*cuNS_droplet) - 1.5*(uTotX1.^2+uTotY1.^2) );
      fluidEq(i,:,:)   = rho_fluid .* W(i) .* ...
                       ( 1 + cuNS_fluid + 0.5*(cuNS_fluid.*cuNS_fluid) - 1.5*(uTotX2.^2+uTotY2.^2) );
      dropletOut(i,:,:)  = dropletIn(i,:,:) - omega_droplet .* (dropletIn(i,:,:)-dropletEq(i,:,:));
      fluidOut(i,:,:)  = fluidIn(i,:,:) - omega_fluid .* (fluidIn(i,:,:)-fluidEq(i,:,:));
   end

    dropletOut = reshape(dropletOut,[9 ly*lx]); dropletOut = dropletOut'; 
    fluidOut = reshape(fluidOut,[9 ly*lx]); fluidOut = fluidOut'; 
    droplet_wall = dropletOut(wall,[1 6 7 8 9 2 3 4 5]); % Bounce-back boundary condition
    fluid_wall = fluidOut(wall,[1 6 7 8 9 2 3 4 5]);     % Bounce-back boundary condition
    dropletOut(wall,:) = droplet_wall;
    fluidOut(wall,:) = fluid_wall;
    dropletOut = dropletOut'; dropletOut = reshape(dropletOut,[9, ly, lx]); 
    fluidOut = fluidOut'; fluidOut = reshape(fluidOut,[9, ly ,lx]);

   for i=1:9
      dropletIn(i,:,:) = circshift(dropletOut(i,:,:), [0,CX(i),CY(i)]);
      fluidIn(i,:,:) = circshift(fluidOut(i,:,:), [0,CX(i),CY(i)]);
   end
    picturewidth = 15;
    hw_ratio = 0.5;
   if(mod(cycle,tplot)==0)
       
       % colorbar limits without altering the colorbar
       largest = 0.4+1e-3;
       smallest = 0;

       rho_fluid     = reshape(rho_fluid,lx,ly);
       Ux = reshape(uTotX2,lx,ly);
       Uy = reshape(uTotY2,lx,ly);
       u_fluid       = sqrt(Ux.^2+Uy.^2);
       u_fluid = u_fluid/max(max(u_fluid))*0.4;
       ufluidt = u_fluid + (wall) * largest; % recoloring walls differently
       imagesc(ufluidt'); 
       tx = ["$t$ = " + num2str(cycle)];
       cb = colorbar; cb.Location = "southoutside";


       indexValue = 0.4;               % value for which to set a particular color
       topColor = [0.2 0.2 0.2];       % color for maximum data value (red = [1 0 0])
       indexColor = [0.2 1 0.2];       % color for indexed data value (white = [1 1 1])
       bottomcolor = [0.2 0.2 1];      % color for minimum data value (blue = [0 0 1])

       L = size(ufluidt,1);
       index = L*abs(indexValue-smallest)/(largest-smallest);
       customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
           linspace(bottomcolor(2),indexColor(2),100*index)',...
           linspace(bottomcolor(3),indexColor(3),100*index)'];
       customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
           linspace(indexColor(2),topColor(2),100*(L-index))',...
           linspace(indexColor(3),topColor(3),100*(L-index))'];
       customCMap = [customCMap1;customCMap2];  % combining colormaps
       h = colormap(customCMap); caxis([smallest largest])

       title(tx);
       axis equal off; 
       set(f, 'Units', 'centimeters', 'Position', [2+5 1+5 picturewidth/2+5 1*picturewidth/2+5])
       set(findall(f,'-property', 'Box'), 'Box', 'on') 
       set(findall(f, '-property', 'Interpreter'), 'Interpreter', 'latex')
       set(findall(f, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
       set(findall(f,'-property', 'FontSize'), 'FontSize', 20) % never change fontsize anymore!
       tx = ["FIG"+cycle];
       drawnow
       if cycle == 10
           frame = getframe(f);
           im = frame2im(frame);
           [imind,cm] = rgb2ind(im,256);
           imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'Loopcount',inf);
       else
           frame = getframe(f);
           im = frame2im(frame);
           [imind,cm] = rgb2ind(im,256);
           imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'WriteMode','append');
       end
   end
end

function out = insert(N,shape)
    small = size(shape,1);
    N((end+1)/2-(small-1)/2:(end+1)/2+(small-1)/2,(end+1)/2-(small-1)/2:(end+1)/2+(small-1)/2) = ...
        N((end+1)/2-(small-1)/2:(end+1)/2+(small-1)/2,(end+1)/2-(small-1)/2:(end+1)/2+(small-1)/2) + shape;
    out = N;
end