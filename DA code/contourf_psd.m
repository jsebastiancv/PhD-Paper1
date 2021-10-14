function [ ax, cb] = contourf_psd( time, Lstar, PSD )
%CONTOURF_FLUX plot Lstar vs time Flux
ax = gca;

% Unconditionaly positive
PSD = abs(PSD);

% Set the "bottom" level of the PSD. For plotting
c_max = -4;
c_min = -8; % number of orders of magnitude to plot
    
PSD(PSD < 10^c_min) = 10^c_min;

contourf(time, Lstar, log10(PSD)', 100,'LineStyle','none');
shading flat
axis tight
hold on;
% colorbar 
cm = colormap(jet);
% cm(1,:) = [0 0 0];
colormap(cm);
% set(ax,'xlim',[0 100],'xtick',0:10:100,'ylim',[1 6.6],'ytick',1:7,'clim',[c_min, c_max]);
ylabel('L^*');
drawnow;

% get position
P = ax.Position;
cb = colorbar; ylabel(cb,{'log_{10} PSD','(c/cm/MeV)^3'});
% draw colorbar;
drawnow;


% Adjust positions;
% cb.Position(1) = P(1)+P(3)+0.01;
% ax.Position = P;
% drawnow;

end


