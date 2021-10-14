clear 
clc
close all

%% Settings
addpath([getenv('HOME'),'/Documents/VERB/3D VERB DA V4/Code/Various_functions/']);
addpath([getenv('HOME'),'/Documents/datalibrary/core/read/']);

%% Plot PSD
psd_rean = load('psdfile_rean.mat');

time = linspace(1,8761,8761);

d_time = time(2) - time(1);
time_range = [time(1) time(end)];

tticks_month = [1,745,1465,2209,2953,3625,4369,5089,5833,6553,7297,8041,8761];
L = linspace(1,6.6,29);
 
figure
subplot(6,1,1)
cmap = jet(256);
% cmap(1,:) = [1 1 1];
colormap(cmap);

[ax1, cb1] = contourf_psd(time(1:end),L,psd_rean.psd_LT(1:8761,:));

set(ax1,'xlim',[time_range(1)-0.5*(time_range(2)-time_range(1)) time_range(end)]);

box on

date = datestr(tticks_month);
set(ax1,'XTick',tticks_month,'layer','top');
axis tight

caxis([-10 -6])

set(ax1,'xticklabel',[]) 
set(ax1,'YTick',[3,4,5,6]);

xlabel('Date (2012-2013)');
hold on
    
ylims = [3 6.6];
ylim(ylims)

ylabel('L*')%,'FontSize',14)%,'FontWeight','bold')
hold on
title({'Electron PSD','3D diffusion + mixed pitch angle - energy diffusion + scattering by EMICs + magnetopause shadowing'},'FontWeight','bold','FontSize',18); 
text(200,4, '\mu = 3500 MeV/G, K = 0.01 G^{0.5} Re','FontWeight','bold','FontSize',14); 
set(ax1,'FontSize',14,'FontWeight','bold')

%% Kp and Dst
sDate = datenum('01-Oct-2012');
eDate = datenum('01-Jan-2013');

subplot(4,1,2)

[tkp, kp] = omniwebdata(sDate, eDate, 38);
[tdst, dst] = omniwebdata(sDate, eDate, 40);
time= tkp;

[hAx,~,~] = plotyy(tkp,kp.arr/10,tdst,dst.arr);
grid minor
ax2 = gca;

ylabel(hAx(1),'Kp','FontSize',14,'FontWeight','bold');
ylabel(hAx(2),{'Dst','(nT)'},'FontSize',14,'FontWeight','bold');
tticks_month_data = time(tticks_month);

set(hAx,'XTick',tticks_month_data','layer','top');
set(gca,'XTickLabel',['01-Oct';'01-Nov';'01-Dec';'01-Jan';'01-Feb';'01-Mar';'01-Apr';'01-May';'01-Jun';'01-Jul';'01-Aug';'01-Sep';'01-Oct']) 

set(hAx,'xlim',[time(1)-0.5*(time(2)-time(1)) time(end)]);

set(hAx(1),'YLim',[0 9],'FontSize',14,'FontWeight','bold');
set(hAx(2),'YLim',[-200 50],'FontSize',14,'FontWeight','bold')

xlabel('Date (2012-2013)');
ax2.Position(3) = ax1.Position(3);