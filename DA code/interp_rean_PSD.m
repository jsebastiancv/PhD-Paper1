clear 
clc
%close all

%% Settings
addpath([getenv('HOME'),'/Documents/VERB/3D VERB DA/Code/Various_functions/']);
addpath([getenv('HOME'),'/Documents/datalibrary/core/read/']);

target_K = 0.11; % G^0.5 Re
target_mu = 300; % MeV/G
angle = pi/180;
mc2 = 0.511;
fastint = false; % will find nearest values only

sDate = datenum('01-Oct-2012'); eDate = datenum('01-Nov-2012');
 
strdate = datestr(sDate, 'yyyymm');
mfm = 'TS07Dmid15';

% Load reanalysis file
fileName = ['reanalysis/reanalysis_',mfm,'_final/Reanalysis_EQE_LatestVERB_noMP_',strdate,'_Gaussian_onera_',mfm,'.mat'];
data = load(fileName);
alpha = data.SimInvAlpha;
energy = data.SimInvEnergy;
K = data.SimInvK;
mu = data.SimInvMu;
L = data.SimL;

time = data.SimTime;
PSD = data.SimPSD;

Kp = data.Kp;
sz = size(PSD);

PSD_new = nan(sz(1),sz(2));

% Interpolation
    if fastint
        fprintf('warning, using coarse interpolation\n')
        linx = floor(size(L,1)/2); % we need to do this on L,mu,alpha grid (approx at middle of grid)
        [~,kinx] = min(abs(target_K - K(linx,1,:)));
        tmpMu = squeeze(mu(:,:,kinx));
        [~,minx] = min(abs(target_mu - tmpMu(linx,:)));
        PSD_new = squeeze(PSD(:,:,minx,kinx));

        % Redefine values
        fprintf('redefining target mu and k to match nearest values\n')
        target_mu = round(tmpMu(linx,minx)); % mu is the same at all L
        target_K = K(linx,1,kinx); %K varies due to fixed alpha
    else
        fprintf('Interpolating reanalysis (this may take some time)...\n');
        t1=tic;

        PSD = log10(PSD);
        K = K.^(1/3);
        mu = log10(mu);

        for it=1:sz(1)
            tpsd = nan(1,sz(2));
            PSDt = squeeze(PSD(it,:,:,:));
            for il=1:sz(2)
                psdvec = squeeze(PSDt(il,:,:));
                muvec= squeeze(mu(il,:,:));
                kvec = squeeze(K(il,:,:));
                test = griddata(muvec,kvec,psdvec,log10(target_mu),target_K.^(1/3),'natural');
                tpsd(il) = test;
            end
            PSD_new(it,:) = tpsd;
            if mod(it,100) == 0
                t2 = toc(t1);
                ltime = t2/it;
                est = (sz(1)-it) * ltime;
                fprintf('%i/%i Estimated Time: %s min\n',it,sz(1),num2str(est/60));
            end
        end
        PSD_new = 10.^PSD_new;
        fprintf('\n');
        fprintf('done\n');
    end    
    
psd_rean = PSD_new ./ 2.997e7; %VERB units to (c/cm/MeV).^3
    
%% Save file
save('psdfile_rean.mat','psd_rean')

