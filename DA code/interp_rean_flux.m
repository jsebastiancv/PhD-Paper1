clear 
clc
%close all

%% Settings
addpath([getenv('HOME'),'/Documents/VERB/3D VERB DA/Code/Various_functions/']);
addpath([getenv('HOME'),'/Documents/datalibrary/core/read/']);

target_energy = 1; % MeV
target_alpha = 50*angle; % degree
angle = pi/180;
mc2 = 0.511;
fastint = true; % will find nearest values only

sDate = datenum('01-Oct-2015'); eDate = datenum('01-Nov-2015'); 

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
        for il=1:sz(2)
            target_L = L(il,1,1);
            target_pc = sqrt((target_energy./mc2+1).^2-1).*mc2;
            target_K = Lalpha2K(target_L,target_alpha); % G^0.5 Re 
            target_mu = pc2mu(target_L, target_pc, target_alpha); % MeV/G
            linx = floor(size(L,1)/2); % we need to do this on L,mu,alpha grid (approx at middle of grid)
            [~,kinx] = min(abs(target_K - K(linx,1,:)));
            tmpMu = squeeze(mu(:,:,kinx));
            [~,minx] = min(abs(target_mu - tmpMu(linx,:)));
            PSD_new = squeeze(PSD(:,il,minx,kinx));  

            % Redefine values
            fprintf('redefining target mu and k to match nearest values\n')
            target_mu = round(tmpMu(linx,minx)); % mu is not the same at different  L
            target_K = K(linx,1,kinx); % K varies due to fixed alpha
            flux_rean(1:sz(1),il) = PSD_new.*target_pc.^2;
        end    
        else
            fprintf('interpolating reanalysis (this may take some time)...\n');
            PSD = log10(PSD);
            K = K.^(1/3);
            mu = log10(mu);

            for il=1:sz(2)
                t1=tic;
                target_L = L(il,1,1);
                target_pc = sqrt((target_energy./mc2+1).^2-1).*mc2;
                target_K = Lalpha2K(target_L,target_alpha); % G^0.5 Re 
                target_mu = pc2mu (target_L, target_pc, target_alpha); % MeV/G
                for it=1:sz(1)
                    tpsd = nan(1,sz(2));
                    PSDt = squeeze(PSD(it,:,:,:));
                    psdvec = squeeze(PSDt(il,:,:));
                    muvec= squeeze(mu(il,:,:));
                    kvec = squeeze(K(il,:,:));
                    test = griddata(muvec,kvec,psdvec,log10(target_mu),target_K.^(1/3),'natural');
                    PSD_new(it,il) = test;
                end
                t2 = toc(t1);
                ltime = t2;
                est = (sz(2)-il) * ltime;
                fprintf('%i/%i Estimated Time: %s sec\n',il,sz(2),num2str(est));
            end
                flux_rean = 10.^PSD_new.*target_pc.^2; 
                fprintf('\n');
 
   
    end

%% Save file
save('fluxfile_rean.mat','flux_rean')
