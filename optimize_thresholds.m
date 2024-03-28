clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kVp = 120; 
frame_time = 300e-6;
mA = 200;

design = '60mmActiveSi/';

% Load energy vectors and spectrum

m_settings = matfile([design 'SimulationSettings.mat']);

E = m_settings.eActVectorkeV;
D = m_settings.eDepVectorkeV;

N_inc = length(E);
N_dep = length(D);

pre_patient_spectrum = m_settings.backgroundSpectrumBeforePatientmm_2';
% normalized_pre_patient_spectrum = pre_patient_spectrum / sum(pre_patient_spectrum); %Commented by MP
% Compute electronic noise variance contributed by each bin (added by MP)
electronic_noise_data=load('./electronic_noise_simulations/electronic_noise_data_1.6keV.mat');
enoise_diff_variance_unbinned=-diff([interp1(electronic_noise_data.thresholds_for_electronic_noise_variance_keV+0.5,... %0.5 takes into account that the bin thresholds are really 0.5 keV below the integer-valued  threshold energy in this code
    electronic_noise_data.electronic_noise_variance,D,'linear',0),0]);
normalized_pre_patient_spectrum= electronic_noise_data.pre_patient_photons_per_mm2*pre_patient_spectrum / sum(pre_patient_spectrum); %"normalized" is a misnomer since this does not sum to one anymore
% Load response

    
eNoiseSigmakeV = 1.6;

%MP commented this since Fano noise is not included, we can load an existing file instead
% halfw = ceil(eNoiseSigmakeV*3)*2 + 1;
% x = -1-halfw:halfw;
% dx = -halfw:halfw;
% kernel = diff(1 + 0.5*erf((x+0.5)/(eNoiseSigmakeV*sqrt(2))));
% 
% % Load PSF
% 
% fprintf('Loading PSF ... ')
% 
% if exist([design 'pixel_psf.mat'],'file')
%     
%     m_psf = matfile([design 'pixel_psf.mat']);
%     PSF = m_psf.PSF;
%     
% else
%     
%     m_PSF = matfile([design 'intermediatePsf.mat']);        
%     
%     PSF = (squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(1,1,:,:,:,:,:)) + ...
%         squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(1,2,:,:,:,:,:)) + ...
%         flip(squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(1,2,:,:,:,:,:)),2) + ...
%         squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,1,:,:,:,:,:)) + ...
%         flip(squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,1,:,:,:,:,:)),1) + ...
%         squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,2,:,:,:,:,:)) + ...
%         flip(squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,2,:,:,:,:,:)),1) + ...
%         flip(squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,2,:,:,:,:,:)),2) + ...
%         flip(flip(squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,2,:,:,:,:,:)),1),2)) / 9;
%     
%     m_psf = matfile([design 'pixel_psf.mat']);
%     m_psf.PSF = PSF;
%     
% end
% 
% 
% fprintf('done!\n')
% 
% % Add electronic noise to PSF
% 
% PSF_E = convn(PSF,reshape(kernel,[1,1,1,length(kernel)]),'same');

%MP: load psf and acf from files instead of re-generating
large_area_psf_file = sprintf('large_area_psf_sigma_%1.1f.mat',eNoiseSigmakeV);
m_large_area_psf = matfile([design large_area_psf_file]);
PSF_E=m_large_area_psf.PSF_E;
large_area_acf_file = sprintf('large_area_acf_sigma_%1.1f.mat',eNoiseSigmakeV);
m_large_area_acf = matfile([design large_area_acf_file]);
large_area_ACF_E = m_large_area_acf.ACF;


response_matrix = squeeze(sum(sum(PSF_E)))';

% Load attenuation

mu_water = m_settings.muWatermm_1' * 10;
mu_iodine = m_settings.muImm_1' * 10;
% mu_iodine = m_settings.muImm_1' * 10;

% Compute detected spectra along water line

N_spectra = 50;

A_water = linspace(0,50,N_spectra);

transmitted_spectra= zeros(N_inc,N_spectra);
detected_spectra = zeros(N_dep,N_spectra);
d_detected_spectra_A = zeros(N_dep,N_spectra,2);

for i = 1:N_spectra
    
    transmitted_spectrum = normalized_pre_patient_spectrum.*exp(-(A_water(i) * mu_water));
    transmitted_spectra(:,i)=transmitted_spectrum; %MP
    detected_spectra(:,i) = response_matrix * transmitted_spectrum;
    
    d_detected_spectra_A(:,i,1) = -response_matrix * ...
        (mu_water .* transmitted_spectrum);
    
    d_detected_spectra_A(:,i,2) = -response_matrix * ...
        (mu_iodine .* transmitted_spectrum);
    
    
end

% Compute Fisher information for ideal silicon detector

lambda = reshape(detected_spectra,[N_dep,N_spectra]);
dlambda_A1 = reshape(d_detected_spectra_A(:,:,1),[N_dep,N_spectra]);
dlambda_A2 = reshape(d_detected_spectra_A(:,:,2),[N_dep,N_spectra]);

X = lambda;
dX_A1 = dlambda_A1;
dX_A2 = dlambda_A2;
X_covariance_unbinned=permute(sum(large_area_ACF_E.*permute(transmitted_spectra,[1 3 4 2])),[2 3 4 1])+diag(enoise_diff_variance_unbinned); %MP

%MP: calculate binned covariance
%MP replaced to use matrix inversion, gives more correct treatment of covariance
fisher = zeros(2,2,N_spectra);
max_ind_spec_nonzero=find(min(X,[],2)>0,1,'last');
for i=1:N_spectra
    dX=[dX_A1(:,i),dX_A2(:,i)];
    fisher(:,:,i)=dX(1:max_ind_spec_nonzero,:)'*(X_covariance_unbinned(1:max_ind_spec_nonzero,1:max_ind_spec_nonzero,i)\dX(1:max_ind_spec_nonzero,:));
end
 
% fisher(1,1,:) = sum(dX_A1 .^2 ./ X,1,'omitnan');
% fisher(1,2,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
% fisher(2,1,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
% fisher(2,2,:) = sum(dX_A2 .^2 ./ X,1,'omitnan');

ideal_fisher = fisher;

% Optimize with coarse grid

step = 4;

q_opt = -inf;

tic

for bin_1 = 5
    for bin_2 = bin_1+step:step:80
        fprintf('b2=%d\n',bin_2) %MP
        for bin_3 = bin_2+step:step:80
            for bin_4 = bin_3+step:step:80
                for bin_5 = bin_4+step:step:80
                    for bin_6 = bin_5+step:step:80
                        for bin_7 = bin_6+step:step:80
                            for bin_8 = bin_7+step:step:80                                                                
                                T = [bin_1,bin_2,bin_3,bin_4,bin_5,bin_6,bin_7,bin_8,200]';                                
                                B_bin = -diff(D >= T);
                                                                
                                %MP: calculate binned covariance
                                X_covariance=zeros(length(T)-1,length(T)-1,N_spectra);
                                for j = 1:length(T)-1
                                    for k = 1:length(T)-1
                                        X_covariance(j,k,:) = sum(X_covariance_unbinned(B_bin(j,:)>0,B_bin(k,:)>0,:),[1 2]);
                                    end
                                end
                                X = B_bin * lambda;
                                dX_A1 = B_bin * dlambda_A1;
                                dX_A2 = B_bin * dlambda_A2;
                                %MP replaced to use matrix inversion, gives more correct treatment of covariance
                                for i=1:N_spectra
                                    dX=[dX_A1(:,i),dX_A2(:,i)];
                                    fisher(:,:,i)=dX'*(X_covariance(:,:,i)\dX);
                                end
%                                 fisher = zeros(2,2,N_spectra);
%                                 
%                                 fisher(1,1,:) = sum(dX_A1 .^2 ./ X,1,'omitnan');
%                                 fisher(1,2,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
%                                 fisher(2,1,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
%                                 fisher(2,2,:) = sum(dX_A2 .^2 ./ X,1,'omitnan');
                                fisher_b = squeeze(fisher(1,2,:))';
                                fisher_g = squeeze(fisher(1,1,:) - fisher(2,2,:))' / 2;
                                fisher_delta = sqrt(fisher_g .^2 + fisher_b.^2);
                                fisher_u = fisher_delta - fisher_g;
                                
                                fisher_v1 = reshape([fisher_b;fisher_u] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
                                fisher_v2 = reshape([-fisher_u;fisher_b] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
                                                                                                
                                q_1 = squeeze(sum(sum(fisher_v1.*fisher.*permute(fisher_v1,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v1.*ideal_fisher.*permute(fisher_v1,[2 1 3]))))';
                                q_2 = squeeze(sum(sum(fisher_v2.*fisher.*permute(fisher_v2,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v2.*ideal_fisher.*permute(fisher_v2,[2 1 3]))))';                                                                                              
                                
                                q = mean(mean([q_1;q_2]));
                                
                                if q > q_opt
                                    q_opt = q;
                                    opt_bins = T;                                    
                                end                                

                            end
                        end
                    end
                end
            end
        end
    end
end

toc
tic

% Optimize with fine grid

range = 2;

for bin_1 = 5
    for bin_2 = opt_bins(2)-range:opt_bins(2)+range
        for bin_3 = opt_bins(3)-range:opt_bins(3)+range
            for bin_4 = opt_bins(4)-range:opt_bins(4)+range
                for bin_5 = opt_bins(5)-range:opt_bins(5)+range
                    for bin_6 = opt_bins(6)-range:opt_bins(6)+range
                        for bin_7 = opt_bins(7)-range:opt_bins(7)+range
                            for bin_8 = opt_bins(8)-range:opt_bins(8)+range
                                
                                i = i+1;
                                
                                T = [bin_1,bin_2,bin_3,bin_4,bin_5,bin_6,bin_7,bin_8,200]';                                
                                B_bin = -diff(D >= T);
                                                                                                                           
                                %MP: calculate binned covariance
                                X_covariance=zeros(length(T)-1,length(T)-1,N_spectra);
                                for j = 1:length(T)-1
                                    for k = 1:length(T)-1
                                        X_covariance(j,k,:) = sum(X_covariance_unbinned(B_bin(j,:)>0,B_bin(k,:)>0,:),[1 2]);
                                    end
                                end
                                X = B_bin * lambda;
                                dX_A1 = B_bin * dlambda_A1;
                                dX_A2 = B_bin * dlambda_A2;
                                %MP replaced to use matrix inversion, gives more correct treatment of covariance
                                for i=1:N_spectra
                                    dX=[dX_A1(:,i),dX_A2(:,i)];
                                    fisher(:,:,i)=dX'*(X_covariance(:,:,i)\dX);
                                end
%                                 fisher = zeros(2,2,N_spectra);
%                                 
%                                 fisher(1,1,:) = sum(dX_A1 .^2 ./ X,1,'omitnan');
%                                 fisher(1,2,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
%                                 fisher(2,1,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
%                                 fisher(2,2,:) = sum(dX_A2 .^2 ./ X,1,'omitnan');
                                                                
                                fisher_b = squeeze(fisher(1,2,:))';
                                fisher_g = squeeze(fisher(1,1,:) - fisher(2,2,:))' / 2;
                                fisher_delta = sqrt(fisher_g .^2 + fisher_b.^2);
                                fisher_u = fisher_delta - fisher_g;
                                
                                fisher_v1 = reshape([fisher_b;fisher_u] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
                                fisher_v2 = reshape([-fisher_u;fisher_b] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
                                                                                                
                                q_1 = squeeze(sum(sum(fisher_v1.*fisher.*permute(fisher_v1,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v1.*ideal_fisher.*permute(fisher_v1,[2 1 3]))))';
                                q_2 = squeeze(sum(sum(fisher_v2.*fisher.*permute(fisher_v2,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v2.*ideal_fisher.*permute(fisher_v2,[2 1 3]))))';                                                                                              
                                
                                q = mean(mean([q_1;q_2]));
                                
                                if q > q_opt
                                    q_opt = q;
                                    opt_bins = T;                                    
                                end    
                                
                            end
                        end
                    end
                end
            end
        end
    end
end

toc


opt_bins'

T = opt_bins;
B_bin = -diff(D >= T);

X = B_bin * lambda;
dX_A1 = B_bin * dlambda_A1;
dX_A2 = B_bin * dlambda_A2;
X_covariance=zeros(length(T)-1,length(T)-1,N_spectra);
for j = 1:length(T)-1
    for k = 1:length(T)-1
        X_covariance(j,k,:) = sum(X_covariance_unbinned(B_bin(j,:)>0,B_bin(k,:)>0,:),[1 2]);
    end
end
%MP: calculate binned covariance
%MP replaced to use matrix inversion, gives more correct treatment of covariance
max_ind_spec_nonzero=find(min(X,[],2)>0,1,'last');
fisher = zeros(2,2,N_spectra);
for i=1:N_spectra
    dX=[dX_A1(:,i),dX_A2(:,i)];
    fisher(:,:,i)=dX'*(X_covariance(:,:,i)\dX);
end
% 
% fisher(1,1,:) = sum(dX_A1 .^2 ./ X,1,'omitnan');
% fisher(1,2,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
% fisher(2,1,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
% fisher(2,2,:) = sum(dX_A2 .^2 ./ X,1,'omitnan');

fisher_b = squeeze(fisher(1,2,:))';
fisher_g = squeeze(fisher(1,1,:) - fisher(2,2,:))' / 2;
fisher_delta = sqrt(fisher_g .^2 + fisher_b.^2);
fisher_u = fisher_delta - fisher_g;

fisher_v1 = reshape([fisher_b;fisher_u] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
fisher_v2 = reshape([-fisher_u;fisher_b] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);

q_1 = squeeze(sum(sum(fisher_v1.*fisher.*permute(fisher_v1,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v1.*ideal_fisher.*permute(fisher_v1,[2 1 3]))))';
q_2 = squeeze(sum(sum(fisher_v2.*fisher.*permute(fisher_v2,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v2.*ideal_fisher.*permute(fisher_v2,[2 1 3]))))';

figure(1)
clf
plot([q_1;q_2]')
xlabel('A_{water} [cm]')
legend('DQE_{density}','DQE_{spectral}')

% 4 bins

% Optimize with coarse grid

step = 2;

q_opt = -inf;

tic

for bin_1 = 5
    for bin_2 = bin_1+step:step:80
        for bin_3 = bin_2+step:step:80
            for bin_4 = bin_3+step:step:80
                                                                                             
                T = [bin_1,bin_2,bin_3,bin_4,200]';
                B_bin = -diff(D >= T);
                

                %MP: calculate binned covariance
                X_covariance=zeros(length(T)-1,length(T)-1,N_spectra);
                for j = 1:length(T)-1
                    for k = 1:length(T)-1
                        X_covariance(j,k,:) = sum(X_covariance_unbinned(B_bin(j,:)>0,B_bin(k,:)>0,:),[1 2]);
                    end
                end
                X = B_bin * lambda;
                dX_A1 = B_bin * dlambda_A1;
                dX_A2 = B_bin * dlambda_A2;
                %MP replaced to use matrix inversion, gives more correct treatment of covariance
                for i=1:N_spectra
                    dX=[dX_A1(:,i),dX_A2(:,i)];
                    fisher(:,:,i)=dX'*(X_covariance(:,:,i)\dX);
                end
                
                fisher_b = squeeze(fisher(1,2,:))';
                fisher_g = squeeze(fisher(1,1,:) - fisher(2,2,:))' / 2;
                fisher_delta = sqrt(fisher_g .^2 + fisher_b.^2);
                fisher_u = fisher_delta - fisher_g;
                
                fisher_v1 = reshape([fisher_b;fisher_u] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
                fisher_v2 = reshape([-fisher_u;fisher_b] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
                
                q_1 = squeeze(sum(sum(fisher_v1.*fisher.*permute(fisher_v1,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v1.*ideal_fisher.*permute(fisher_v1,[2 1 3]))))';
                q_2 = squeeze(sum(sum(fisher_v2.*fisher.*permute(fisher_v2,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v2.*ideal_fisher.*permute(fisher_v2,[2 1 3]))))';
                
                q = mean(mean([q_1;q_2]));
                
                if q > q_opt
                    q_opt = q;
                    opt_bins = T;
                end
                
            end
        end
    end
end


toc
tic

% Optimize with fine grid

range = 2;

for bin_1 = 5
    for bin_2 = opt_bins(2)-range:opt_bins(2)+range
        for bin_3 = opt_bins(3)-range:opt_bins(3)+range
            for bin_4 = opt_bins(4)-range:opt_bins(4)+range
                
                i = i+1;
                
                T = [bin_1,bin_2,bin_3,bin_4,200]';
                B_bin = -diff(D >= T);
                
                
                %MP: calculate binned covariance
                X_covariance=zeros(length(T)-1,length(T)-1,N_spectra);
                for j = 1:length(T)-1
                    for k = 1:length(T)-1
                        X_covariance(j,k,:) = sum(X_covariance_unbinned(B_bin(j,:)>0,B_bin(k,:)>0,:),[1 2]);
                    end
                end
                X = B_bin * lambda;
                dX_A1 = B_bin * dlambda_A1;
                dX_A2 = B_bin * dlambda_A2;
                %MP replaced to use matrix inversion, gives more correct treatment of covariance
                for i=1:N_spectra
                    dX=[dX_A1(:,i),dX_A2(:,i)];
                    fisher(:,:,i)=dX'*(X_covariance(:,:,i)\dX);
                end
                
                fisher_b = squeeze(fisher(1,2,:))';
                fisher_g = squeeze(fisher(1,1,:) - fisher(2,2,:))' / 2;
                fisher_delta = sqrt(fisher_g .^2 + fisher_b.^2);
                fisher_u = fisher_delta - fisher_g;
                
                fisher_v1 = reshape([fisher_b;fisher_u] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
                fisher_v2 = reshape([-fisher_u;fisher_b] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
                
                q_1 = squeeze(sum(sum(fisher_v1.*fisher.*permute(fisher_v1,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v1.*ideal_fisher.*permute(fisher_v1,[2 1 3]))))';
                q_2 = squeeze(sum(sum(fisher_v2.*fisher.*permute(fisher_v2,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v2.*ideal_fisher.*permute(fisher_v2,[2 1 3]))))';
                
                q = mean(mean([q_1;q_2]));
                
                if q > q_opt
                    q_opt = q;
                    opt_bins = T;
                end
                
            end
        end
    end
end

toc


opt_bins'

T = opt_bins;
B_bin = -diff(D >= T);

X = B_bin * lambda;
dX_A1 = B_bin * dlambda_A1;
dX_A2 = B_bin * dlambda_A2;
X_covariance=zeros(length(T)-1,length(T)-1,N_spectra);
for j = 1:length(T)-1
    for k = 1:length(T)-1
        X_covariance(j,k,:) = sum(X_covariance_unbinned(B_bin(j,:)>0,B_bin(k,:)>0,:),[1 2]);
    end
end
%MP: calculate binned covariance
%MP replaced to use matrix inversion, gives more correct treatment of covariance
fisher = zeros(2,2,N_spectra);
max_ind_spec_nonzero=find(min(X,[],2)>0,1,'last');
for i=1:N_spectra
    dX=[dX_A1(:,i),dX_A2(:,i)];
    fisher(:,:,i)=dX'*(X_covariance(:,:,i)\dX);
end
% 
% fisher = zeros(2,2,N_spectra);
% 
% fisher(1,1,:) = sum(dX_A1 .^2 ./ X,1,'omitnan');
% fisher(1,2,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
% fisher(2,1,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
% fisher(2,2,:) = sum(dX_A2 .^2 ./ X,1,'omitnan');

fisher_b = squeeze(fisher(1,2,:))';
fisher_g = squeeze(fisher(1,1,:) - fisher(2,2,:))' / 2;
fisher_delta = sqrt(fisher_g .^2 + fisher_b.^2);
fisher_u = fisher_delta - fisher_g;

fisher_v1 = reshape([fisher_b;fisher_u] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
fisher_v2 = reshape([-fisher_u;fisher_b] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);

q_1 = squeeze(sum(sum(fisher_v1.*fisher.*permute(fisher_v1,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v1.*ideal_fisher.*permute(fisher_v1,[2 1 3]))))';
q_2 = squeeze(sum(sum(fisher_v2.*fisher.*permute(fisher_v2,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v2.*ideal_fisher.*permute(fisher_v2,[2 1 3]))))';

figure(2)
clf
plot([q_1;q_2]')
xlabel('A_{water} [cm]')
legend('DQE_{density}','DQE_{spectral}')

% 2 bins

% Optimize with coarse grid

step = 2;

q_opt = -inf;

tic

for bin_1 = 5
    for bin_2 = bin_1+step:step:80        
                                                                                             
        T = [bin_1,bin_2,200]';
        B_bin = -diff(D >= T);
        
        %MP: calculate binned covariance
        X_covariance=zeros(length(T)-1,length(T)-1,N_spectra);
        for j = 1:length(T)-1
            for k = 1:length(T)-1
                X_covariance(j,k,:) = sum(X_covariance_unbinned(B_bin(j,:)>0,B_bin(k,:)>0,:),[1 2]);
            end
        end
        X = B_bin * lambda;
        dX_A1 = B_bin * dlambda_A1;
        dX_A2 = B_bin * dlambda_A2;
        %MP replaced to use matrix inversion, gives more correct treatment of covariance
        for i=1:N_spectra
            dX=[dX_A1(:,i),dX_A2(:,i)];
            fisher(:,:,i)=dX'*(X_covariance(:,:,i)\dX);
        end
        
        fisher_b = squeeze(fisher(1,2,:))';
        fisher_g = squeeze(fisher(1,1,:) - fisher(2,2,:))' / 2;
        fisher_delta = sqrt(fisher_g .^2 + fisher_b.^2);
        fisher_u = fisher_delta - fisher_g;
        
        fisher_v1 = reshape([fisher_b;fisher_u] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
        fisher_v2 = reshape([-fisher_u;fisher_b] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
        
        q_1 = squeeze(sum(sum(fisher_v1.*fisher.*permute(fisher_v1,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v1.*ideal_fisher.*permute(fisher_v1,[2 1 3]))))';
        q_2 = squeeze(sum(sum(fisher_v2.*fisher.*permute(fisher_v2,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v2.*ideal_fisher.*permute(fisher_v2,[2 1 3]))))';
        
        q = mean(mean([q_1;q_2]));
        
        if q > q_opt
            q_opt = q;
            opt_bins = T;
        end
        
    end
end


toc
tic

% Optimize with fine grid

range = 2;

for bin_1 = 5
    for bin_2 = opt_bins(2)-range:opt_bins(2)+range
               
        i = i+1;
        
        T = [bin_1,bin_2,200]';
        B_bin = -diff(D >= T);
        
        %MP: calculate binned covariance
        X_covariance=zeros(length(T)-1,length(T)-1,N_spectra);
        for j = 1:length(T)-1
            for k = 1:length(T)-1
                X_covariance(j,k,:) = sum(X_covariance_unbinned(B_bin(j,:)>0,B_bin(k,:)>0,:),[1 2]);
            end
        end
        X = B_bin * lambda;
        dX_A1 = B_bin * dlambda_A1;
        dX_A2 = B_bin * dlambda_A2;
        %MP replaced to use matrix inversion, gives more correct treatment of covariance
        for i=1:N_spectra
            dX=[dX_A1(:,i),dX_A2(:,i)];
            fisher(:,:,i)=dX'*(X_covariance(:,:,i)\dX);
        end

        
        fisher_b = squeeze(fisher(1,2,:))';
        fisher_g = squeeze(fisher(1,1,:) - fisher(2,2,:))' / 2;
        fisher_delta = sqrt(fisher_g .^2 + fisher_b.^2);
        fisher_u = fisher_delta - fisher_g;
        
        fisher_v1 = reshape([fisher_b;fisher_u] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
        fisher_v2 = reshape([-fisher_u;fisher_b] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
        
        q_1 = squeeze(sum(sum(fisher_v1.*fisher.*permute(fisher_v1,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v1.*ideal_fisher.*permute(fisher_v1,[2 1 3]))))';
        q_2 = squeeze(sum(sum(fisher_v2.*fisher.*permute(fisher_v2,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v2.*ideal_fisher.*permute(fisher_v2,[2 1 3]))))';
        
        q = mean(mean([q_1;q_2]));
        
        if q > q_opt
            q_opt = q;
            opt_bins = T;
        end
        
    end
end

toc


opt_bins'
% opt_bins = [5;33;200]

T = opt_bins;
B_bin = -diff(D >= T);

X = B_bin * lambda;
dX_A1 = B_bin * dlambda_A1;
dX_A2 = B_bin * dlambda_A2;
X_covariance=zeros(length(T)-1,length(T)-1,N_spectra);
for j = 1:length(T)-1
    for k = 1:length(T)-1
        X_covariance(j,k,:) = sum(X_covariance_unbinned(B_bin(j,:)>0,B_bin(k,:)>0,:),[1 2]);
    end
end
%MP: calculate binned covariance
%MP replaced to use matrix inversion, gives more correct treatment of covariance
fisher = zeros(2,2,N_spectra);
max_ind_spec_nonzero=find(min(X,[],2)>0,1,'last');
for i=1:N_spectra
    dX=[dX_A1(:,i),dX_A2(:,i)];
    fisher(:,:,i)=dX'*(X_covariance(:,:,i)\dX);
end
% fisher(1,1,:) = sum(dX_A1 .^2 ./ X,1,'omitnan');
% fisher(1,2,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
% fisher(2,1,:) = sum(dX_A1 .* dX_A2 ./ X,1,'omitnan');
% fisher(2,2,:) = sum(dX_A2 .^2 ./ X,1,'omitnan');

fisher_b = squeeze(fisher(1,2,:))';
fisher_g = squeeze(fisher(1,1,:) - fisher(2,2,:))' / 2;
fisher_delta = sqrt(fisher_g .^2 + fisher_b.^2);
fisher_u = fisher_delta - fisher_g;

fisher_v1 = reshape([fisher_b;fisher_u] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);
fisher_v2 = reshape([-fisher_u;fisher_b] ./ sqrt(2 * fisher_delta .* fisher_u),[2,1,N_spectra]);

q_1 = squeeze(sum(sum(fisher_v1.*fisher.*permute(fisher_v1,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v1.*ideal_fisher.*permute(fisher_v1,[2 1 3]))))';
q_2 = squeeze(sum(sum(fisher_v2.*fisher.*permute(fisher_v2,[2 1 3]))))' ./ squeeze(sum(sum(fisher_v2.*ideal_fisher.*permute(fisher_v2,[2 1 3]))))';

figure(3)
clf
plot([q_1;q_2]')
xlabel('A_{water} [cm]')
legend('DQE_{density}','DQE_{spectral}')













 