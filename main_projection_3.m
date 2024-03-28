clear
fprintf('Setup ... ')
fig_dir = '../../Paper/Latex/Revision III/';

print_figs = true;

design = '60mmActiveSi/'; 

dqeFilesToAverageR150mm={'si_cnr_dqe_%d_bin_%s_interpolation_ROIsize7.50_ROIspacing4.mat','si_cnr_dqe_%d_bin_%s_interpolation_ROIsize11.25_ROIspacing4.mat','si_cnr_dqe_%d_bin_%s_interpolation_ROIsize15.00_ROIspacing4.mat'}; %The one just named _interpolation.mat is ROIsize7.50_ROIspacing4.
dqeFilesToAverageR200mm={'si_cnr_dqe_%d_bin_%s_interpolation_ROIsize7.50_ROIspacing4.mat','si_cnr_dqe_%d_bin_%s_interpolation_ROIsize11.25_ROIspacing4.mat','si_cnr_dqe_%d_bin_%s_interpolation_ROIsize15.00_ROIspacing4.mat'}; %The one just named _interpolation.mat is ROIsize7.50_ROIspacing4.
% Load energy vectors and spectrum

m_settings = matfile([design 'SimulationSettings.mat']);

E = m_settings.eActVectorkeV;
D = m_settings.eDepVectorkeV;

post_patient_spectrum = m_settings.backgroundSpectrumAfterPatientmm_2';

fprintf('done!\n')

% Load PSF

fprintf('Loading PSF ... ')

if exist([design 'pixel_psf.mat'],'file')
    
    m_pixel_psf = matfile([design 'pixel_psf.mat']);
    PSF = m_pixel_psf.PSF;
    
else
    
    m_PSF = matfile([design 'intermediatePsf.mat']);        
    
    PSF = (squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(1,1,:,:,:,:,:)) + ...
        squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(1,2,:,:,:,:,:)) + ...
        flip(squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(1,2,:,:,:,:,:)),2) + ...
        squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,1,:,:,:,:,:)) + ...
        flip(squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,1,:,:,:,:,:)),1) + ...
        squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,2,:,:,:,:,:)) + ...
        flip(squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,2,:,:,:,:,:)),1) + ...
        flip(squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,2,:,:,:,:,:)),2) + ...
        flip(flip(squeeze(m_PSF.depositedSpectraInPixelsForBeamPositionsPerIncPh(2,2,:,:,:,:,:)),1),2)) / 9;
    
    m_pixel_psf = matfile([design 'pixel_psf.mat']);
    m_pixel_psf.PSF = PSF;
end


fprintf('done!\n')

% Compute Fano noise
    
e_Si = 3.6e-3;
F_Si = 0.115;

for eNoiseSigmakeV = [1.6] %3.2 %MP Note that it is the last value in this list that will be used for the plots
FanoSigmakeV = sqrt(F_Si * e_Si * D);


halfw = ceil(eNoiseSigmakeV*3)*2 + 1;
x = -1-halfw:halfw;
dx = -halfw:halfw;

kernels = zeros(length(D),length(dx));

for i = 1:length(D)
    sigma = sqrt(eNoiseSigmakeV^2 + FanoSigmakeV(i)^2);
    kernels(i,:) = diff(1 + 0.5*erf((x+0.5)/(sigma*sqrt(2))));
end

kernel_mu = mean(kernels);
kernel_Sigma = cov(kernels);

[V,Lambda] = eig(kernel_Sigma);
Lambda = diag(Lambda);
[~,idx] = sort(Lambda,'descend');

V = V(:,idx);
Lambda = Lambda(idx);

kernel_U = V(:,1);
kernel_T = (kernels - kernel_mu) * kernel_U;

% Add electronic noise and Fano noise to PSF

fprintf('Adding noise to PSF ... ')

% Add electronic noise to PSF

PSF_E = ...
    convn(PSF,reshape(kernel_mu,[1,1,1,length(dx)]),'same') + ...
    convn(PSF,reshape(kernel_U,[1,1,1,length(dx)]),'same') .* ...
    reshape(kernel_T,[1,1,1,length(D)]);

fprintf('done!\n')

%MP save the large area PSF
large_area_psf_file = sprintf('large_area_psf_sigma_%1.1f.mat',eNoiseSigmakeV);
if ~exist([design large_area_psf_file],'file')
    m_large_area_psf = matfile([design large_area_psf_file]);
    m_large_area_psf.PSF_E = PSF_E;
end
% Load acf

fprintf('Loading large area ACF ... ')

large_area_acf_file = sprintf('large_area_acf_sigma_%1.1f.mat',eNoiseSigmakeV);
large_area_1d_acf_file = sprintf('large_area_1d_acf_sigma_%1.1f.mat',eNoiseSigmakeV);

if exist([design large_area_acf_file],'file')
    
    m_large_area_acf = matfile([design large_area_acf_file]);
    large_area_ACF_E = m_large_area_acf.ACF;
    
else
    
    m_acf = matfile([design 'acf.mat']);

    large_area_ACF_E = zeros(length(E),length(D),length(D));
    large_area_1d_ACF_E = zeros(m_settings.nPixelsx,length(E),length(D),length(D)); %MP

    for i = 1:length(E)        

        i

        ACF_i = squeeze(m_acf.acfPerIncPh(:,:,:,:,i,:,:));

        large_area_ACF_E_i = 0;

        for j = 1:m_settings.nPixelsx
            large_area_ACF_E_ij = 0; %MP
            
            for k = 1:m_settings.nPixelsy

                % If central pixel:

                if j == (m_settings.nPixelsx+1)/2 && k == (m_settings.nPixelsy+1)/2

                    % Add 1-d electronic noise convolution of diagonal                    
                    
                    ACF_i_diag = diag(squeeze(ACF_i(j,k,:,:)));
                                                                                   
                    ACF_ijk_2 = diag(... #Modified by MP to be able to compute both 1d and low-frequency ACF
                            conv(ACF_i_diag, kernel_mu', 'same') + ...
                            conv(ACF_i_diag, kernel_U', 'same') .* kernel_T);                                            
                else
                    
                    % Add 2-d electronic noise convolution
                    
                    ACF_ijk = squeeze(ACF_i(j,k,:,:));
                    
                    ACF_ijk_1 = convn(ACF_ijk, kernel_mu', 'same') + ...  #Modified by MP to be able to compute both 1d and low-frequency ACF
                              convn(ACF_ijk, kernel_U', 'same') .* kernel_T;
                          
                    ACF_ijk_2 = convn(ACF_ijk_1, kernel_mu, 'same') + ... 
                              convn(ACF_ijk_1, kernel_U, 'same') .* kernel_T';                                                                                                                      
                end                                
                large_area_ACF_E_i = large_area_ACF_E_i + ACF_ijk_2;
                large_area_ACF_E_ij = large_area_ACF_E_ij + ACF_ijk_2;
            end
            large_area_1d_ACF_E(j,i,:,:) = large_area_ACF_E_ij;  %MP. j,i is the correct index order.
        end        

        large_area_ACF_E(i,:,:) = large_area_ACF_E_i;        

    end
    
    m_large_area_acf = matfile([design large_area_acf_file]);
    m_large_area_acf.ACF = large_area_ACF_E;
    m_acf = matfile([design large_area_1d_acf_file]); %Added by MP
    m_acf.ACF = large_area_1d_ACF_E;    %added by MP
end

fprintf('done!\n')

end

% Compute large area response

R_large_area = squeeze(sum(sum(PSF_E,1),2))';
R_large_row = squeeze(sum(sum(PSF_E(21,:,:,:),1),2))';
% Compute large area variance

R_var = zeros(length(D),length(E));

for i = 1:length(E)
    for j = 1:length(D) 
        R_var(j,i) = sum(sum(large_area_ACF_E(i,j:end,j:end),2),3);
    end
end

% Compute Compton/Photo thresholding

photo_idx = D' > E/2;
compton_idx = ~photo_idx;

R_primary = squeeze((PSF_E(21,21,:,:)))';
R_secondary = R_large_area - R_primary;


R_primary_photo = R_primary.*photo_idx;
R_primary_compton = R_primary.*compton_idx;


% Load unique number of photons

m = matfile([design 'acf.mat']);

R_unique = squeeze(m.uniquePhotonsAboveThresholdAfterChShPerIncPh)';

R_unique_E = ...
    [convn(-diff(R_unique),kernel_mu','same');zeros(1,length(E))] + ...
    [convn(-diff(R_unique),kernel_U','same');zeros(1,length(E))] .* kernel_T ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gray-scale and spectral dose efficiency simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_8_bin = [5     7    14    33    38    47    56    73   200]';
T_4_bin = [5    10    34    56   200]';
T_2_bin = [5    32   200]';

pre_patient_spectrum = m_settings.backgroundSpectrumBeforePatientmm_2';

%MP addition: calculate spectrum of electronic noise counts:
pre_patient_photons_per_mm2_and_mAs=1.8e6; %Persson et al Med. Phys. 43, 4398, 2016
mA=800;
measurement_time_seconds=0.25; %Should not affect the results, only gives a scale factor that is divided away in the end.
mAs=measurement_time_seconds*mA;
detector_pixel_area_mm2=0.25;
pre_patient_photons_per_mm2=pre_patient_photons_per_mm2_and_mAs*mAs;
% pre_patient_photons_per_pixel=pre_patient_photons_per_mm2*detector_pixel_area_mm2;
pre_patient_spectrum_rescaling_factor=pre_patient_photons_per_mm2/sum(pre_patient_spectrum); 
pre_patient_spectrum=pre_patient_spectrum_rescaling_factor*pre_patient_spectrum; %photons per measurement_time_seconds per mm2
post_patient_spectrum=pre_patient_spectrum_rescaling_factor*post_patient_spectrum; % photons per measurement_time_seconds per mm2
loaded_noise_sweep=load('./electronic_noise_simulations/corr1p6keV.mat');
electronic_noise_variance=measurement_time_seconds*loaded_noise_sweep.countsStd.^2./...
    (10^-9*loaded_noise_sweep.frametimens)/detector_pixel_area_mm2; %variance measured as equivalent quanta per measurement_time_seconds and mm2
thresholds_for_electronic_noise_variance_keV=loaded_noise_sweep.thresholdskeV;
% save(sprintf('electronic_noise_simulations/electronic_noise_data_%1.1fkeV.mat',eNoiseSigmakeV),'thresholds_for_electronic_noise_variance_keV','electronic_noise_variance',...
%     'detector_pixel_area_mm2','pre_patient_photons_per_mm2','pre_patient_photons_per_mm2_and_mAs','mAs')
% End of MP addition

normalized_post_patient_spectrum = post_patient_spectrum/sum(post_patient_spectrum);

mu_water = m_settings.muWatermm_1' * 10;
mu_iodine = m_settings.muDiluteImm_1' * 10;

R = R_large_area;

lowest_thresholds = 1:118;
A_waters = [30;40]; 

dqe_density_unbinned = nan(length(A_waters),length(lowest_thresholds));
dqe_spectral_unbinned = nan(length(A_waters),length(lowest_thresholds));

dqe_density_8_bin = nan(length(A_waters),length(lowest_thresholds));
dqe_spectral_8_bin = nan(length(A_waters),length(lowest_thresholds));

dqe_density_4_bin = nan(length(A_waters),length(lowest_thresholds));
dqe_spectral_4_bin = nan(length(A_waters),length(lowest_thresholds));

dqe_density_2_bin = nan(length(A_waters),length(lowest_thresholds));
dqe_spectral_2_bin = nan(length(A_waters),length(lowest_thresholds));

V_unbinned = zeros(2, 2, length(lowest_thresholds));
V_8_bin = zeros(2, 2, length(lowest_thresholds));
V_4_bin = zeros(2, 2, length(lowest_thresholds));
V_2_bin = zeros(2, 2, length(lowest_thresholds));

for i_task = 1:length(A_waters)
    for i_threshold = 1:length(lowest_thresholds)
        
        T_0 = lowest_thresholds(i_threshold)
        
        T = (T_0:200)';
        N_bins = length(T) - 1;
        
        % Compute bins
        
        B_unbinned = -diff(D >= T);
        
        B_8_bin = -diff(D >= max([T_0, T_8_bin(2:end)'],T_0)');
        B_4_bin = -diff(D >= max([T_0, T_4_bin(2:end)'],T_0)');
        B_2_bin = -diff(D >= max([T_0, T_2_bin(2:end)'],T_0)');
        
        % Compute PSFs
        
        R_unbinned = B_unbinned*R;
        R_8_bin = B_8_bin*R;
        R_4_bin = B_4_bin*R;
        R_2_bin = B_2_bin*R;
        
        %
        % Compute electronic noise variance contributed by each bin (added by MP)
        enoise_diff_variance_unbinned=-diff(interp1(thresholds_for_electronic_noise_variance_keV+0.5,electronic_noise_variance,max([T_0, T(2:end)'],T_0),'linear',0)); %Only works with 1 keV steps! (otherwise will be incorrectly normalized)
        enoise_diff_variance_8_bin=-diff(interp1(thresholds_for_electronic_noise_variance_keV+0.5,electronic_noise_variance,max([T_0, T_8_bin(2:end)'],T_0),'linear',0));
        enoise_diff_variance_4_bin=-diff(interp1(thresholds_for_electronic_noise_variance_keV+0.5,electronic_noise_variance,max([T_0, T_4_bin(2:end)'],T_0),'linear',0));
        enoise_diff_variance_2_bin=-diff(interp1(thresholds_for_electronic_noise_variance_keV+0.5,electronic_noise_variance,max([T_0, T_2_bin(2:end)'],T_0),'linear',0)); %increase threshold by 0.5 in order to match the viewpoint that a threeshold of for example 20 actually includes everything from 0.5 keV and up.

        % Compute ACFs
        
        ACF_unbinned = zeros(131,N_bins,N_bins);
        ACF_8_bin = zeros(131,8,8);
        ACF_4_bin = zeros(131,4,4);
        ACF_2_bin = zeros(131,2,2);
        
        % Unbinned ACF
        
        for j = 1:length(T)-1
            for k = 1:length(T)-1
                ACF_unbinned(:,j,k) = sum(large_area_ACF_E(:,find(B_unbinned(j,:)),find(B_unbinned(k,:))),[2 3]);
            end
        end
        
        % Binned ACFs
        
        for j = 1:8
            for k = 1:8
                ACF_8_bin(:,j,k) = sum(large_area_ACF_E(:,find(B_8_bin(j,:)),find(B_8_bin(k,:))),[2 3]);
            end
        end
        
        for j = 1:4
            for k = 1:4
                ACF_4_bin(:,j,k) = sum(large_area_ACF_E(:,find(B_4_bin(j,:)),find(B_4_bin(k,:))),[2 3]);
            end
        end
        
        for j = 1:2
            for k = 1:2
                ACF_2_bin(:,j,k) = sum(large_area_ACF_E(:,find(B_2_bin(j,:)),find(B_2_bin(k,:))),[2 3]);
            end
        end
        
        % Simulate task
        
        A_water = A_waters(i_task);
        A_iodine = 0;
        
        % Unbinned
        
        lambda_unbinned = R_unbinned*(pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        dlambda_unbinned_A1 = R_unbinned*(-mu_water.*pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        dlambda_unbinned_A2 = R_unbinned*(-mu_iodine.*pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        
        dlambda_unbinned = [dlambda_unbinned_A1';dlambda_unbinned_A2']';
        Sigma_unbinned = squeeze(sum(ACF_unbinned.*(pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)))))+diag(enoise_diff_variance_unbinned);
        
        % 8 bin
        
        lambda_8_bin = R_8_bin*(pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        dlambda_8_bin_A1 = R_8_bin*(-mu_water.*pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        dlambda_8_bin_A2 = R_8_bin*(-mu_iodine.*pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        
        dlambda_8_bin = [dlambda_8_bin_A1';dlambda_8_bin_A2']';
        Sigma_8_bin = squeeze(sum(ACF_8_bin.*(pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)))))...
            +diag(enoise_diff_variance_8_bin);
        
        % 4 bin
        
        lambda_4_bin = R_4_bin*(pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        dlambda_4_bin_A1 = R_4_bin*(-mu_water.*pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        dlambda_4_bin_A2 = R_4_bin*(-mu_iodine.*pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        
        dlambda_4_bin = [dlambda_4_bin_A1';dlambda_4_bin_A2']';
        Sigma_4_bin = squeeze(sum(ACF_4_bin.*(pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)))))...
            +diag(enoise_diff_variance_4_bin);
        
        % 2 bin
        
        lambda_2_bin = R_2_bin*(pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        dlambda_2_bin_A1 = R_2_bin*(-mu_water.*pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        dlambda_2_bin_A2 = R_2_bin*(-mu_iodine.*pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        
        dlambda_2_bin = [dlambda_2_bin_A1';dlambda_2_bin_A2']';
        Sigma_2_bin = squeeze(sum(ACF_2_bin.*(pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)))))...
            +diag(enoise_diff_variance_2_bin);
        
        % Ideal
        
        lambda_ideal = (pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        dlambda_A1_ideal = (-mu_water.*pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        dlambda_A2_ideal = (-mu_iodine.*pre_patient_spectrum.*exp(-(mu_water*A_water + mu_iodine*A_iodine)));
        
        dlambda_ideal = [dlambda_A1_ideal';dlambda_A2_ideal']';
        Sigma_ideal = diag(lambda_ideal);
        
        % Compute Fisher information
        
        idx = lambda_unbinned > 1e-12;
        Fisher_unbinned = dlambda_unbinned(idx,:)'*(Sigma_unbinned(idx,idx)\dlambda_unbinned(idx,:));
        
        idx = lambda_8_bin > 1e-12;
        Fisher_8_bin = dlambda_8_bin(idx,:)'*(Sigma_8_bin(idx,idx)\dlambda_8_bin(idx,:));
        
        idx = lambda_4_bin > 1e-12;
        Fisher_4_bin = dlambda_4_bin(idx,:)'*(Sigma_4_bin(idx,idx)\dlambda_4_bin(idx,:));
        
        idx = lambda_2_bin > 1e-12;
        Fisher_2_bin = dlambda_2_bin(idx,:)'*(Sigma_2_bin(idx,idx)\dlambda_2_bin(idx,:));
        
        idx = lambda_ideal > 0;
        Fisher_ideal = dlambda_ideal(idx,:)'*(Sigma_ideal(idx,idx)\dlambda_ideal(idx,:));
        
        % Compute CRLB
        
        CRLB_unbinned = Fisher_unbinned\eye(2);
        CRLB_8_bin = Fisher_8_bin\eye(2);
        CRLB_4_bin = Fisher_4_bin\eye(2);
        CRLB_2_bin = Fisher_2_bin\eye(2);
        CRLB_ideal = Fisher_ideal\eye(2);
        
        % Compute DQEs
        
        [V,Eig] = eig(CRLB_unbinned);
        [~,sort_idx] = sort(diag(Eig));
        V = V(:,sort_idx);
        V_unbinned(:, :, i_threshold) = V;
        
        DQE_unbinned = diag(V'*CRLB_ideal*V) ./ (diag(V'*CRLB_unbinned*V));
        DQE_unbinned = diag(V'*(CRLB_unbinned\V)) ./ diag(V'*(CRLB_ideal\V));
        
        [V,Eig] = eig(CRLB_8_bin);
        [~,sort_idx] = sort(diag(Eig));
        V = V(:,sort_idx);
        V_8_bin(:, :, i_threshold) = V;
                        
        DQE_8_bin = diag(V'*CRLB_ideal*V) ./ (diag(V'*CRLB_8_bin*V));      
        DQE_8_bin = diag(V'*(CRLB_8_bin\V)) ./ diag(V'*(CRLB_ideal\V));
        
        [V,Eig] = eig(CRLB_4_bin);
        [~,sort_idx] = sort(diag(Eig));
        V = V(:,sort_idx);
        V_4_bin(:, :, i_threshold) = V;
        
        DQE_4_bin = diag(V'*CRLB_ideal*V) ./ (diag(V'*CRLB_4_bin*V));
        DQE_4_bin = diag(V'*(CRLB_4_bin\V)) ./ diag(V'*(CRLB_ideal\V));
                
        [V,Eig] = eig(CRLB_2_bin);
        [~,sort_idx] = sort(diag(Eig));
        V = V(:,sort_idx);
        V_2_bin(:, :, i_threshold) = V;
        
        DQE_2_bin = diag(V'*CRLB_ideal*V) ./ (diag(V'*CRLB_2_bin*V));
        DQE_2_bin = diag(V'*(CRLB_2_bin\V)) ./ diag(V'*(CRLB_ideal\V));
        
        % Store
        
        dqe_density_unbinned(i_task,i_threshold) = DQE_unbinned(1);
        dqe_spectral_unbinned(i_task,i_threshold) = DQE_unbinned(2);
        
        if sum(sum(B_8_bin,2) > 0) >= 2            
            dqe_density_8_bin(i_task,i_threshold) = DQE_8_bin(1);
            dqe_spectral_8_bin(i_task,i_threshold) = DQE_8_bin(2);
        end
        
        if sum(sum(B_4_bin,2) > 0) >= 2        
            dqe_density_4_bin(i_task,i_threshold) = DQE_4_bin(1);
            dqe_spectral_4_bin(i_task,i_threshold) = DQE_4_bin(2);
        end
        
        if sum(sum(B_2_bin,2) > 0) >= 2
            dqe_density_2_bin(i_task,i_threshold) = DQE_2_bin(1);
            dqe_spectral_2_bin(i_task,i_threshold) = DQE_2_bin(2);
        end
        
        
    end
    
end

save('weight_matrices.mat', ...
    'V_unbinned',...
    'V_8_bin',...
    'V_4_bin',...
    'V_2_bin',...
    'lowest_thresholds')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

font_size = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 2a  Forward response and distribution plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Forward response and distribution plots

forward_energy = 70;

figure(1)
clf
hold on
plot(D,R_primary(:,E == forward_energy),'-','LineWidth',1.5)
plot(D,R_large_row(:,E == forward_energy),'--','LineWidth',1.5)
plot(D,R_large_area(:,E == forward_energy),'-.','LineWidth',1.5)

xlim([0 140])
ylim([0 0.075])

text(forward_energy - 13,1.05*R_large_area(D == forward_energy,E == forward_energy),'Photo peak','FontSize',12)
text(forward_energy/2 - 15,0.012,'Charge sharing tail','FontSize',12)
text(5,1.05*R_large_area(D == forward_energy/7,E == forward_energy)+0.01,'Compton hump','FontSize',12)
xlabel('E_{dep} [keV]')
ylabel('Interactions per inc. photon')
set(gca,'FontSize',font_size)
pause(.5)
legend('Single pixel','Large row','Large area')

if print_figs
    print('-depsc',[fig_dir 'Figure_2a.eps'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 2b  Distribution of events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
clf
hold on

plot(D,[zeros(19,1);normalized_post_patient_spectrum],'-','LineWidth',1.5)
plot(D,R_large_area*normalized_post_patient_spectrum,'-','LineWidth',1.5)
plot(D,R_unique_E*normalized_post_patient_spectrum,'-.','LineWidth',1.5)

xlim([0 140])
ylim([0 0.05])
ylabel('Interactions per inc. photon')
legend('Incident spectrum','Total deposited spectrum','Primary deposited spectrum')
xlabel('E_{dep} [keV]')
set(gca,'FontSize',font_size)

if print_figs
    print('-depsc',[fig_dir 'Figure_2b.eps'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 3a Intrinsic bin response bins 1-3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute binning matrix

B = -diff(D >= T_8_bin);

figure(3)
clf
hold on

plot(E,(B(1:3,:)*R_large_area),'LineWidth',1.5)

xlim([20 140])
ylim([0 1])
xlabel('E_{in} [keV]')

leg_strs = cell(8,1);
for bin = 1:8
    leg_strs{bin} = sprintf('Bin %d: %d-%d keV',bin,T_8_bin(bin),T_8_bin(bin+1));
end

leg_strs{8} = sprintf('Bin 8:   > %d keV',T_8_bin(bin));



legend(leg_strs{1:3})
ylabel('Interactions per inc. photon')
set(gca,'FontSize',font_size)

if print_figs
    print('-depsc',[fig_dir 'Figure_3a.eps'])    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 3b Intrinsic bin response bins 4-8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
clf
hold on

set(gca,'ColorOrderIndex',4)
plot(E,(B(4:8,:)*R_large_area),'LineWidth',1.5)

xlim([20 140])
ylim([0 1])
xlabel('E_{in} [keV]')
legend(leg_strs{4:end})
ylabel('Interactions per inc. photon')
set(gca,'FontSize',font_size)

if print_figs
    print('-depsc',[fig_dir 'Figure_3b.eps'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 3c Spectrum weighted bin response bins 1-3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
clf
hold on

plot(E,(B(1:3,:)*R_large_area*diag(normalized_post_patient_spectrum)),'LineWidth',1.5)

xlim([20 140])
ylim([0 0.011])
xlabel('E_{in} [keV]')
legend(leg_strs{1:3})
ylabel(sprintf('Interactions per energy\n and inc. photon $[\\rm{keV}^{-1}]$'),'interpreter','latex')
set(gca,'FontSize',font_size)

if print_figs
    print('-depsc',[fig_dir 'Figure_3c.eps'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 3d Spectrum weighted bin response bins 4-8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
clf
hold on

set(gca,'ColorOrderIndex',4)
plot(E,(B(4:8,:)*imgaussfilt(R_large_area,[1.5,0.001])*diag(normalized_post_patient_spectrum)),'LineWidth',1.5)

xlim([20 140])
ylim([0 0.011])
xlabel('E_{in} [keV]')
legend(leg_strs{4:8})
ylabel(sprintf('Interactions per energy\n and inc. photon $[\\rm{keV}^{-1}]$'),'interpreter','latex')
set(gca,'FontSize',font_size)

if print_figs
    print('-depsc',[fig_dir 'Figure_3d.eps'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 4ab Bin PSFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_psf = matfile([design 'psf.mat']);
addpath([design 'source']) 

X = squeeze(m_psf.weightedPsfmm_2wf(:,:,:,5,:));
X_bin = squeeze(sum(dimtools.insdim(B',[1 2]).*X,3));


m_settings = matfile([design 'SimulationSettings.mat']);
area = (m_settings.pixelPitchxmm) * (m_settings.pixelPitchymm);

m_fourier = matfile([design 'FourierMetrics.mat']);

freqSamplePointsxmm_1 = m_fourier.freqSamplePointsxmm_1;
freqSamplePointsymm_1 = m_fourier.freqSamplePointsymm_1;


squaredApertureTransform = ...
    (sinc(m_settings.psfSampleIntervalxmm*freqSamplePointsxmm_1).^2)'.*...
     sinc(m_settings.psfSampleIntervalymm*freqSamplePointsymm_1).^2; %Size nFrequenciesx, nFrequenciesy
 

unbinnedMaterialBasisTransferMatrixForCurrSpectrummm_2 = ...
    m_fourier.unbinnedMaterialBasisTransferMatrixmm_2;

unbinnedCrossSpectralDensityForCurrSpectrummm2 = ...
    m_fourier.unbinnedCrossSpectralDensityForCurrentSpectrummm2;

binnedMaterialBasisTransferMatrixmm_2 = ...
    binWithENoise(m_fourier.transferFunctionWeightedBySpectrummm_2,m_settings.eDepVectorkeV,T_8_bin(1:8)',1.6,5);

OTF = squeeze(abs(binnedMaterialBasisTransferMatrixmm_2./sqrt(squaredApertureTransform)));


MTF = OTF./OTF(m_fourier.frequencyMidindexx,m_fourier.frequencyMidindexy,:);

B = -diff(D >= T_8_bin);

figure(7)
clf
hold on
bin_psfs = squeeze(sum(X_bin,2)/sum(post_patient_spectrum)*area)/3;
normalized_bin_psfs = bin_psfs./sum(bin_psfs)*3;

plot(-61:61,normalized_bin_psfs(:,1:3),'LineWidth',1.5)
grid on

plot([-10.5 -10.5],[0 1.1],'k--','LineWidth',1)
plot([-7.5 -7.5],[0 1.1],'k--','LineWidth',1)
plot([-4.5 -4.5],[0 1.1],'k--','LineWidth',1)
plot([-1.5 -1.5],[0 1.1],'k--','LineWidth',1)
plot([1.5 1.5],[0 1.1],'k--','LineWidth',1)
plot([4.5 4.5],[0 1.1],'k--','LineWidth',1)
plot([7.5 7.5],[0 1.1],'k--','LineWidth',1)
plot([10.5 10.5],[0 1.1],'k--','LineWidth',1)


hleg=legend(leg_strs{1:3});
hleg.Position=hleg.Position+[0.03 0 0 0];
set(gca,'FontSize',font_size)
ylim([0 1.1])
xlim([-10.5 10.5])  
xticks([-12 -9 -6 -3 0 3 6 9 12])
xticklabels([-4 -3 -2 -1 0 1 2 3 4])
xlabel('Pixel number')
ylabel('Normalized bin PSF')

if print_figs
    print('-depsc',[fig_dir 'Figure_4a.eps'])
end


figure(8)
clf
hold on

set(gca,'ColorOrderIndex',4)
plot(-61:61,normalized_bin_psfs(:,4:8),'LineWidth',1.5)
grid on
% plot(-61:61,bin_psfs(:,8),'b')

plot([-10.5 -10.5],[0 1.1],'k--','LineWidth',1)
plot([-7.5 -7.5],[0 1.1],'k--','LineWidth',1)
plot([-4.5 -4.5],[0 1.1],'k--','LineWidth',1)
plot([-1.5 -1.5],[0 1.1],'k--','LineWidth',1)
plot([1.5 1.5],[0 1.1],'k--','LineWidth',1)
plot([4.5 4.5],[0 1.1],'k--','LineWidth',1)
plot([7.5 7.5],[0 1.1],'k--','LineWidth',1)
plot([10.5 10.5],[0 1.1],'k--','LineWidth',1)


hleg=legend(leg_strs{4:8});
hleg.Position=hleg.Position+[0.03 0 0 0];
set(gca,'FontSize',font_size)
ylim([0 1.1])
xlim([-10.5 10.5])  
xticks([-12 -9 -6 -3 0 3 6 9 12])
xticklabels([-4 -3 -2 -1 0 1 2 3 4])
xlabel('Pixel number')
ylabel('Normalized bin PSF')

if print_figs
    print('-depsc',[fig_dir 'Figure_4b.eps'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 4cd Bin MTFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_mtfs = fftshift(abs(fft(ifftshift(normalized_bin_psfs/3))))./sinc(m_settings.psfSampleIntervalxmm*freqSamplePointsxmm_1)';

figure(9)
clf
hold on
plot(m_fourier.freqSamplePointsxmm_1,bin_mtfs(:,1:3),'LineWidth',1.5)
plot(m_fourier.freqSamplePointsxmm_1,abs(sinc(m_settings.psfSampleIntervalxmm*freqSamplePointsxmm_1*3)),'k--','LineWidth',1.5)
xlabel('Spatial frequency [lp/mm]')
xlim([0 2.5])
ylim([0 1])
grid on
legend(leg_strs{1:3},'Square aperture')
set(gca,'FontSize',font_size)
ylabel('Bin MTFs')

if print_figs
    print('-depsc',[fig_dir 'Figure_4c.eps'])
end

figure(10)
clf
hold on
set(gca,'ColorOrderIndex',4)
plot(m_fourier.freqSamplePointsxmm_1,bin_mtfs(:,4:8),'LineWidth',1.5)
plot(m_fourier.freqSamplePointsxmm_1,abs(sinc(m_settings.psfSampleIntervalxmm*freqSamplePointsxmm_1*3)),'k--','LineWidth',1.5)

xlabel('Spatial frequency [lp/mm]')
xlim([0 2.5])
ylim([0 1])
grid on
legend(leg_strs{4:8},'Square aperture')
set(gca,'FontSize',font_size)
ylabel('Bin MTFs')

if print_figs
    print('-depsc',[fig_dir 'Figure_4d.eps'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 4ed Bin NPSs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist([design large_area_1d_acf_file],'file') %MP
    
    m_acf = matfile([design large_area_1d_acf_file]); %MP
    large_area_1d_ACF_E = m_acf.ACF;
    
end

B = -diff(D >= T_8_bin);

row_PSF = squeeze(sum(PSF_E,2));

row_bin_PSF = zeros(41,length(E),8);

row_bin_ACF = zeros(41,length(E),8,8);

N_bins = size(B,1);

for j = 1:N_bins
            
    row_bin_PSF(:,:,j) = sum(row_PSF(:,:,B(j,:) == 1),3);
    
    for k = 1:N_bins
        
        row_bin_ACF(:,:,j,k) = sum(large_area_1d_ACF_E(:,:,B(j,:) == 1, B(k,:) == 1),[3 4]);
        
    end
end

weighted_row_bin_PSF = squeeze(sum(row_bin_PSF.* post_patient_spectrum',2))*area;

differential_noise_spectrum_for_nps_plot=-diff(interp1(thresholds_for_electronic_noise_variance_keV+0.5,electronic_noise_variance,T_8_bin,'linear',0)); %increase threshold by 0.5 in order to match the viewpoint that a threeshold of for example 20 actually includes everything from 0.5 keV and up.

weighted_row_bin_ACF = squeeze(sum(row_bin_ACF.* post_patient_spectrum',2))*area;
weighted_row_bin_ACF(round((size(weighted_row_bin_ACF,1)+1)/2),:,:)=...
    weighted_row_bin_ACF(round((size(weighted_row_bin_ACF,1)+1)/2),:,:)+...
    permute(diag(differential_noise_spectrum_for_nps_plot),[3 1 2])*area; %Add electr noise counts, rescaled to be counts (per measurement) per pixel instead of per area 
bin_npss = fftshift(abs(fft(ifftshift(weighted_row_bin_ACF)))) ./ sum(weighted_row_bin_PSF);
NPSFreqSamplePointsSpacingmm_1=1/(m_settings.pixelPitchxmm)/length(weighted_row_bin_ACF);
NPSFreqSamplePointsmm_1=NPSFreqSamplePointsSpacingmm_1*(-floor(length(weighted_row_bin_ACF)/2):floor(length(weighted_row_bin_ACF)/2));
% NPSFreqSamplePointsmm_1=NPSFreqSamplePointsmm_1(2:end);
figure(11)
clf
hold on
for bin = 1:3
    plot(NPSFreqSamplePointsmm_1,bin_npss(:,bin,bin),'LineWidth',1.5)
end
xlabel('Spatial frequency [lp/mm]')
xlim([0 1])
ylim([0.8 2.3])
grid on
legend(leg_strs{1:3},'Location','East')
set(gca,'FontSize',font_size)
ylabel('Normalized bin NPS')
if print_figs
    print('-depsc',[fig_dir 'Figure_4e.eps'])
end

figure(12)
clf
hold on
set(gca,'ColorOrderIndex',4)
for bin = 4:8
    plot(NPSFreqSamplePointsmm_1,bin_npss(:,bin,bin),'LineWidth',1.5)
end
xlabel('Spatial frequency [lp/mm]')
xlim([0 1])
ylim([0.8 1.2])
grid on
legend(leg_strs{4:8})
set(gca,'FontSize',font_size)
ylabel('Normalized bin NPS')
if print_figs
    print('-depsc',[fig_dir 'Figure_4f.eps'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 5 SPR DQE-Penalty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MP note: electronic noise counts are not included in this particular plot - explain this in the paper.
unique_spectrum = cumsum(R_unique_E*normalized_post_patient_spectrum,'reverse');

total_spectrum = cumsum(R_large_area*normalized_post_patient_spectrum,'reverse');

var_spectrum = R_var*normalized_post_patient_spectrum;

secondary_spectrum = (total_spectrum - unique_spectrum);
rho = secondary_spectrum./total_spectrum;

spr = secondary_spectrum ./ unique_spectrum;

figure(13)
clf
hold on
grid on
plot(D,unique_spectrum,'LineWidth',1.5)
plot(D,spr./(1 + spr),'-v','LineWidth',1.5)
plot(D,spr,'-o','LineWidth',1.5)
plot(D,(1 + spr)./(1 + 2*spr),'--','LineWidth',1.5)
plot(D,unique_spectrum.*(1 + spr)./(1 + 2*spr),'-.','LineWidth',1.5)
plot(D,total_spectrum.^2./var_spectrum,':','LineWidth',1.5)

xlabel('Lowest threshold energy [keV]')
xlim([0 35])

leg = legend('Primary counts','\rho (Scattering probability)','SPR','Scatter DQE penalty','DQE: Geometric-Poisson','DQE: Monte Carlo','Location','NorthEast');

rect = [0.59, 0.56, .25, .25];
set(leg, 'Position', rect)

set(gca,'FontSize',font_size)

if print_figs
    print('-depsc',[fig_dir 'Figure_5.eps'])
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figs 6. Gray-scale & spectral dose-efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% print_figs = false
x = lowest_thresholds; %Changed by MP, started on 0 before

figure(14)
clf
hold on
grid on
plot(x,dqe_density_unbinned(1,:),'-','LineWidth',1.5)
plot(x,dqe_density_8_bin(1,:),':','LineWidth',1.5)
plot(x,dqe_density_4_bin(1,:),'-.','LineWidth',1.5)
plot(x,dqe_density_2_bin(1,:),'--','LineWidth',1.5)

xlabel('Lowest threshold energy [keV]')

leg = legend(...
    '120-bin DQE^{density}_{projection}',...
    '8-bin     DQE^{density}_{projection}',...
    '4-bin     DQE^{density}_{projection}',...
    '2-bin     DQE^{density}_{projection}',...    
    'Location','NorthEast');

text(19.1,0.45,'30 cm water background','FontSize',12)

set(gca,'FontSize',font_size)

xlim([0 35])
ylim([0.3 0.7])

if print_figs
    print('-depsc',[fig_dir 'Figure_6a.eps'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(15)
clf
hold on
grid on
plot(x,dqe_density_unbinned(2,:),'-','LineWidth',1.5)
plot(x,dqe_density_8_bin(2,:),':','LineWidth',1.5)
plot(x,dqe_density_4_bin(2,:),'-.','LineWidth',1.5)
plot(x,dqe_density_2_bin(2,:),'--','LineWidth',1.5)

xlabel('Lowest threshold energy [keV]')
set(gca,'FontSize',font_size)

xlim([0 35])
ylim([0.3 0.7])

leg = legend(...
    '120-bin DQE^{density}_{projection}',...
    '8-bin     DQE^{density}_{projection}',...
    '4-bin     DQE^{density}_{projection}',...
    '2-bin     DQE^{density}_{projection}',...    
    'Location','NorthEast');

text(19.1,0.45,'40 cm water background','FontSize',12)

if print_figs
    print('-depsc',[fig_dir 'Figure_6b.eps'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(16)
clf
hold on
grid on
plot(x,dqe_spectral_unbinned(1,:),'-','LineWidth',1.5)
plot(x,dqe_spectral_8_bin(1,:),':','LineWidth',1.5)
plot(x,dqe_spectral_4_bin(1,:),'-.','LineWidth',1.5)
plot(x,dqe_spectral_2_bin(1,:),'--','LineWidth',1.5)


xlim([0 35])
ylim([0 0.45])

set(gca,'FontSize',font_size)
xlabel('Lowest threshold energy [keV]')

leg = legend(...
    '120-bin DQE^{spectral}_{projection}',...
    '8-bin     DQE^{spectral}_{projection}',...
    '4-bin     DQE^{spectral}_{projection}',...
    '2-bin     DQE^{spectral}_{projection}',...    
    'Location','NorthEast');
leg.Position=leg.Position+[0.05 0.03 0 0];

text(2.5,0.35,'30 cm water background','FontSize',12)

if print_figs
    print('-depsc',[fig_dir 'Figure_6c.eps'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(17)
clf
hold on
grid on
plot(x,dqe_spectral_unbinned(2,:),'-','LineWidth',1.5)
plot(x,dqe_spectral_8_bin(2,:),':','LineWidth',1.5)
plot(x,dqe_spectral_4_bin(2,:),'-.','LineWidth',1.5)
plot(x,dqe_spectral_2_bin(2,:),'--','LineWidth',1.5)

xlim([0 35])
ylim([0 0.45])

set(gca,'FontSize',font_size)
xlabel('Lowest threshold energy [keV]')

leg = legend(...
    '120-bin DQE^{spectral}_{projection}',...
    '8-bin     DQE^{spectral}_{projection}',...
    '4-bin     DQE^{spectral}_{projection}',...
    '2-bin     DQE^{spectral}_{projection}',...    
    'Location','NorthEast');
leg.Position=leg.Position+[0.05 0.03 0 0];

text(2.5,0.35,'40 cm water background','FontSize',12)

if print_figs
    print('-depsc',[fig_dir 'Figure_6d.eps'])
end

%%

interpolation = 'linear'

dqe_density_8_bin_image = nan(2,length(lowest_thresholds));
dqe_density_4_bin_image = nan(2,length(lowest_thresholds));
dqe_density_2_bin_image = nan(2,length(lowest_thresholds));

dqe_spectral_8_bin_image = nan(2,length(lowest_thresholds));
dqe_spectral_4_bin_image = nan(2,length(lowest_thresholds));
dqe_spectral_2_bin_image = nan(2,length(lowest_thresholds));

averageDqes=[];
dqes=[];
for PHANTOM_RADIUS_MM = 150
    
    for number_of_bins = [8 4 2]
        dqe =0;
        for dqe_roi_config_no=1:length(dqeFilesToAverageR150mm)
            m_dqe = matfile([design sprintf('sinogram_data_%d_cm',2*PHANTOM_RADIUS_MM/10) filesep sprintf(dqeFilesToAverageR150mm{dqe_roi_config_no},number_of_bins,interpolation)]);
            dqes(dqe_roi_config_no,round(log(number_of_bins)/log(2)),1:size(m_dqe.dqe,1),:)=m_dqe.dqe;
            dqe = dqe +m_dqe.dqe;
        end
        dqe = dqe / length(dqeFilesToAverageR150mm);
        averageDqes(round(log(number_of_bins)/log(2)),1:size(dqe,1),:)=dqe;
        switch number_of_bins
            case 8
                dqe_density_8_bin_image(1,1:length(dqe)) = dqe(:,1);
                dqe_spectral_8_bin_image(1,1:length(dqe)) = dqe(:,2);
                
            case 4
                dqe_density_4_bin_image(1,1:length(dqe)) = dqe(:,1);
                dqe_spectral_4_bin_image(1,1:length(dqe)) = dqe(:,2);
                
            case 2
                dqe_density_2_bin_image(1,1:length(dqe)) = dqe(:,1);
                dqe_spectral_2_bin_image(1,1:length(dqe)) = dqe(:,2);
        end
    end    
   
end


for PHANTOM_RADIUS_MM = 200
    
    for number_of_bins = [8 4 2]
        dqe =0;
        for dqe_roi_config_no=1:length(dqeFilesToAverageR200mm)
            m_dqe = matfile([design sprintf('sinogram_data_%d_cm',2*PHANTOM_RADIUS_MM/10) filesep sprintf(dqeFilesToAverageR200mm{dqe_roi_config_no},number_of_bins,interpolation)]);
            dqe = dqe +m_dqe.dqe;
        end
        dqe = dqe / length(dqeFilesToAverageR200mm);        
        switch number_of_bins
            case 8
                dqe_density_8_bin_image(2,1:length(dqe)) = dqe(:,1);
                dqe_spectral_8_bin_image(2,1:length(dqe)) = dqe(:,2);
                
            case 4
                dqe_density_4_bin_image(2,1:length(dqe)) = dqe(:,1);
                dqe_spectral_4_bin_image(2,1:length(dqe)) = dqe(:,2);
                
            case 2
                dqe_density_2_bin_image(2,1:length(dqe)) = dqe(:,1);
                dqe_spectral_2_bin_image(2,1:length(dqe)) = dqe(:,2);
        end
    end    
   
end
%Maximum deviation between one of the DQEs and the average
max(abs(squeeze(1-dqes(:,:,:,1)./permute(averageDqes(:,:,1),[4 1 2 3]))),[],[1,2,3,4])
%Spectral
max(abs(squeeze(1-dqes(:,:,:,2)./permute(averageDqes(:,:,2),[4 1 2 3]))),[],[1,2,3,4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = lowest_thresholds; %Changed by MP, started on 0 before

figure(18)
clf
hold on
grid on
plot(x,dqe_density_unbinned(1,:),'-','LineWidth',1.5)
plot(x,dqe_density_8_bin_image(1,:),'o','LineWidth',1.5)
plot(x,dqe_density_4_bin_image(1,:),'x','LineWidth',1.5)
plot(x,dqe_density_2_bin_image(1,:),'s','LineWidth',1.5)
set(gca,'ColorOrderIndex',2)
plot(x,dqe_density_8_bin(1,:),':','LineWidth',1.5)
plot(x,dqe_density_4_bin(1,:),'-.','LineWidth',1.5)
plot(x,dqe_density_2_bin(1,:),'--','LineWidth',1.5)

xlabel('Lowest threshold energy [keV]')

leg = legend(...
    '120-bin DQE^{density}_{projection}',...
    '8-bin  DQE^{density}_{image}',...
    '4-bin  DQE^{density}_{image}',...
    '2-bin  DQE^{density}_{image}',...    
    'Location','NorthEast');

text(19.1,0.45,'30 cm water background','FontSize',12)

set(gca,'FontSize',font_size)

xlim([0 35])
ylim([0.3 0.7])

if print_figs
    print('-depsc',[fig_dir 'Figure_7a.eps'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(19)
clf
hold on
grid on
plot(x,dqe_density_unbinned(2,:),'-','LineWidth',1.5)
plot(x,dqe_density_8_bin_image(2,:),'o','LineWidth',1.5)
plot(x,dqe_density_4_bin_image(2,:),'x','LineWidth',1.5)
plot(x,dqe_density_2_bin_image(2,:),'s','LineWidth',1.5)
set(gca,'ColorOrderIndex',2)
plot(x,dqe_density_8_bin(2,:),':','LineWidth',1.5)
plot(x,dqe_density_4_bin(2,:),'-.','LineWidth',1.5)
plot(x,dqe_density_2_bin(2,:),'--','LineWidth',1.5)

xlabel('Lowest threshold energy [keV]')

leg = legend(...
    '120-bin DQE^{density}_{projection}',...
    '8-bin  DQE^{density}_{image}',...
    '4-bin  DQE^{density}_{image}',...
    '2-bin  DQE^{density}_{image}',...    
    'Location','NorthEast');

xlabel('Lowest threshold energy [keV]')
set(gca,'FontSize',font_size)

xlim([0 35])
ylim([0.3 0.7])

leg = legend(...
    '120-bin DQE^{density}_{projection}',...
    '8-bin  DQE^{density}_{image}',...
    '4-bin  DQE^{density}_{image}',...
    '2-bin  DQE^{density}_{image}',...    
    'Location','NorthEast')

text(19.1,0.45,'40 cm water background','FontSize',12)

if print_figs
    print('-depsc',[fig_dir 'Figure_7b.eps'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(20)
clf
hold on
grid on
plot(x,dqe_spectral_unbinned(1,:),'-','LineWidth',1.5)
plot(x,dqe_spectral_8_bin_image(1,:),'o','LineWidth',1.5)
plot(x,dqe_spectral_4_bin_image(1,:),'x','LineWidth',1.5)
plot(x,dqe_spectral_2_bin_image(1,:),'s','LineWidth',1.5)
set(gca,'ColorOrderIndex',2)
plot(x,dqe_spectral_8_bin(1,:),':','LineWidth',1.5)
plot(x,dqe_spectral_4_bin(1,:),'-.','LineWidth',1.5)
plot(x,dqe_spectral_2_bin(1,:),'--','LineWidth',1.5)


xlim([0 35])
ylim([0 0.45])

set(gca,'FontSize',font_size)
xlabel('Lowest threshold energy [keV]')

leg = legend(...
    '120-bin DQE^{spectral}_{projection}',...
    '8-bin  DQE^{spectral}_{image}',...
    '4-bin  DQE^{spectral}_{image}',...
    '2-bin  DQE^{spectral}_{image}',...    
    'Location','NorthEast')
leg.Position=leg.Position+[0.05 0 0 0];
text(2.5,0.35,'30 cm water background','FontSize',12)

if print_figs
    print('-depsc',[fig_dir 'Figure_7c.eps'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(21)
clf
hold on
grid on
plot(x,dqe_spectral_unbinned(2,:),'-','LineWidth',1.5)
plot(x,dqe_spectral_8_bin_image(2,:),'o','LineWidth',1.5)
plot(x,dqe_spectral_4_bin_image(2,:),'x','LineWidth',1.5)
plot(x,dqe_spectral_2_bin_image(2,:),'s','LineWidth',1.5)
set(gca,'ColorOrderIndex',2)
plot(x,dqe_spectral_8_bin(2,:),':','LineWidth',1.5)
plot(x,dqe_spectral_4_bin(2,:),'-.','LineWidth',1.5)
plot(x,dqe_spectral_2_bin(2,:),'--','LineWidth',1.5)

xlim([0 35])
ylim([0 0.45])

set(gca,'FontSize',font_size)
xlabel('Lowest threshold energy [keV]')

leg = legend(...
    '120-bin DQE^{spectral}_{projection}',...
    '8-bin  DQE^{spectral}_{image}',...
    '4-bin  DQE^{spectral}_{image}',...
    '2-bin  DQE^{spectral}_{image}',...    
    'Location','NorthEast')
leg.Position=leg.Position+[0.05 0 0 0];
text(2.5,0.35,'40 cm water background','FontSize',12)

if print_figs
    print('-depsc',[fig_dir 'Figure_7d.eps'])
end

%%

DECOMPOSITION = 'WATER_IODINE';

m_settings = matfile([design 'SimulationSettings.mat']);


basis_pairs = strsplit(DECOMPOSITION,'_');

E = m_settings.eActVectorkeV;
N_inc = length(E);
F = zeros(N_inc,2);

for i = 1:2
    switch basis_pairs{i}
        case 'WATER'
            F(:,i) = m_settings.muWatermm_1' * 10;
        case 'BONE'
            F(:,i) = m_settings.muBonemm_1' * 10;
        case 'IODINE'
            F(:,i) = m_settings.muDiluteImm_1'' * 10;
    end
end

X = F ./ F(:,1) * 1000;


m_cov = matfile('ensemble_covariances_ROIsize11.25_ROIspacing4.mat')

real_ensemble_covariance_30_cm_5_keV = m_cov.real_ensemble_covariance_30_cm_5_keV
real_ensemble_covariance_30_cm_35_keV = m_cov.real_ensemble_covariance_30_cm_35_keV
real_ensemble_covariance_40_cm_5_keV = m_cov.real_ensemble_covariance_40_cm_5_keV
real_ensemble_covariance_40_cm_35_keV = m_cov.real_ensemble_covariance_40_cm_35_keV

variance_30_cm_5_keV = zeros(N_inc,1);
variance_30_cm_35_keV = zeros(N_inc,1);
variance_40_cm_5_keV = zeros(N_inc,1);
variance_40_cm_35_keV = zeros(N_inc,1);

for i = 1:N_inc
    
    V = X(i,:)';
    
    variance_30_cm_5_keV(i) = (V'*real_ensemble_covariance_30_cm_5_keV*V);
    variance_30_cm_35_keV(i) = (V'*real_ensemble_covariance_30_cm_35_keV*V);
    variance_40_cm_5_keV(i) = (V'*real_ensemble_covariance_40_cm_5_keV*V);
    variance_40_cm_35_keV(i) = (V'*real_ensemble_covariance_40_cm_35_keV*V);
    
end


figure(22)
clf
semilogy(E,(variance_30_cm_5_keV),'-','LineWidth',1.5)
hold on
semilogy(E,(variance_30_cm_35_keV),'-.','LineWidth',1.5)
semilogy(E,(variance_40_cm_5_keV),':','LineWidth',1.5)
semilogy(E,(variance_40_cm_35_keV),'--','LineWidth',1.5)
grid on
xlim([35 120])
ylim([0.02 100]) %[0.1 500] for 7.5 mm ROIs
xlabel('Synthetic monoenergy [keV]')
ylabel('ROI ensemble variance [$\rm{HU}^2$]','interpreter','latex')
set(gca,'FontSize',font_size)

hl=legend(...
    '30 cm water, 8-bin with $T_1$ = 5 keV',... %Subscript not working due to bug in matlab
    '30 cm water, 8-bin with $T_1$ = 35 keV',...
    '40 cm water, 8-bin with $T_1$ = 5 keV',...
    '40 cm water, 8-bin with $T_1$ = 35 keV')
set(hl, 'Interpreter','latex')
if print_figs
    print('-depsc',[fig_dir 'Figure_8a.eps'])
end

figure(23)
clf
hold on
plot(E,variance_30_cm_5_keV ./ variance_30_cm_35_keV,'-', 'LineWidth',1.5)
plot(E,variance_40_cm_5_keV ./ variance_40_cm_35_keV,'-.', 'LineWidth',1.5)
grid on
xlim([35 120])
ylim([0 1.5])

xlabel('Synthetic monoenergy [keV]')
ylabel('ROI ensemble variance ratio')
set(gca,'FontSize',font_size)

hl=legend(...
    '30 cm water, 8-bin with $T_1$ = 5 keV / 35 keV',...
    '40 cm water, 8-bin with $T_1$ = 5 keV / 35 keV')
set(hl, 'Interpreter','latex')
if print_figs
    print('-depsc',[fig_dir 'Figure_8b.eps'])
end


%%
1
max(max(1 - 1.02*dqe_spectral_8_bin_image ./ dqe_spectral_8_bin))
max(max(1 - dqe_density_8_bin_image ./ dqe_density_8_bin))

max(max(1 - 1.02*dqe_spectral_4_bin_image ./ dqe_spectral_4_bin))
max(max(1 - dqe_density_4_bin_image ./ dqe_density_4_bin))

% max(max(1 - dqe_spectral_2_bin_image ./ dqe_spectral_2_bin))
% max(max(1 - dqe_density_2_bin_image ./ dqe_density_2_bin))


% dqe_gray(1) / dqe_gray(31)
% dqe_spectral(1) / dqe_spectral(31)
% 
% load spectral_image_dqe
% 
% dqe(6,:) ./ dqe(31,:)
% 
% image_dqe = dqe
% image_dqe_gray = image_dqe(:,1);
% image_dqe_spectral = image_dqe(:,2);
% 
% %%
% 
% agreement = (image_dqe_gray - dqe_gray_8(1:41)')./(dqe_gray_8(1:41)');
% agreement = agreement(1:31);
% max(abs(agreement))
% 
% agreement = (image_dqe_spectral - dqe_spectral_8(1:41)')./(dqe_spectral_8(1:41)');
% agreement = agreement(1:31);
% max(abs(agreement))

plot(x,dqe_spectral_unbinned(2,:),'-','LineWidth',1.5)
plot(x,dqe_spectral_8_bin_image(2,:),'o','LineWidth',1.5)
plot(x,dqe_spectral_4_bin_image(2,:),'x','LineWidth',1.5)
plot(x,dqe_spectral_2_bin_image(2,:),'s','LineWidth',1.5)
set(gca,'ColorOrderIndex',2)
plot(x,dqe_spectral_8_bin(2,:),':','LineWidth',1.5)
plot(x,dqe_spectral_4_bin(2,:),'-.','LineWidth',1.5)
plot(x,dqe_spectral_2_bin(2,:),'--','LineWidth',1.5)

%% MP additions:
% Maximum relative error (without sign) between image-domain and projection-domain metrics
disp('density')
max(abs(((1-dqe_density_2_bin_image ./ dqe_density_2_bin))'))
max(abs(((1-dqe_density_4_bin_image ./ dqe_density_4_bin))'))
max(abs(((1-dqe_density_8_bin_image ./ dqe_density_8_bin))'))
 
disp('spectral')
max(abs(((1-dqe_spectral_2_bin_image ./ dqe_spectral_2_bin))'))
max(abs(((1-dqe_spectral_4_bin_image ./ dqe_spectral_4_bin))'))
max(abs(((1-dqe_spectral_8_bin_image ./ dqe_spectral_8_bin))'))
 