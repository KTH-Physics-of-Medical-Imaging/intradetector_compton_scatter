clear

% Reconstruction parameters

filter_name = 'hann';
frequency_scaling = 1;
interpolation = 'linear';

print_figs = false;

design = '60mmActiveSi/';

% Load energy vectors and spectrum

m_settings = matfile([design 'SimulationSettings.mat']);

real_ensemble_covariance_30_cm_5_keV = nan;
real_ensemble_covariance_30_cm_35_keV = nan;
real_ensemble_covariance_40_cm_5_keV = nan;
real_ensemble_covariance_40_cm_35_keV = nan;

N_ensemble_realizations = 50;

for PHANTOM_RADIUS_MM = [150 200]
    
    % Image parameters
    
    N_pixels = 1024;
    FOV_mm = 420;
    image_pixel_size_mm = FOV_mm/N_pixels;
    
    % Fan beam parameters
    
    N_views = 2000;
    
    SID_mm = 500;
    SDD_mm = 1000;
    fan_beam_D = SID_mm / image_pixel_size_mm;
    
    fan_rotation_increment_deg = 360/N_views;
    fan_sensor_spacing_deg = atan(m_settings.pixelPitchxmm/SDD_mm)*180/pi;
    
    x = ((1:N_pixels) - (N_pixels + 1)/2)*image_pixel_size_mm;
    y = ((1:N_pixels) - (N_pixels + 1)/2)*image_pixel_size_mm;
    [X,Y] = meshgrid(x,y);
    
    % Build phantom
    
    phantom_radius_mm = PHANTOM_RADIUS_MM;
    phantom = X.^2 + Y.^2 < phantom_radius_mm^2;
    
    insert_radius_mm = 15;
    insert = X.^2 + Y.^2 < insert_radius_mm^2;
    
    insert_mask = X.^2 + Y.^2 < 0.85*insert_radius_mm^2;
    phantom_mask = X.^2 + Y.^2 < 0.9*phantom_radius_mm^2 & X.^2 + Y.^2 > 1.15*insert_radius_mm^2;
    
    
    [A_mm,FposArcDeg,Fangles] = fanbeam(...
        phantom,...
        fan_beam_D,...
        'FanSensorGeometry','arc',...
        'FanSensorSpacing',fan_sensor_spacing_deg,...
        'FanRotationIncrement',fan_rotation_increment_deg);
    
    A_mm = A_mm * image_pixel_size_mm;
    A_cm = A_mm / 10;
    
    ROI_size = 7.5;
    ROI_spacing = 4;
    
    water_mask = X.^2 + Y.^2 < 0.95*phantom_radius_mm^2;    
    
    % Create ROIs
    
    ROI_masks = cell(11);
    
    mask = 0;
    
    K = 6;
    
    for i = -K:K
        for j = -K:K
            
            ROI_masks{i+K+1,j+K+1} = ...
                ( max(abs(X - ROI_spacing*i*ROI_size),abs(Y - ROI_spacing*j*ROI_size)) <= ROI_size & water_mask);
            
            
            if sum(ROI_masks{i+K+1,j+K+1}(:)) < 0.9 * (2*ROI_size/image_pixel_size_mm)^2
                ROI_masks{i+K+1,j+K+1} = ROI_masks{i+K+1,j+K+1} & 0;
            end
            
            mask = mask | ROI_masks{i+K+1,j+K+1} ;
            
        end
    end      
    
    % Compute DQE
   
    for number_of_bins = 8
        
        switch number_of_bins
            case 8
                N_thresholds = 36;
            case 4
                N_thresholds = 36;
            case 2
                N_thresholds = 31;
        end
        
        real_ROI_means = nan([size(ROI_masks),N_ensemble_realizations,2,N_thresholds]);
        
        dqe = zeros(N_thresholds,2);
        
        real_I_1_tot = zeros(N_pixels,N_pixels,N_ensemble_realizations);
        real_I_2_tot = zeros(N_pixels,N_pixels,N_ensemble_realizations);                
        
        for k = [5 35]
            
            m_sino = matfile([design ...
                sprintf('sinogram_data_%d_cm',2*phantom_radius_mm/10) ...
                filesep sprintf('si_sinogram_%d_bin_T_%d.mat',[number_of_bins, k])]);
            
            A_est = m_sino.A;
            
            A_1_bias = A_cm - mean(A_est(:,:,1,:),[2 4]);
            A_2_bias = - mean(A_est(:,:,2,:),[2 4]);
            
            for n = 1:N_ensemble_realizations
                
                n
                
                tic
                
                A_n = A_est(:,:,:,n);
                
                I_1 = ifanbeam(...
                    A_n(:,:,1) + A_1_bias,...
                    fan_beam_D,...
                    'FanSensorGeometry','arc',...
                    'FanSensorSpacing',fan_sensor_spacing_deg,...
                    'FanRotationIncrement',fan_rotation_increment_deg,...
                    'OutputSize',N_pixels,...
                    'Filter',filter_name,...
                    'FrequencyScaling',frequency_scaling,...
                    'Interpolation',interpolation)/image_pixel_size_mm * 10;
                
                I_2 = ifanbeam(...
                    A_n(:,:,2) + A_2_bias,...
                    fan_beam_D,...
                    'FanSensorGeometry','arc',...
                    'FanSensorSpacing',fan_sensor_spacing_deg,...
                    'FanRotationIncrement',fan_rotation_increment_deg,...
                    'OutputSize',N_pixels,...
                    'Filter',filter_name,...
                    'FrequencyScaling',frequency_scaling,...
                    'Interpolation',interpolation)/image_pixel_size_mm * 10;
                
                real_ROI_means(:,:,n,1,k) = cellfun(@(x) mean(I_1(x)),ROI_masks);
                real_ROI_means(:,:,n,2,k) = cellfun(@(x) mean(I_2(x)),ROI_masks);
                
                real_I_1_tot(:,:,n) = I_1;
                real_I_2_tot(:,:,n) = I_2;
                
                toc                               
                
            end
            
            real_ensemble_ROI_covariances = nan([size(ROI_masks),2,2]);
            
            for i = 1:2*K+1
                for j = 1:2*K+1
                    
                    ROI_means = squeeze(real_ROI_means(i,j,:,:,k));
                    
                    ensemble_ROI_mean = mean(ROI_means);
                    
                    if any(isnan(ensemble_ROI_mean))
                        continue
                    end
                    real_ensemble_ROI_covariances(i,j,:,:) = cov(ROI_means);
                end
                
            end
            
            real_ensemble_covariance = squeeze(nanmean(real_ensemble_ROI_covariances,[1,2]));
            
            if k == 5 && PHANTOM_RADIUS_MM == 150
                real_ensemble_covariance_30_cm_5_keV = real_ensemble_covariance;
            end
            
            if k == 35 && PHANTOM_RADIUS_MM == 150
                real_ensemble_covariance_30_cm_35_keV = real_ensemble_covariance;
            end
            
            if k == 5 && PHANTOM_RADIUS_MM == 200
                real_ensemble_covariance_40_cm_5_keV = real_ensemble_covariance;
            end
            
            if k == 35 && PHANTOM_RADIUS_MM == 200
                real_ensemble_covariance_40_cm_35_keV = real_ensemble_covariance;
            end                                               
             
        end                
        
    end
    
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

X = F ./ F(:,1) * 1000

real_ensemble_covariance_30_cm_5_keV
real_ensemble_covariance_30_cm_35_keV
real_ensemble_covariance_40_cm_5_keV
real_ensemble_covariance_40_cm_35_keV

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
    
    
%     variance_ratio_40_cm(i) = (V'*real_ensemble_covariance_40_cm_5_keV*V) / ((V'*real_ensemble_covariance_40_cm_35_keV*V))
    
    
end

font_size = 14;
fig_dir = '..\..\Paper\Latex\Revision I\';
print_figs = false;

figure(1)
clf
semilogy(E,(variance_30_cm_5_keV),'-','LineWidth',1.5)
hold on
semilogy(E,(variance_30_cm_35_keV),'-.','LineWidth',1.5)
semilogy(E,(variance_40_cm_5_keV),':','LineWidth',1.5)
semilogy(E,(variance_40_cm_35_keV),'--','LineWidth',1.5)
grid on
xlim([35 120])
ylim([0.1 500])
xlabel('Synthetic monoenergy [keV]')
ylabel('ROI ensemble variance [HU^2]')
set(gca,'FontSize',font_size)

legend(...
    '30 cm water, 8-bin with T_1 = 5 keV',...
    '30 cm water, 8-bin with T_1 = 35 keV',...
    '40 cm water, 8-bin with T_1 = 5 keV',...
    '40 cm water, 8-bin with T_1 = 35 keV')



figure(2)
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

legend(...
    '30 cm water, 8-bin with T_1 = 5 keV / 35 keV',...
    '40 cm water, 8-bin with T_1 = 5 keV / 35 keV')



m_cov = matfile(sprintf('ensemble_covariances_ROIsize%1.2f_ROIspacing%d.mat',ROI_size,ROI_spacing)');

m_cov.real_ensemble_covariance_30_cm_5_keV = real_ensemble_covariance_30_cm_5_keV;
m_cov.real_ensemble_covariance_30_cm_35_keV = real_ensemble_covariance_30_cm_35_keV;
m_cov.real_ensemble_covariance_40_cm_5_keV = real_ensemble_covariance_40_cm_5_keV;
m_cov.real_ensemble_covariance_40_cm_35_keV = real_ensemble_covariance_40_cm_35_keV;



















