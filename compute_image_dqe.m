clear

INCLUDE_INSERTS = false;
if INCLUDE_INSERTS
    insert_string='_inserts';
else
    insert_string='';
end

% Reconstruction parameters

filter_name = 'hann';
frequency_scaling = 1;
interpolation = 'linear';

print_figs = false;

design = '60mmActiveSi/';

% Load energy vectors and spectrum

m_settings = matfile([design 'SimulationSettings.mat']);

for PHANTOM_RADIUS_MM = [150] %150 200 %Do one at a time
    
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
%     insert = X.^2 + Y.^2 < insert_radius_mm^2;
    
%     insert_mask = X.^2 + Y.^2 < 0.85*insert_radius_mm^2;
    phantom_mask = X.^2 + Y.^2 < 0.9*phantom_radius_mm^2 & X.^2 + Y.^2 > 1.15*insert_radius_mm^2;
    
    
    [A_mm,FposArcDeg,Fangles] = fanbeam(...
        phantom,...
        fan_beam_D,...
        'FanSensorGeometry','arc',...
        'FanSensorSpacing',fan_sensor_spacing_deg,...
        'FanRotationIncrement',fan_rotation_increment_deg);
    
    A_mm = A_mm * image_pixel_size_mm;
    
    water_mask = X.^2 + Y.^2 < 0.95*phantom_radius_mm^2;
    
    % Create ROIs
    
    
    mask = 0;
    if ~INCLUDE_INSERTS %Grid of square ROIs
        K = 6;
        ROI_size = 7.5; %half side of square in  mm
        ROI_spacing = 4;
        circular_insert=false;
    else %Grid of circular ROIS, four of which match the locations of the inserts
        K = 4;
        ROI_size = 0.85*insert_radius_mm; %Radius of circle in mm
        ROI_spacing = (75/2*sqrt(2))/ROI_size ;
        circular_insert=true;

    end    
    
    ROI_masks = cell(2*K+1);
    
    for i = -K:K
        for j = -K:K
            if circular_insert
                full_ROI_area=(ROI_size/image_pixel_size_mm)^2*pi;

                ROI_masks{i+K+1,j+K+1} = ...
                    ( (X - ROI_spacing*i*ROI_size).^2+(Y - ROI_spacing*j*ROI_size).^2 <= ROI_size.^2 & water_mask);                
            else
                full_ROI_area=(2*ROI_size/image_pixel_size_mm)^2;
                ROI_masks{i+K+1,j+K+1} = ...
                    ( max(abs(X - ROI_spacing*i*ROI_size),abs(Y - ROI_spacing*j*ROI_size)) <= ROI_size & water_mask);
            end
            
            if sum(ROI_masks{i+K+1,j+K+1}(:)) < 0.9 * full_ROI_area
                ROI_masks{i+K+1,j+K+1} = ROI_masks{i+K+1,j+K+1} & 0;
            end
            
            mask = mask | ROI_masks{i+K+1,j+K+1} ;
            
        end
    end
    % Compute ensemble covariance for ideal detector
    
    N_ensemble_realizations = 50;
    
    ideal_ROI_means = zeros([size(ROI_masks),N_ensemble_realizations,2]);
    
    m_sino = matfile([design ...
        sprintf('sinogram_data_%d_cm%s',2*phantom_radius_mm/10,insert_string) ...
        filesep 'ideal_detector_sinogram.mat']);
    
    A_est = m_sino.A;
    
    ideal_I_1_tot = zeros(N_pixels,N_pixels,N_ensemble_realizations);
    ideal_I_2_tot = zeros(N_pixels,N_pixels,N_ensemble_realizations);
    
    
    A_cm = A_mm * image_pixel_size_mm/10;
    
    A_1_bias = 0; %A_cm - mean(A_est(:,:,1,:),[2 4]); %Removed by MP
    A_2_bias = 0; %- mean(A_est(:,:,2,:),[2 4]); %Removed by MP
    
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
        
        ideal_I_1_tot(:,:,n) = I_1;
        ideal_I_2_tot(:,:,n) = I_2;
        
        ideal_ROI_means(:,:,n,1) = cellfun(@(x) mean(I_1(x)),ROI_masks);
        ideal_ROI_means(:,:,n,2) = cellfun(@(x) mean(I_2(x)),ROI_masks);
        
        toc
        %
        %     figure(1)
        %     clf
        %     imagesc(I_1)
        %     axis tight equal
        %     colormap gray
        % %     caxis([0.9 1.1])
        %     pause
        
    end
    
    % Compute ensemble covariance
    
    ideal_ensemble_ROI_covariances = nan([size(ROI_masks),2,2]);
    
    for i = 1:2*K+1
        for j = 1:2*K+1
            
            ROI_means = squeeze(ideal_ROI_means(i,j,:,:));
            
            ensemble_ROI_mean = mean(ROI_means);
            
            if any(isnan(ensemble_ROI_mean))
                continue
            end
            ideal_ensemble_ROI_covariances(i,j,:,:) = cov(ROI_means);
            
        end
    end
    
    ideal_ensemble_covariance = squeeze(nanmean(ideal_ensemble_ROI_covariances,[1,2]));
    
    % Compute DQE
   
    for number_of_bins = [8 4 2]
        
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
        
        for k = 1:N_thresholds
            
            m_sino = matfile([design ...
                sprintf('sinogram_data_%d_cm%s',2*phantom_radius_mm/10,insert_string) ...
                filesep sprintf('si_sinogram_%d_bin_T_%d.mat',[number_of_bins, k])]);
            
            A_est = m_sino.A;
            
            A_1_bias = 0; %A_cm - mean(A_est(:,:,1,:),[2 4]); %Removed by MP
            A_2_bias = 0; %- mean(A_est(:,:,2,:),[2 4]); %Removed by MP
            
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
                
                %                     figure(2)
                %                     clf
                %                     imagesc(I_1)
                %                     axis tight equal
                %                     colormap gray
                %             %         caxis([0.9 1.1])
                %                     pause
                
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
            
            real_ensemble_covariance = squeeze(nanmean(real_ensemble_ROI_covariances,[1,2]))
            
            
            [V,~] = eig(real_ensemble_covariance)
            
            dqe(k,:) = diag(V'*(real_ensemble_covariance \ V)) ./ diag(V'*(ideal_ensemble_covariance \ V));
            
            figure(1)
            set(gca,'ColorOrderIndex',5 -log2(number_of_bins));            
            plot(0:k-1,dqe(1:k,1))
            ylim([0.3 0.7])
            xlim([0 35])
            
            drawnow
            
            figure(2)
            set(gca,'ColorOrderIndex',5 -log2(number_of_bins));            
            plot(0:k-1,dqe(1:k,2))
            ylim([0 0.45])
            xlim([0 35])
            
            drawnow
             
        end
        
        m_dqe = matfile([design ...
                sprintf('sinogram_data_%d_cm%s',2*phantom_radius_mm/10,insert_string) ...
            filesep sprintf('si_cnr_dqe_%d_bin_%s_interpolation_ROIsize%1.2f_ROIspacing%d.mat',number_of_bins,interpolation,ROI_size,ROI_spacing)]);
        
        m_dqe.dqe = dqe;
        m_dqe.real_ROI_means = real_ROI_means;
        m_dqe.ideal_ROI_means = ideal_ROI_means;
        m_dqe.real_ensemble_ROI_covariances = real_ensemble_ROI_covariances;
        m_dqe.ideal_ensemble_ROI_covariances = ideal_ensemble_ROI_covariances;
        
    end
    
end





