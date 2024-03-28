clear
IDEAL_DETECTOR=false;
INCLUDE_INSERTS = true;
if INCLUDE_INSERTS
    insert_string='_inserts';
else
    insert_string='';
end
INCLUDE_LARGE_AREA_COV_ESTIMATION=false; %also performs calculation using a large-area covariance matrix. This gives a large area variance that agrees with the CRLB. Useful for debugging.

iodine_concentrations_relative_to_basis=[0.1 0.2 0.5 1]; %corresponds to 1 2 5 10 mg/ml since the basis conc is 10 mg/ml.

%%
for PHANTOM_RADIUS_MM = [150 200]
    DECOMPOSITION = 'WATER_IODINE';


    fprintf('Setup ... ')

    fig_dir = '../../Paper/Latex/Revision II/';

    print_figs = false;

    design = '60mmActiveSi/';

    % Load energy vectors and spectrum

    m_settings = matfile([design 'SimulationSettings.mat']);

    E = m_settings.eActVectorkeV;
    D = m_settings.eDepVectorkeV;

    N_inc = length(E);
    N_dep = length(D);

    post_patient_spectrum = m_settings.backgroundSpectrumAfterPatientmm_2';

    N_neighbors = m_settings.nPixelsx;

    fprintf('done!\n')

    % Compute gaussian kernel for adding electronic noise and Fano noise

    eNoiseSigmakeV = 1.6;

    % % Compute Fano noise 
    % 
    % fano = 0.115;
    % e_Si = 0.0036;
    % sigma_fano = sqrt(fano*e_Si*E);
    % 
    % halfw = ceil(eNoiseSigmakeV*3)*2 + 1;
    % x = -1-halfw:halfw;
    % dx = -halfw:halfw;
    % 
    % kernels = zeros(length(dx),length(D))
    % 
    % 
    % kernel = diff(1 + 0.5*erf((x+0.5)/(eNoiseSigmakeV*sqrt(2))));
    % 
    % 
    % 
    % 
    % 
    % 
    % 
    % %%
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
    % % Add electronic noise and Fano noise to PSF
    % 
    % PSF_E = convn(PSF,reshape(kernel,[1,1,1,length(kernel)]),'same');
    % 
    %Above commented and replaced with the following by MP
    large_area_psf_file = sprintf('large_area_psf_sigma_%1.1f.mat',eNoiseSigmakeV);
    m_large_area_psf = matfile([design large_area_psf_file]);
    PSF_E = m_large_area_psf.PSF_E;



    %%

    % Load 1d ACF with electronic noise

    fprintf('Loading large area ACF ... ')

    large_area_1d_acf_file = sprintf('large_area_1d_acf_sigma_%1.1f.mat',eNoiseSigmakeV);
    if exist([design large_area_1d_acf_file],'file')

        m_acf = matfile([design large_area_1d_acf_file]);
        large_area_1d_ACF_E = m_acf.ACF;

    else

        m_acf = matfile([design 'acf.mat']);
        large_area_1d_ACF_E = zeros(m_settings.nPixelsx,length(E),length(D),length(D));

        for i = 1:length(E)        

            i

            ACF_i = squeeze(m_acf.acfPerIncPh(:,:,:,:,i,:,:));        

            for j = 1:m_settings.nPixelsx

                large_area_ACF_E_ij = 0;

                for k = 1:m_settings.nPixelsy

                    % If central pixel:

                    if j == (m_settings.nPixelsx+1)/2 && k == (m_settings.nPixelsy+1)/2

                        % Add 1-d electronic noise convolution of diagonal

                        large_area_ACF_E_ij = large_area_ACF_E_ij + ...
                            diag(conv(diag(squeeze(ACF_i(j,k,:,:))),kernel,'same'));                                        

                    else

                        % Add 2-d electronic noise convolution

                        large_area_ACF_E_ij = large_area_ACF_E_ij + ...
                            conv2(kernel,kernel,squeeze(ACF_i(j,k,:,:)),'same');
                    end                                

                end

                large_area_1d_ACF_E(j,i,:,:) = large_area_ACF_E_ij;  
            end        

        end

        m_acf = matfile([design 'large_area_1d_acf.mat']);
        m_acf.ACF = large_area_1d_ACF_E;

    end

    fprintf('done!\n')

    % Image parameters

    N_pixels = 1024; % Change to 
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

    %% Build phantom

    phantom_radius_mm = PHANTOM_RADIUS_MM;
    phantom = X.^2 + Y.^2 < phantom_radius_mm^2;

    insert_radius_mm = 15;
    insert_placement_radius = 75;
    if INCLUDE_INSERTS %MP
        insert_1 = ...
            (X - insert_placement_radius/sqrt(2)).^2 + ...
            (Y - insert_placement_radius/sqrt(2)).^2 < insert_radius_mm^2;

        insert_2 = ...
            (X - insert_placement_radius/sqrt(2)).^2 + ...
            (Y + insert_placement_radius/sqrt(2)).^2 < insert_radius_mm^2;

        insert_3 = ...
            (X + insert_placement_radius/sqrt(2)).^2 + ...
            (Y + insert_placement_radius/sqrt(2)).^2 < insert_radius_mm^2;

        insert_4 = ...
            (X + insert_placement_radius/sqrt(2)).^2 + ...
            (Y - insert_placement_radius/sqrt(2)).^2 < insert_radius_mm^2;
    end
    % Compute sinograms

    [A_mm,FposArcDeg,Fangles] = fanbeam(...
        phantom,...
        fan_beam_D,...
        'FanSensorGeometry','arc',...
        'FanSensorSpacing',fan_sensor_spacing_deg,...
        'FanRotationIncrement',fan_rotation_increment_deg);

    if INCLUDE_INSERTS %MP
        [A_insert_1_mm,FposArcDeg,Fangles] = fanbeam(...
            insert_1,...
            fan_beam_D,...
            'FanSensorGeometry','arc',...
            'FanSensorSpacing',fan_sensor_spacing_deg,...
            'FanRotationIncrement',fan_rotation_increment_deg);

        [A_insert_2_mm,FposArcDeg,Fangles] = fanbeam(...
            insert_2,...
            fan_beam_D,...
            'FanSensorGeometry','arc',...
            'FanSensorSpacing',fan_sensor_spacing_deg,...
            'FanRotationIncrement',fan_rotation_increment_deg);

        [A_insert_3_mm,FposArcDeg,Fangles] = fanbeam(...
            insert_3,...
            fan_beam_D,...
            'FanSensorGeometry','arc',...
            'FanSensorSpacing',fan_sensor_spacing_deg,...
            'FanRotationIncrement',fan_rotation_increment_deg);

        [A_insert_4_mm,FposArcDeg,Fangles] = fanbeam(...
            insert_4,...
            fan_beam_D,...
            'FanSensorGeometry','arc',...
            'FanSensorSpacing',fan_sensor_spacing_deg,...
            'FanRotationIncrement',fan_rotation_increment_deg);
    end
    % Convert to units of mm

    A_mm = A_mm * image_pixel_size_mm;
    if INCLUDE_INSERTS
        A_insert_1_mm = iodine_concentrations_relative_to_basis(1) * A_insert_1_mm * image_pixel_size_mm;
        A_insert_2_mm = iodine_concentrations_relative_to_basis(2) * A_insert_2_mm * image_pixel_size_mm;
        A_insert_3_mm = iodine_concentrations_relative_to_basis(3) * A_insert_3_mm * image_pixel_size_mm;
        A_insert_4_mm = iodine_concentrations_relative_to_basis(4) * A_insert_4_mm * image_pixel_size_mm;
        A_all_bases_mm=cat(3,A_mm,A_insert_1_mm+A_insert_2_mm+A_insert_3_mm+A_insert_4_mm);
    else
        A_all_bases_mm=cat(3,A_mm,zeros(size(A_mm)));
    end

    %%

%     figure(1)
%     clf
%     imagesc(A_insert_1_mm)


    %%





    N_dels = size(A_mm,1);

    % Detector parameters

    pixel_area_mm2 = m_settings.pixelPitchxmm * m_settings.pixelPitchymm;

    % Source parameters

    mAs = 200;

    bowtie_mm = 2*phantom_radius_mm - mean(A_mm,2); %MP note: let the bowtie be adapted to the water phantom excluding inserts
    flux_per_mAs_mm2 = 1.8e6; 
    mAs_per_view = mAs / N_views;
    
    pre_bowtie_spectrum = m_settings.backgroundSpectrumBeforePatientmm_2;
    normalized_pre_patient_spectrum = pre_bowtie_spectrum / sum(pre_bowtie_spectrum );

    pre_bowtie_spectrum_per_view = normalized_pre_patient_spectrum * ...
        flux_per_mAs_mm2 * pixel_area_mm2 * mAs_per_view;

    enoise_data=load(sprintf('./electronic_noise_simulations/electronic_noise_data_%1.1fkeV.mat',eNoiseSigmakeV));
    
    enoise_normalization_factor=sum(pre_bowtie_spectrum_per_view)/enoise_data.pre_patient_photons_per_mm2;
    electronic_noise_variance_per_view=enoise_data.electronic_noise_variance*enoise_normalization_factor; %Note: the normalization factor carries a unit mm2 so this is now equivalent quanta per pixel and view
	
    % Filter spectra

    mu_water = m_settings.muWatermm_1';
    mu_bone = m_settings.muBonemm_1';
    mu_iodine =	m_settings.muDiluteImm_1';
    pre_patient_spectrum_per_view = pre_bowtie_spectrum_per_view.*...
        exp(-bowtie_mm.*mu_water');

    post_patient_spectrum_per_view = permute(pre_patient_spectrum_per_view,[1 3 2]).*...
        exp(-A_all_bases_mm(:,:,1).*permute(mu_water,[3 2 1])-A_all_bases_mm(:,:,2).*permute(mu_iodine,[3 2 1]));

    fprintf('done!\n')

    T_8_bin = [5     7    14    33    38    47    56    73   200]';
    T_4_bin = [5    10    34    56   200]';
    T_2_bin = [5    32   200]';

    N_bin_combinations = 3;
    bin_combinations = [8 4 2];

    T_all_bin = unique([T_8_bin;T_4_bin;T_2_bin]);

    N_bins = length(T_all_bin) - 1;

    B_all_to_8 = -diff(T_all_bin(1:end-1)' >= T_8_bin);
    B_all_to_4 = -diff(T_all_bin(1:end-1)' >= T_4_bin);
    B_all_to_2 = -diff(T_all_bin(1:end-1)' >= T_2_bin);




    %% Compute row PSF and row ACF

    slice_thickness_mm = 10;

    m_PSF = matfile([design 'intermediatePsf.mat']);   

    u = 2*diff(max(min(m_PSF.pixelBordersymm,slice_thickness_mm),-slice_thickness_mm));

    row_PSF = squeeze(sum(PSF_E,2)) * sum(u);
    row_ACF = large_area_1d_ACF_E * sum(u);

    % Loop over thresholds

    N_ensemble_realizations = 50;

    lowest_thresholds = 1:36;
    N_thresholds = length(lowest_thresholds);

    % Material decomposition real detector

    basis_pairs = strsplit(DECOMPOSITION,'_');

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

    if ~IDEAL_DETECTOR
        for thresholds_idx = 1:N_thresholds

            % tic                

            % Set thresholds

            T_0 = lowest_thresholds(thresholds_idx);
            T = max([T_0;T_all_bin(2:end)], T_0);
            enoise_diff_variance=-diff(interp1(enoise_data.thresholds_for_electronic_noise_variance_keV+0.5,electronic_noise_variance_per_view,T,'linear',0)); %increase threshold by 0.5 in order to match the viewpoint that a threeshold of for example 20 actually includes everything from 0.5 keV and up.
            fprintf('Simulating lowest threshold T_0 = %d\n',T_0)

            % Compute binned PSF and ACF

            N_bins = length(T) - 1;
            B = -diff(D >= T);

            row_bin_PSF = zeros(N_neighbors,length(E),N_bins);
            row_bin_ACF = zeros(N_neighbors,length(E),N_bins,N_bins);

            for j = 1:N_bins

                row_bin_PSF(:,:,j) = sum(row_PSF(:,:,B(j,:) == 1),3);

                for k = 1:N_bins

                    row_bin_ACF(:,:,j,k) = sum(row_ACF(:,:,B(j,:) == 1, B(k,:) == 1),[3 4]);

                end
            end                

            % Build detector PSF matrix

            I = zeros(N_dels*N_bins*N_neighbors*N_inc,1);
            J = zeros(N_dels*N_bins*N_neighbors*N_inc,1);
            S = zeros(N_dels*N_bins*N_neighbors*N_inc,1);

            idx = 1;

            for i = 1:N_dels

                for j = -(N_neighbors - 1)/2:(N_neighbors - 1)/2

                    j_pix = i + j;
                    j_psf = j + (N_neighbors - 1)/2+1;

                    if (j_pix - 1) < 0 || (j_pix - 1) >= N_dels
                        continue
                    end

                    for k = 1:N_bins

                        for l = 1:N_inc

                            I(idx) = (i - 1) * N_bins + k;
                            J(idx) = (j_pix - 1) * N_inc + l;
                            S(idx) = row_bin_PSF(j_psf,l,k);

                            idx = idx + 1;

                        end
                    end
                end
            end

            R_lambda = sparse(I(1:idx-1),J(1:idx-1),S(1:idx-1),N_dels*N_bins,N_dels*N_inc);

            clear I J S

            % Build detector covariance matrix for post patient spectrum

            weighted_row_bin_ACF = permute(sum(row_bin_ACF.*permute(mean(post_patient_spectrum_per_view,2),[4 3 2 5 1]),2),[5 1 3 4 2]);
            weighted_row_bin_ACF(:,round((size(weighted_row_bin_ACF,2)+1)/2),:,:)=weighted_row_bin_ACF(:,round((size(weighted_row_bin_ACF,2)+1)/2),:,:)+...
                permute(diag(enoise_diff_variance),[3 4 1 2])*sum(u);
            %MP this is an approximation of the noise covariance matrix due to the angle averaging, which needs to be mentioned in the paper

            idx = 1;

            I = zeros(N_dels*N_bins*N_neighbors*N_bins,1);
            J = zeros(N_dels*N_bins*N_neighbors*N_bins,1);
            S = zeros(N_dels*N_bins*N_neighbors*N_bins,1);

            for i = 1:N_dels

                for j = -(N_neighbors - 1)/2:(N_neighbors - 1)/2

                    j_pix = i + j;
                    j_psf = j + (N_neighbors - 1)/2+1;

                    for k = 1:N_bins

                        for l = 1:N_bins

                            if (j_pix - 1) < 0 || (j_pix - 1) >= N_dels
                                continue
                            end

                            I(idx) = (i - 1) * N_bins + k;
                            J(idx) = (j_pix - 1) * N_bins + l;
                            S(idx) = weighted_row_bin_ACF(i,j_psf,k,l);

                            idx = idx + 1;

                        end
                    end
                end
            end

            R_Sigma = sparse(I(1:idx-1),J(1:idx-1),S(1:idx-1),N_dels*N_bins,N_dels*N_bins);

            clear I J S

            % Compute projections and derivatives

            lambda = reshape(R_lambda * reshape(permute(post_patient_spectrum_per_view,[3 1 2]),[N_inc*N_dels,N_views]),[N_bins,N_dels,N_views]);

        %     asd


            dlambda_A1 = dimtools.insdim(reshape(R_lambda * reshape(-permute(post_patient_spectrum_per_view,[3 1 2]).*F(:,1),[N_inc*N_dels,N_views]),[N_bins,N_dels,N_views]),2);
            dlambda_A2 = dimtools.insdim(reshape(R_lambda * reshape(-permute(post_patient_spectrum_per_view,[3 1 2]).*F(:,2),[N_inc*N_dels,N_views]),[N_bins,N_dels,N_views]),2);

            dlambda_A = cat(2,dlambda_A1,dlambda_A2);

            % Compute noise realizations        
            tic
            nonzero_bins = repmat(~all(B == 0,2),N_dels,1);

            R = zeros(N_dels*N_bins);
            R(nonzero_bins,nonzero_bins) = chol(R_Sigma(nonzero_bins,nonzero_bins))';

            % Resparsify

            R = sparse(R);
            toc
            % Compute MD

        %     A_0 = [A_mm(:,1)' / 10; zeros(1,N_dels)];
            A_0 = permute(A_all_bases_mm,[3 2 1])/10; %MP
            A_est_real_detector = zeros(N_dels,N_views,2,N_ensemble_realizations,N_bin_combinations);
            if INCLUDE_LARGE_AREA_COV_ESTIMATION
                A_est_real_detector_large_area_cov = zeros(N_dels,N_views,2,N_ensemble_realizations,N_bin_combinations);
            end
            for n = 1:N_ensemble_realizations

                n
                tic;
                rng(n);

                white_noise_realization = randn(N_bins*N_dels,N_views);

                correlated_noise_realization = reshape(R*white_noise_realization,[N_bins,N_dels,N_views]);

                % Create noisy sinogram data

                Y = lambda + correlated_noise_realization;
                toc;

                for m = 1:N_bin_combinations
                    sprintf('n=%d, m=%d',n,m)
                    tic

                    number_of_bins = bin_combinations(m);

                    switch number_of_bins
                        case 8
                            B_bin = B_all_to_8;
                        case 4
                            B_bin = B_all_to_4;
                        case 2
                            B_bin = B_all_to_2;
                    end
                    
                    binned_enoise_diff_variance=B_bin*enoise_diff_variance;
                    if INCLUDE_LARGE_AREA_COV_ESTIMATION
                        large_area_covariance=B_bin*squeeze(mean(sum(weighted_row_bin_ACF,2),1))*B_bin'; %Added by MP to be able to assess the impact of the weighting matrix
                        large_area_inv_covarianceCA=cell(1,N_views);
                        [large_area_inv_covarianceCA{:}]=deal(sparse(inv(large_area_covariance)));
                        large_area_SigmaInv_sparse=blkdiag(large_area_inv_covarianceCA{:});
                    end
                    % Compute MD

                    A_est = zeros(N_dels,N_views,2);
                    A_est_large_area_cov = zeros(N_dels,N_views,2);
                    for i = 1:N_dels
                        Jtemp=reshape(B_bin*reshape(dlambda_A(:,:,i,:),N_bins, []),size(B_bin,1),size(dlambda_A,2),N_views);
                        col_in_J=repmat((1:size(Jtemp,1))',[1, size(Jtemp,2),size(Jtemp,3)])+permute(size(Jtemp,1)*(0:size(Jtemp,3)-1),[1 3 2]);
                        row_in_J=repmat(1:size(Jtemp,2),[size(Jtemp,1), 1,size(Jtemp,3)])+permute(size(Jtemp,2)*(0:size(Jtemp,3)-1),[1 3 2]);
                        J=sparse(col_in_J(:),row_in_J(:),Jtemp(:));
                        lambda_i = B_bin*permute(lambda(:,i,:),[1 3 2]);
                        Y_i = B_bin*permute(Y(:,i,:),[1 3 2]);

                        if sum(sum(lambda_i,2) > 0) < 2  %Don't run decomposition if there are less than two bins registering counts
                            continue
                        end
                        row_bin_variance_in_ml_model=lambda_i+binned_enoise_diff_variance*sum(u); %Does not have correlations since ML model does not include these.
                        SigmaInv = spdiags(1 ./ row_bin_variance_in_ml_model(:),0,numel(lambda_i),numel(lambda_i));
                        SigmaInv(isinf(SigmaInv)) = 0;

                        A_est(i,:,:) = (A_0(:,:,i) + reshape((J'*SigmaInv*J)\(J'*SigmaInv*(Y_i(:) - lambda_i(:))),[],N_views))';
                        if INCLUDE_LARGE_AREA_COV_ESTIMATION
                            A_est_large_area_cov(i,:,:) = (A_0(:,:,i) + reshape((J'*large_area_SigmaInv_sparse*J)\(J'*large_area_SigmaInv_sparse*(Y_i(:) - lambda_i(:))),[],N_views))';
                        end
                    end
                    toc
                    A_est_real_detector(:,:,:,n,m) = A_est;
                    if INCLUDE_LARGE_AREA_COV_ESTIMATION
                        A_est_real_detector_large_area_cov(:,:,:,n,m) = A_est_large_area_cov;
                    end
                end

            end        

            % Save sinograms to file

            for m = 1:N_bin_combinations

                number_of_bins = bin_combinations(m);

                m_sino = matfile([design ...
                    sprintf('sinogram_data_%d_cm%s',2*phantom_radius_mm/10,insert_string) ...
                    filesep sprintf('si_sinogram_%d_bin_T_%d.mat',[number_of_bins, T_0])]);

                m_sino.A = A_est_real_detector(:,:,:,:,m);
                if INCLUDE_LARGE_AREA_COV_ESTIMATION
                    m_sino.A_LA = A_est_real_detector_large_area_cov(:,:,:,:,m);
                end
            end

            % toc

        end
    else
        % Material decomposition ideal detector. changed by MP

        lambda = permute(post_patient_spectrum_per_view,[3 1 2]) * sum(u);

        dlambda_A1 = dimtools.insdim(-reshape(-permute(post_patient_spectrum_per_view,[3 1 2]).*F(:,1),[N_inc,N_dels,N_views]),2) * sum(u);
        dlambda_A2 = dimtools.insdim(-reshape(-permute(post_patient_spectrum_per_view,[3 1 2]).*F(:,2),[N_inc,N_dels,N_views]),2) * sum(u);

        dlambda_A = cat(2,dlambda_A1,dlambda_A2);

        % Preallocate MD

        A_est_ideal_detector = zeros(N_dels,N_views,2,N_ensemble_realizations);

        rng(0)

        for n = 1:N_ensemble_realizations

            % Create noisy sinogram data

            Y = zeros(N_inc,N_dels,N_views);

            tic

            parfor view = 1:N_views
                Y(:,:,view) = poissrnd(lambda(:,:,view));
            end    

            % Compute MD    

            A_0 = permute(A_all_bases_mm,[3 2 1])/10; %MP
            A_est = zeros(N_dels,N_views,2);        

            for i = 1:N_dels

                Jtemp=permute(dlambda_A(:,:,i,:),[1 2 4 3]);
                col_in_J=repmat((1:size(Jtemp,1))',[1, size(Jtemp,2),size(Jtemp,3)])+permute(size(Jtemp,1)*(0:size(Jtemp,3)-1),[1 3 2]);
                row_in_J=repmat(1:size(Jtemp,2),[size(Jtemp,1), 1,size(Jtemp,3)])+permute(size(Jtemp,2)*(0:size(Jtemp,3)-1),[1 3 2]);
                J=sparse(col_in_J(:),row_in_J(:),Jtemp(:));
                lambda_i = permute(lambda(:,i,:),[1 3 2]);
                Y_i = squeeze(Y(:,i,:));

                SigmaInv = spdiags(1 ./ lambda_i(:),0,numel(lambda_i),numel(lambda_i));
                SigmaInv(isinf(SigmaInv)) = 0;

                A_est(i,:,:) = (A_0(:,:,i) + reshape((J'*SigmaInv*J)\(J'*SigmaInv*(Y_i(:) - lambda_i(:))),[],N_views))';

            end        

            A_est_ideal_detector(:,:,:,n) = A_est;

            toc

        end

        m_sino = matfile([design ...
            sprintf('sinogram_data_%d_cm%s',2*phantom_radius_mm/10,insert_string) ...
            filesep 'ideal_detector_sinogram.mat']);

        m_sino.A = A_est_ideal_detector;
    end
end









