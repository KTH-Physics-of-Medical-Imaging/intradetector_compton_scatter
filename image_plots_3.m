%MP oct 2023. Modification of image_plots_3 where I have removed the
%"debiasing" step that actually turned out to introduce bias.
clear

pauseTimeseconds=1;
reconstruct_all_realizations=true;
font_size = 14;
fig_dir = '..\..\Paper\Latex\Revision II\';
print_figs = true;
imageDescrTextx=-20;
imageDescrTexty=-23;
thresholdTextx=7.5;
thresholdTexty=-23;
thresholdText2x=8.7;
thresholdText2y=-23;
cax = [0.9 1.1];
icax = [-0.5 0.5];
textColor='black';

filter_name = 'hann';
frequency_scaling = 1;
interpolation = 'linear';


design = '60mmActiveSi/';

m_settings = matfile([design 'SimulationSettings.mat']);

PHANTOM_RADIUS_MM = 200;
DECOMPOSITION = 'WATER_IODINE';

T_0 = 5;

m_sino = matfile([design ...
    sprintf('sinogram_data_%d_cm_inserts',2*PHANTOM_RADIUS_MM/10) ...
    filesep sprintf('si_sinogram_%d_bin_T_%d.mat',[8, T_0])]);
                
A_est_5_keV = m_sino.A;

T_0 = 35;

m_sino = matfile([design ...
    sprintf('sinogram_data_%d_cm_inserts',2*PHANTOM_RADIUS_MM/10) ...
    filesep sprintf('si_sinogram_%d_bin_T_%d.mat',[8, T_0])]);
                
A_est_35_keV = m_sino.A;

V_5 = [ 0.3571   0.9341; 0.9341    -0.3571] ./ [0.3571; 0.9341];
V_35 = [0.3352   0.9421; 0.9421    -0.3352] ./ [0.3352; 0.9421];


tic

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

% Build phantom

phantom_radius_mm = PHANTOM_RADIUS_MM;
phantom = X.^2 + Y.^2 < phantom_radius_mm^2;

insert_radius_mm = 15;
% insert = X.^2 + Y.^2 < insert_radius_mm^2;

% insert_mask = X.^2 + Y.^2 < 0.85*insert_radius_mm^2; %obsolete
phantom_mask = X.^2 + Y.^2 < 0.9*phantom_radius_mm^2 & X.^2 + Y.^2 > 1.15*insert_radius_mm^2;


ROI_size = 7.5;
ROI_spacing = 4;

water_mask = X.^2 + Y.^2 < 0.95*phantom_radius_mm^2;

% Create ROIs

mask = 0;
% if ~INCLUDE_INSERTS %Grid of square ROIs
%     K = 6;
%     ROI_size = 7.5; %half side of square in  mm
%     ROI_spacing = 4;
%     circular_insert=false;
%     Xoffset=0;
%     Yoffset=0;
% else %Grid of circular ROIS, matching the locations of the inserts
    K = 1;
    ROI_size = 0.8*insert_radius_mm; %Radius of circle in mm
    ROI_spacing = (75*sqrt(2))/ROI_size ;
    circular_insert=true;
    XOffset=-ROI_spacing*ROI_size/2;
    YOffset=-ROI_spacing*ROI_size/2;
% end

ROI_masks = cell(K+1);

for i = 0:K%-K:K
    for j = 0:K%-K:K
        if circular_insert
            full_ROI_area=(ROI_size/image_pixel_size_mm)^2*pi;
            
            ROI_masks{i+1,j+1} = ...
                ( (X - XOffset -ROI_spacing*i*ROI_size).^2+(Y - YOffset- ROI_spacing*j*ROI_size).^2 <= ROI_size.^2 & water_mask);
        else
            full_ROI_area=(2*ROI_size/image_pixel_size_mm)^2;
            ROI_masks{i+1,j+1} = ...
                ( max(abs(X - XOffset - ROI_spacing*i*ROI_size),abs(Y - YOffset - ROI_spacing*j*ROI_size)) <= ROI_size & water_mask);
        end
        
        if sum(ROI_masks{i+1,j+1}(:)) < 0.9 * full_ROI_area
            ROI_masks{i+1,j+1} = ROI_masks{i+1,j+1} & 0;
        end
        
        mask = mask | ROI_masks{i+1,j+1} ;
        
    end
end

%Extra ROI mask for water in the phantom center.
ROI_mask_center=( (X).^2+(Y).^2 <= ROI_size.^2 & water_mask);
mask = mask | ROI_mask_center ;

ROI_masks=[ROI_masks(:)' ROI_mask_center];

[A_mm,FposArcDeg,Fangles] = fanbeam(...
    phantom,...
    fan_beam_D,...
    'FanSensorGeometry','arc',...
    'FanSensorSpacing',fan_sensor_spacing_deg,...
    'FanRotationIncrement',fan_rotation_increment_deg);


A_mm = A_mm * image_pixel_size_mm;
A_cm = A_mm / 10;


A_5_1_bias = 0; %Previously A_cm - mean(A_est_5_keV(:,:,1,:),[2 4]);
A_5_2_bias = 0; %previously - mean(A_est_5_keV(:,:,2,:),[2 4]);

A_35_1_bias = 0; %previously A_cm - mean(A_est_35_keV(:,:,1,:),[2 4]);
A_35_2_bias = 0; %Previously- mean(A_est_35_keV(:,:,2,:),[2 4]);

% Choose realization

n = 50;

I_1_5_keV = ifanbeam(...
    A_est_5_keV(:,:,1,n) + A_5_1_bias,...
    fan_beam_D,...
    'FanSensorGeometry','arc',...
    'FanSensorSpacing',fan_sensor_spacing_deg,...
    'FanRotationIncrement',fan_rotation_increment_deg,...
    'OutputSize',N_pixels,...
    'Filter',filter_name,...
    'FrequencyScaling',frequency_scaling,...
    'Interpolation',interpolation)/image_pixel_size_mm * 10;

I_2_5_keV = ifanbeam(...
    A_est_5_keV(:,:,2,n) + A_5_2_bias,...
    fan_beam_D,...
    'FanSensorGeometry','arc',...
    'FanSensorSpacing',fan_sensor_spacing_deg,...
    'FanRotationIncrement',fan_rotation_increment_deg,...
    'OutputSize',N_pixels,...
    'Filter',filter_name,...
    'FrequencyScaling',frequency_scaling,...
    'Interpolation',interpolation)/image_pixel_size_mm * 10;

I_1_35_keV = ifanbeam(...
    A_est_35_keV(:,:,1,n) + A_35_1_bias,...
    fan_beam_D,...
    'FanSensorGeometry','arc',...
    'FanSensorSpacing',fan_sensor_spacing_deg,...
    'FanRotationIncrement',fan_rotation_increment_deg,...
    'OutputSize',N_pixels,...
    'Filter',filter_name,...
    'FrequencyScaling',frequency_scaling,...
    'Interpolation',interpolation)/image_pixel_size_mm * 10;

I_2_35_keV = ifanbeam(...
    A_est_35_keV(:,:,2,n) + A_35_2_bias,...
    fan_beam_D,...
    'FanSensorGeometry','arc',...
    'FanSensorSpacing',fan_sensor_spacing_deg,...
    'FanRotationIncrement',fan_rotation_increment_deg,...
    'OutputSize',N_pixels,...
    'Filter',filter_name,...
    'FrequencyScaling',frequency_scaling,...
    'Interpolation',interpolation)/image_pixel_size_mm * 10;
        

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

  
X = F./F(:,1);

c_40 = X(E == 40,:);
c_70 = X(E == 70,:);
c_100 = X(E == 100,:);

x = ((1:N_pixels) - (N_pixels + 1)/2)*image_pixel_size_mm;
y = ((1:N_pixels) - (N_pixels + 1)/2)*image_pixel_size_mm;


figure(1)
clf
imagesc(x/10,x/10,I_1_5_keV * c_40(1) + I_2_5_keV * c_40(2))
colormap gray
axis tight equal
xlabel('x [cm]')
ylabel('y [cm]')
caxis(cax)
text(imageDescrTextx,imageDescrTexty,'40 keV mono','Color',textColor,'FontSize',14)
text(thresholdText2x,thresholdText2y,'T_1 = 5 keV','Color',textColor,'FontSize',14)


set(gca,'FontSize',font_size)
pause(pauseTimeseconds) % MP, to give figures time to refresh
if print_figs
    switch PHANTOM_RADIUS_MM
        case 150
            print('-depsc',[fig_dir 'Figure_9a.eps'])
        case 200
            print('-depsc',[fig_dir 'Figure_10a.eps'])
    end
end

figure(2)
clf
imagesc(x/10,x/10,I_1_35_keV * c_40(1) + I_2_35_keV * c_40(2))
colormap gray
axis tight equal
xlabel('x [cm]')
ylabel('y [cm]')
caxis(cax)
text(imageDescrTextx,imageDescrTexty,'40 keV mono','Color',textColor,'FontSize',14)
text(thresholdTextx,thresholdTexty,'T_1 = 35 keV','Color',textColor,'FontSize',14)

set(gca,'FontSize',font_size)

pause(pauseTimeseconds) % MP, to give figures time to refresh
if print_figs
    switch PHANTOM_RADIUS_MM
        case 150
            print('-depsc',[fig_dir 'Figure_9b.eps'])
        case 200
            print('-depsc',[fig_dir 'Figure_10b.eps'])
    end
end

figure(3)
clf
imagesc(x/10,x/10,I_1_5_keV * c_70(1) + I_2_5_keV * c_70(2))
colormap gray
axis tight equal
xlabel('x [cm]')
ylabel('y [cm]')
caxis(cax)

text(imageDescrTextx,imageDescrTexty,'70 keV mono','Color',textColor,'FontSize',14)
text(thresholdText2x,thresholdText2y,'T_1 = 5 keV','Color',textColor,'FontSize',14)

set(gca,'FontSize',font_size)

pause(pauseTimeseconds) % MP, to give figures time to refresh
if print_figs
    switch PHANTOM_RADIUS_MM
        case 150
            print('-depsc',[fig_dir 'Figure_9c.eps'])
        case 200
            print('-depsc',[fig_dir 'Figure_10c.eps'])
    end
end

figure(4)
clf
imagesc(x/10,x/10,I_1_35_keV * c_70(1) + I_2_35_keV * c_70(2))
colormap gray
axis tight equal
xlabel('x [cm]')
ylabel('y [cm]')
caxis(cax)

text(imageDescrTextx,imageDescrTexty,'70 keV mono','Color',textColor,'FontSize',14)
text(thresholdTextx,thresholdTexty,'T_1 = 35 keV','Color',textColor,'FontSize',14)

set(gca,'FontSize',font_size)

pause(pauseTimeseconds) % MP, to give figures time to refresh
if print_figs
    switch PHANTOM_RADIUS_MM
        case 150
            print('-depsc',[fig_dir 'Figure_9d.eps'])
        case 200
            print('-depsc',[fig_dir 'Figure_10d.eps'])
    end
end

figure(5)
clf
imagesc(x/10,x/10,I_1_5_keV * c_100(1) + I_2_5_keV * c_100(2))
colormap gray
axis tight equal
xlabel('x [cm]')
ylabel('y [cm]')
caxis(cax)

text(imageDescrTextx,imageDescrTexty,'100 keV mono','Color',textColor,'FontSize',14)
text(thresholdText2x,thresholdText2y,'T_1 = 5 keV','Color',textColor,'FontSize',14)

set(gca,'FontSize',font_size)

pause(pauseTimeseconds) % MP, to give figures time to refresh
if print_figs
    switch PHANTOM_RADIUS_MM
        case 150
            print('-depsc',[fig_dir 'Figure_9e.eps'])
        case 200
            print('-depsc',[fig_dir 'Figure_10e.eps'])
    end
end

figure(6)
clf
imagesc(x/10,x/10,I_1_35_keV * c_100(1) + I_2_35_keV * c_100(2))
colormap gray
axis tight equal
xlabel('x [cm]')
ylabel('y [cm]')
caxis(cax)

text(imageDescrTextx,imageDescrTexty,'100 keV mono','Color',textColor,'FontSize',14)
text(thresholdTextx,thresholdTexty,'T_1 = 35 keV','Color',textColor,'FontSize',14)

set(gca,'FontSize',font_size)

pause(pauseTimeseconds) % MP, to give figures time to refresh
if print_figs
    switch PHANTOM_RADIUS_MM
        case 150
            print('-depsc',[fig_dir 'Figure_9f.eps'])
        case 200
            print('-depsc',[fig_dir 'Figure_10f.eps'])
    end
end

figure(7)
clf
imagesc(x/10,x/10,I_1_5_keV)
colormap gray
axis tight equal
xlabel('x [cm]')
ylabel('y [cm]')
caxis(cax)

text(thresholdText2x,thresholdText2y,'T_1 = 5 keV','Color',textColor,'FontSize',14)
text(imageDescrTextx,imageDescrTexty,'Water image','Color',textColor,'FontSize',14)

pause(pauseTimeseconds) % MP, to give figures time to refresh
if print_figs
    switch PHANTOM_RADIUS_MM
        case 150
            print('-depsc',[fig_dir 'Figure_11a.eps'])
        case 200
            print('-depsc',[fig_dir 'Figure_12a.eps'])
    end
end

figure(8)
clf
imagesc(x/10,x/10,I_1_35_keV)
colormap gray
axis tight equal
xlabel('x [cm]')
ylabel('y [cm]')
caxis(cax)

text(thresholdText2x,thresholdText2y,'T_1 = 35 keV','Color',textColor,'FontSize',14)
text(imageDescrTextx,imageDescrTexty,'Water image','Color',textColor,'FontSize',14)

pause(pauseTimeseconds) % MP, to give figures time to refresh
if print_figs
    switch PHANTOM_RADIUS_MM
        case 150
            print('-depsc',[fig_dir 'Figure_11b.eps'])
        case 200
            print('-depsc',[fig_dir 'Figure_12b.eps'])
    end
end

figure(9)
clf
imagesc(x/10,x/10,I_2_5_keV)
colormap gray
axis tight equal
xlabel('x [cm]')
ylabel('y [cm]')
caxis(icax)

text(thresholdText2x,thresholdText2y,'T_1 = 5 keV','Color',textColor,'FontSize',14)
text(imageDescrTextx,imageDescrTexty,'Iodine image','Color',textColor,'FontSize',14)

pause(pauseTimeseconds) % MP, to give figures time to refresh
if print_figs
    switch PHANTOM_RADIUS_MM
        case 150
            print('-depsc',[fig_dir 'Figure_11c.eps'])
        case 200
            print('-depsc',[fig_dir 'Figure_12c.eps'])
    end
end

figure(10)
clf
imagesc(x/10,x/10,I_2_35_keV)
colormap gray
axis tight equal
xlabel('x [cm]')
ylabel('y [cm]')
caxis(icax)

text(thresholdText2x,thresholdText2y,'T_1 = 35 keV','Color',textColor,'FontSize',14)
text(imageDescrTextx,imageDescrTexty,'Iodine image','Color',textColor,'FontSize',14)

pause(pauseTimeseconds) % MP, to give figures time to refresh
if print_figs
    switch PHANTOM_RADIUS_MM
        case 150
            print('-depsc',[fig_dir 'Figure_11d.eps'])
        case 200
            print('-depsc',[fig_dir 'Figure_12d.eps'])
    end
end
%%

I_1_5_keV_all=zeros(size(I_2_35_keV,1),size(I_2_35_keV,2),size(A_est_5_keV,4));
I_1_35_keV_all=zeros(size(I_2_35_keV,1),size(I_2_35_keV,2),size(A_est_5_keV,4));
I_2_5_keV_all=zeros(size(I_2_35_keV,1),size(I_2_35_keV,2),size(A_est_5_keV,4));
I_2_35_keV_all=zeros(size(I_2_35_keV,1),size(I_2_35_keV,2),size(A_est_5_keV,4));
if reconstruct_all_realizations
    for n = 1:size(A_est_5_keV,4)
        I_1_5_keV(:,:,n) = ifanbeam(...
            A_est_5_keV(:,:,1,n) + A_5_1_bias,...
            fan_beam_D,...
            'FanSensorGeometry','arc',...
            'FanSensorSpacing',fan_sensor_spacing_deg,...
            'FanRotationIncrement',fan_rotation_increment_deg,...
            'OutputSize',N_pixels,...
            'Filter',filter_name,...
            'FrequencyScaling',frequency_scaling,...
            'Interpolation',interpolation)/image_pixel_size_mm * 10;

        I_2_5_keV_all(:,:,n) = ifanbeam(...
            A_est_5_keV(:,:,2,n) + A_5_2_bias,...
            fan_beam_D,...
            'FanSensorGeometry','arc',...
            'FanSensorSpacing',fan_sensor_spacing_deg,...
            'FanRotationIncrement',fan_rotation_increment_deg,...
            'OutputSize',N_pixels,...
            'Filter',filter_name,...
            'FrequencyScaling',frequency_scaling,...
            'Interpolation',interpolation)/image_pixel_size_mm * 10;

        I_1_35_keV_all(:,:,n) = ifanbeam(...
            A_est_35_keV(:,:,1,n) + A_35_1_bias,...
            fan_beam_D,...
            'FanSensorGeometry','arc',...
            'FanSensorSpacing',fan_sensor_spacing_deg,...
            'FanRotationIncrement',fan_rotation_increment_deg,...
            'OutputSize',N_pixels,...
            'Filter',filter_name,...
            'FrequencyScaling',frequency_scaling,...
            'Interpolation',interpolation)/image_pixel_size_mm * 10;

        I_2_35_keV_all(:,:,n) = ifanbeam(...
            A_est_35_keV(:,:,2,n) + A_35_2_bias,...
            fan_beam_D,...
            'FanSensorGeometry','arc',...
            'FanSensorSpacing',fan_sensor_spacing_deg,...
            'FanRotationIncrement',fan_rotation_increment_deg,...
            'OutputSize',N_pixels,...
            'Filter',filter_name,...
            'FrequencyScaling',frequency_scaling,...
            'Interpolation',interpolation)/image_pixel_size_mm * 10;
    end
    save(sprintf('./60mmActiveSi/sinogram_data_%d_cm_inserts/reconstruced_realizations_5_35_keV',round(PHANTOM_RADIUS_MM*2/10)),'I_1_5_keV_all', 'I_2_5_keV_all','I_1_35_keV_all', 'I_2_35_keV_all')
else
    load(sprintf('./60mmActiveSi/sinogram_data_%d_cm_inserts/reconstruced_realizations_5_35_keV',round(PHANTOM_RADIUS_MM*2/10)))
end
%% Table of iodine insert measurements.
basisIConcmgml_1=10;
nRealizations=size(I_2_5_keV_all,3);
ROIMaskIndexOrder=[5 4 2 1 3]; %Order of increasing I concentration
ROIMeansT5mgml_1=zeros(1,nRealizations);
ROIMeansT35mgml_1=zeros(1,nRealizations);

disp('mean I mg/ml')
for  ROIMaskIndex=ROIMaskIndexOrder 
    for realizationNo=1:nRealizations
        I_2_5_keV_curr=I_2_5_keV_all(:,:,realizationNo);
        I_2_35_keV_curr=I_2_35_keV_all(:,:,realizationNo);
        ROIMeansT5mgml_1(ROIMaskIndex,realizationNo)=squeeze(basisIConcmgml_1*mean(I_2_5_keV_curr(ROI_masks{ROIMaskIndex})));
        ROIMeansT35mgml_1(ROIMaskIndex,realizationNo)=squeeze(basisIConcmgml_1*mean(I_2_35_keV_curr(ROI_masks{ROIMaskIndex})));
        ROIStdsT5mgml_1(ROIMaskIndex,realizationNo)=squeeze(basisIConcmgml_1*std(I_2_5_keV_curr(ROI_masks{ROIMaskIndex})));
        ROIStdsT35mgml_1(ROIMaskIndex,realizationNo)=squeeze(basisIConcmgml_1*std(I_2_35_keV_curr(ROI_masks{ROIMaskIndex})));
    end
    fprintf('%1.3f $\\pm$ %1.3f & ',mean(ROIMeansT5mgml_1(ROIMaskIndex,:)),std(ROIMeansT5mgml_1(ROIMaskIndex,:)))    
    fprintf('%1.3f $\\pm$ %1.3f & ',mean(ROIMeansT35mgml_1(ROIMaskIndex,:)),std(ROIMeansT35mgml_1(ROIMaskIndex,:)))    
    fprintf('\n')
end

fprintf('\n')
disp('std dev I mg/ml')
for ROIMaskIndex=[5 4 2 1 3] %Order of increasing I concentration
    fprintf('%1.3f $\\pm$ %1.3f & ',mean(ROIStdsT5mgml_1(ROIMaskIndex,:)),std(ROIStdsT5mgml_1(ROIMaskIndex,:)))    
    fprintf('%1.3f $\\pm$ %1.3f & ',mean(ROIStdsT35mgml_1(ROIMaskIndex,:)),std(ROIStdsT35mgml_1(ROIMaskIndex,:)))    
    fprintf('\n')
end

%% Sanity check: standard deviation in the center of the phantom
centerROIVectors=basisIConcmgml_1*squeeze((reshape(I_2_5_keV_all(512-70:512+70,512-70:512+70,:),[],nRealizations)));
mean(std(centerROIVectors))
std(std(centerROIVectors))
