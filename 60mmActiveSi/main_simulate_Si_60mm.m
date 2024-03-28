%Note on units: variables representing physical quantities have their units included as suffixes (_1 means ^-1 etc.) However, weighted quantities may
%have different dimensions in different parts of the matrix, if the different weight functions have different dimensions. For this reason the
%suffix wf indicates that the specified dimension needs to be multiplied by the dimension of the weight function in order to give the complete dimension of the variable in question.

addpath '../own_general_tools/spektr_wrappers'
addpath '../from_others/Spektr3/Spektr Code'
addpath ./source/

calculatePsf=true;
calculateAcf=true;
postprocessPsf=true;
transformPsfAndAcf=true;

enforcePsfAndAcfSymmetry=true; 
symmetrizePsfAndAcfByTranspose=false; %If x and y are interchangeable. No effect if enforcePsfSymmetry is false.
%Note that these symmetries do not hold for all detectors, for example if the scatter is directional or there are layers in one direction such as in a Si-strip detector.

%Careful - these parameters have to agree with the Monte Carlo simulation.
%One way to get rid of this risk of errors would be to save metadata in a certain .mat file.
nPsfSubsamplesPerPixelx=3; %Note that this can be even or odd. The first beam position is always in the pixel center. If odd, the last sub-beam just fits within the +x boundary of the pixel. If even, the last sub-beam position is centered on the last boundary of the pixel.
nPsfSubsamplesPerPixely=3; %Same comment applies as for nPsfSubsamplesPerPixelx
useSeparateSimulationsForShiftsInx=false;
useSeparateSimulationsForShiftsIny=true;
pixelPitchxmm=0.5;
pixelPitchymm=0.5;

intermediatePsfMatfileName='./Si_for_comparison_paper/60mmActiveSi/Si120kVp/intermediatePsf';
psfMatfileName='./Si_for_comparison_paper/60mmActiveSi/Si120kVp/psf';
acfMatfileName='./Si_for_comparison_paper/60mmActiveSi/Si120kVp/acf';
fourierMetricsMatfileName='./Si_for_comparison_paper/60mmActiveSi/Si120kVp/FourierMetrics';
simulationSettingsFilename='./Si_for_comparison_paper/60mmActiveSi/Si120kVp/SimulationSettings';
inputFilenamesAcfCA={'../../simulations/GATE_simulations_materialcomp_2018/Si_1keVsteps/Si_61mm_fullpixel/packaged_simulation_step1.mat'};

%To avoid rerunning the entire simulation in cases where only the postpatient spectrum has been changed,
%we can specify input files containing simulated intermediate psf and acf where the full spectral
%response has been saved. These will be weighted with the new spectrum.
reweightingInputIntermediatePsfMatfileName='';
reweightingInputAcfMatfileName='';

if ~useSeparateSimulationsForShiftsIny && ~useSeparateSimulationsForShiftsIny
    nInputFilenames=1;
    inputFilenamesPsfCA=cell(1,1);    
elseif useSeparateSimulationsForShiftsInx && ~useSeparateSimulationsForShiftsIny
    nInputFilenames=ceil(nPsfSubsamplesPerPixelx/2);
    inputFilenamesPsfCA=cell(nInputFilenames,1);
elseif ~useSeparateSimulationsForShiftsInx && useSeparateSimulationsForShiftsIny
    nInputFilenames=ceil(nPsfSubsamplesPerPixely/2);
    inputFilenamesPsfCA=cell(1,nInputFilenames);
else
    error('Using a 2D grid of simulation input files is currently not supported.');
end
[inputFilenamesPsfCA{1:nInputFilenames}]=deal('../../simulations/GATE_simulations_materialcomp_2018/Si_1keVsteps/Si_61mm_subpixel/packaged_simulation_step%d'); %For psf (note: changes below)
inputFilenamesPsfCA=cellfun(@(string,i)sprintf(string,i),inputFilenamesPsfCA,num2cell(1:nInputFilenames),'UniformOutput',false);

nPixelsx=41; %Should be odd so that there is a center pixel
nPixelsy=41; %Should be odd so that there is a center pixel

extPixelBordersxmm=-pixelPitchxmm*(nPixelsx/2+1):pixelPitchxmm:pixelPitchxmm*(nPixelsx/2+1);
extPixelBordersymm=-pixelPitchymm*(nPixelsy/2+1):pixelPitchymm:pixelPitchymm*(nPixelsy/2+1);
depthSegmentLengthsmm=[60];
deadLayerPathLengthInFrontOfActiveMaterial=0.5;
extPixelBorderszmm=[-inf deadLayerPathLengthInFrontOfActiveMaterial+[0 cumsum(depthSegmentLengthsmm)] inf];

waferNormalDirection='y'; %This is the Z coordinate in the Gate coordinate system

radiusRefmm=0.001997139816603;
% radiusRefmm=1e-6; %No charge sharing
ERefkeV=1;
chargeCloudRadiusmmFunction=@(EDepInInteractionkeV,distanceToBorderTopmm)radiusRefmm*(EDepInInteractionkeV/ERefkeV).^(0.528578007613418); %Fitted to synchrotron data Liu et al. TNS
chargeFractionCollectedFunction=@calculateChargeFractionGaussian;
eDepVectorkeV=1:1:150; %Must be equispaced
nPrimaryPhotPerSubset=10000;
saveFullSpectralResponse=true;

eActVectorkeV=20:1:150;

muWatermm_1=spektrMuRhoWrapper(eActVectorkeV,'Water');
muBonemm_1=spektrMuRhoWrapper(eActVectorkeV,'Bone');
muImm_1=spektrMuRhoWrapper(eActVectorkeV,'I');
iodineConcForBasismg_ml=10;
muDiluteImm_1=(iodineConcForBasismg_ml/4930)*muImm_1;
muGdmm_1=spektrMuRhoWrapper(eActVectorkeV,'Gd');
gadoliniumConcForBasismg_ml=1;
muDiluteGdmm_1=(gadoliniumConcForBasismg_ml/7900)*muGdmm_1;

muBasisFunctionsForPsfAcfmm_1 = [muWatermm_1; muBonemm_1; muDiluteImm_1; muDiluteGdmm_1];
nBasisFunctionsForPsfAcf = size(muBasisFunctionsForPsfAcfmm_1,1);

backgroundMaterialThicknessesmm=[300 0 0 0];
kVp=120;
mmAlFiltration=3; %Al & Ti filtration to mimic Siemens Force 
mmFiltrationInAdditionToAl=[0.9];
filtrationMaterialsInAdditionToAl={'Ti'};
muFilterMaterialsmm_1=zeros(length(filtrationMaterialsInAdditionToAl),length(eActVectorkeV));
for filtrationMaterialNo=1:length(filtrationMaterialsInAdditionToAl)
    muFilterMaterialsmm_1(filtrationMaterialNo,:)=spektrMuRhoWrapper(eActVectorkeV,filtrationMaterialsInAdditionToAl{filtrationMaterialNo});
end
percentRipple=0;
nPrepatientPhotonsPermm2=1e6; %Note: same fluence per area gives different numbers of prepatient photons per pixel depending on the pixel size

%Upsampling factors in the fourier domain. If 1, the number of sample points between -Nyq and +Nyq will be acfMatrixSize{x,y}.
%(Expected to be equal to the size of the pixel grid that psfmm_2 was simulated on, before including the beam shifts).
%If an integer >1, the spacing between frequency sample points is decreased by this factor.
%To keep the code simple I implemented it in a way that requires the upsampling factor to be an integer.
fourDomainUpsamplingFactorx = 1;
fourDomainUpsamplingFactory = 1;

nWorkersPsf=4; %Use fewer workers to avoid running out of memory
simulateDetectorMetrics