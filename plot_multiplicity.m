baseFolder='/home/persson6/Documents/MATLAB/spectral_detector_evaluation/Si_for_comparison_paper/60mmActiveSi/Si120kVp_with_multiplicity/';
load([baseFolder 'SimulationSettings.mat'])
load([baseFolder 'largeAreaMetrics.mat'])
load([baseFolder 'acf.mat'])
intermPsf=load([baseFolder 'intermediatePsf.mat']); %Avoid overwriting things loaded from the other variables

%%
lineWidth=1.5;
markersMonteCarlo={'o', 'x', 'v'};
markersGeomPoiss={'s', 'p', '^'};
fontSize=14;
%% 
%Note: electronic noise is not taken into account here, since we want the
%number of raw deposited events over a given threshold. Therefore the SPR
%and rho are not identical with the values in FIg 5b in the paper.
weightedMultiplicityPdfForCurrentSpectrumwf=permute(weightedMultiplicityPdfwf(1,1,end,:,:),[4 5 1 2 3]); 
weightedMultiplicityPdfForCurrentSpectrumPerIncPhoton=weightedMultiplicityPdfForCurrentSpectrumwf/sum(backgroundSpectrumAfterPatientmm_2);
thresholdsToPlotkeV= [5 15];
eDepIndicesForThresholds=find(ismember(eDepVectorkeV,thresholdsToPlotkeV));
multiplicities=1:size(weightedMultiplicityPdfForCurrentSpectrumPerIncPhoton,1);
totalEventsForThresholds=flip(cumsum(flip(permute(sum(largeAreaSpectralResponse.*backgroundSpectrumAfterPatientmm_2,2),[1 3 2]),2)),2);
uniquePhotonsForThresholds=squeeze(weightedUniquePhotonsAboveThresholdAfterChShwf(:,:,end,:))';
SPR=(totalEventsForThresholds-uniquePhotonsForThresholds)./uniquePhotonsForThresholds;
rho=SPR./(1+SPR);
legendStrings={};
for thresholdNo=1:length(thresholdsToPlotkeV)
    currThresholdIndex=eDepIndicesForThresholds(thresholdNo);
    currDetectionProbability=uniquePhotonsForThresholds(eDepIndicesForThresholds(thresholdNo))/sum(backgroundSpectrumAfterPatientmm_2);
    plot(multiplicities,weightedMultiplicityPdfForCurrentSpectrumPerIncPhoton(:,eDepIndicesForThresholds(thresholdNo))/currDetectionProbability,[markersMonteCarlo{thresholdNo} '-'],'lineWidth', lineWidth); %Probability given detected photon
%     plot(multiplicities,weightedMultiplicityPdfForCurrentSpectrumPerIncPhoton(:,eDepIndicesForThresholds(thresholdNo)),[markersMonteCarlo{thresholdNo} '-'],'lineWidth', lineWidth); %Probability given incident photon
    hold all
    currentGeomPoDistr=rho(eDepIndicesForThresholds(thresholdNo)).^(multiplicities-1)*(1-rho(eDepIndicesForThresholds(thresholdNo)));
    plot(multiplicities,currentGeomPoDistr,[markersGeomPoiss{thresholdNo} '--'],'lineWidth', lineWidth) %Probability given detected photon
%     plot(multiplicities,currDetectionProbability*currentGeomPoDistr,[markersGeomPoiss{thresholdNo} '--'],'lineWidth', lineWidth) %Probability given incident photon
    legendStrings=[legendStrings, sprintf('T= %d keV, Monte Carlo',thresholdsToPlotkeV(thresholdNo)),...
        sprintf('T= %d keV, Geometric',thresholdsToPlotkeV(thresholdNo))];
end
xlim([0 6])
hleg=legend(legendStrings); %Check if the keVs should be offset by 0.5
set(hleg,'Position',get(hleg,'Position')+[0.00 0.01 0 0])
set(hleg, 'FontSize',fontSize)
set(gca, 'FontSize',fontSize)
xlabel('Multiplicity')
ylabel('Probability')
grid on
%%
print('-depsc','../../Paper/Latex/Revision III/multiplicity.eps')
%%
figure
%FOr testing Sanity check of boundary conditions. 
%The sum of multiplicities should be the total number of unique counts. 
plot(sum(weightedMultiplicityPdfForCurrentSpectrumwf),'--','lineWidth', lineWidth)
hold all,
plot(squeeze(weightedUniquePhotonsAboveThresholdAfterChShwf(:,:,end,:)))
max(abs((squeeze(weightedUniquePhotonsAboveThresholdAfterChShwf(:,:,end,:))-sum(weightedMultiplicityPdfForCurrentSpectrumwf)')))
%This should give the total number of interactions. 
%%
weightedMultiplicityPdfForCurrentSpectrumFromIntermPsfSimwf=permute(sum(intermPsf.weightedMultiplicityPdfwf(:,:,end,:,:).*[1 2; 2 4]/9,[1, 2]),[4 5 1 2 3]);

plot(sum(weightedMultiplicityPdfForCurrentSpectrumwf.*multiplicities'),'--','lineWidth', 2)
%small discrepancy with totalEventsForThresholds due to stochastic differences PSF-ACF sim.
hold all
plot(sum(weightedMultiplicityPdfForCurrentSpectrumFromIntermPsfSimwf.*multiplicities'),':','lineWidth', 2) %Agrees well with totalEventsForThresholds but not totally perfect, why?
plot(eDepVectorkeV,totalEventsForThresholds,'k-')
%%
% plot(sum(squeeze(multiplicityPdf(:,:,100,:,:)))) 