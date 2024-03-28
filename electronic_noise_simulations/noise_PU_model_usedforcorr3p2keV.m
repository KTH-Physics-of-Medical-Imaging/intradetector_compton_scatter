savefigs=false;
lineWidth=1.5;
fontSize=14;
%%
rng(1013) %To be able to regenerate the output
deadtimens=100;
clockTimens=10; %COnsistent with 1 ns discretization step and downsampling x10 in Gr?nberg et al Medical Physics, 45 (8), August 2018
eStepkeV=0.5;
thresholdskeV=-4:eStepkeV:15; %Extend the energy range compared to the 1.6 keV simulation
sigmakeV=3.2;
frametimens=100000; %As in Gr?nberg et al Medical Physics, 45 (8), August 2018
nCycles=round(frametimens/clockTimens);
cyclesPerDeadtime=round(deadtimens/clockTimens);
filterSignal =true;
%%
nChannels=10000;
tDiscr=1:20;
tauPns=40;
contInputSignal=normrnd(zeros(nChannels,nCycles),sigmakeV*ones(nChannels,nCycles));
if filterSignal %Non-white noise, for testing
    filterFunction=(tDiscr/(tauPns/clockTimens)).^2.*exp(2*(1-tDiscr./(tauPns/clockTimens))); %From Gr?nberg et al. 2018
    contInputSignal=conv2(contInputSignal,filterFunction);
    contInputSignal=contInputSignal(:,1:end-length(filterFunction)+1);
    contInputSignal=contInputSignal./mean(std(contInputSignal,[],2))*sigmakeV;
end
%%
discrInputSignal=int8(permute(contInputSignal,[1 3 2])>thresholdskeV);
outputSignal=int8(zeros(size(discrInputSignal)));
remainingDeadtimeCycles=int8(zeros(nChannels,length(thresholdskeV)));
tic
for cycleNo=1:nCycles
    remainingDeadtimeCycles=max(remainingDeadtimeCycles-1,0);
    triggered=(remainingDeadtimeCycles==0)&(discrInputSignal(:,:,cycleNo)>0);
    remainingDeadtimeCycles(triggered)=cyclesPerDeadtime;
    outputSignal(:,:,cycleNo)=triggered;
end
toc
%%
tic
registeredCounts=double(sum(outputSignal,3));
countsMean=mean(registeredCounts,1);
countsStd=std(registeredCounts,1);
toc
%%
p=0.5*(1-erf(thresholdskeV/(sqrt(2)*sigmakeV))); %Probability of signal above threshold
inputCountRateErfPU=p/clockTimens;
outputCountRateErfPU=inputCountRateErfPU./(1+inputCountRateErfPU*deadtimens);
outputCountsErfPU=outputCountRateErfPU*frametimens;
% outputCountRateErfPUMod=inputCountRateErfPU./(1+inputCountRateErfPU*deadtimens-inputCountRateErfPU*deadtimens/clockTimens);
outputCountRateErfPUMod=inputCountRateErfPU./(1+inputCountRateErfPU*(deadtimens-clockTimens));
outputCountsErfPUMod=outputCountRateErfPUMod*frametimens;
outputCountsVarErfPU=frametimens*inputCountRateErfPU./(1+inputCountRateErfPU*deadtimens).^3;
outputCountsVarErfPUMod=frametimens*inputCountRateErfPU./(1+inputCountRateErfPU*(deadtimens-clockTimens)).^3;
%%
%Renewal process model: mean renewal time = dead time + (cycle time)*(expected value of geometric distribution (defined on Z+) with parameter p)
%=> mean renewal time = (deadtime-cycle time) + (cycle time)/p Note the - cycle time since the geometric distribution is >=1 always and the minimum renewal time equals the deadtime.
%and total counts =frametimens/ mean renewal time
%variance of renewal time = (cycle time)^2*(1-p)/p^2
%variance of counts=frametimens*var(renewal time)/(mean renewal time)^3
outputCountRateRenProc=1./(deadtimens-clockTimens+clockTimens./p);
outputCountsRenProc=outputCountRateRenProc*frametimens; %=p*frametimens/(clockTimens+p*(deadtimes-clockTimens))=frametimens*inputCountRate/(1+inputCountRate*(deadtimes-clockTimens))
outputCountsVarRenProc=frametimens*clockTimens.^2*(1-p)./(p.^2.*(deadtimens-clockTimens+clockTimens./p).^3);

%%
plot(thresholdskeV, outputCountsRenProc/(frametimens*1e-9),'-','linewidth',lineWidth)
hold all,
plot(thresholdskeV,countsMean/(frametimens*1e-9),'o','linewidth',lineWidth,'MarkerSize',5)%'markerEdgeColor','none','markerFaceColor','k')

plot(thresholdskeV, outputCountsVarRenProc/(frametimens*1e-9),'--','linewidth',lineWidth)
plot(thresholdskeV,countsStd.^2/(frametimens*1e-9),'x','linewidth',lineWidth,'MarkerSize',5)%,'markerEdgeColor','none','markerFaceColor','y')

plot(thresholdskeV, 1/(deadtimens/1e9)*ones(size(thresholdskeV)),':','linewidth',lineWidth)

hleg=legend('Analytic mean','Simulated mean','Analytic variance','Simulated variance','1/\tau');
set(hleg,'fontSize',fontSize,'Location','East')
set(gca,'fontSize',fontSize)

set(hleg,'fontSize',fontSize,'Location','East')
set(hleg,'position',get(hleg,'position')+[0 0.15 0 0])
xlabel('Threshold (keV)','fontSize',fontSize)
ylabel({'Registered count rate (cps)','Registered count variance/frame time (s^{-1})'},'fontSize',fontSize)
ylim([-0.05 1.1]*1/(deadtimens/1e9))
xlim([min(thresholdskeV),max(thresholdskeV)])
grid on
if savefigs
    print('-depsc','Figure_noisemodel_a.eps')
end
%%
ylim([-0.005 0.08]*1/(deadtimens/1e9))
xlim([3,8])

if savefigs
    print('-depsc','Figure_noisemodel_b.eps')
end
%% OLD
% discrInputSignalConv=convn(discrInputSignal,ones(1,cyclesPerDeadtime),'full');
% discrInputSignalConv=discrInputSignalConv(:,1:end-cyclesPerDeadtime+1,:);
% resetSignal=convn(discrInputSignal,[zeros(1,cyclesPerDeadtime-1) 1],'full');
% resetSignal=resetSignal(:,1:end-cyclesPerDeadtime+1,:);
% currentPiledUpEvents=cumsum(discrInputSignal-resetSignal,2);
