projdqes_0p8keV=load('./60mmActiveSi/projdqes_0p8keV');
projdqes_1p6keV=load('./60mmActiveSi/projdqes_1p6keV');
projdqes_3p2keV=load('./60mmActiveSi/projdqes_3p2keV');
%These files were saved as follows
% save projdqes_3p2keV x eNoiseSigmakeV dqe_density_unbinned dqe_density_8_bin dqe_density_4_bin dqe_density_2_bin dqe_spectral_unbinned dqe_spectral_8_bin dqe_spectral_4_bin dqe_spectral_2_bin
fig_dir = '../../Paper/Latex/Revision II/';
font_size = 14;
print_figs=true;
lineColor1=[0    0.4470    0.7410];
lineColor2=[.8500    0.3250    0.0980];
lineColor3=[0.9290    0.6940    0.1250];

%%
figure
grid on
plot(projdqes_0p8keV.x,projdqes_0p8keV.dqe_density_unbinned(1,:),'-o','Color',lineColor1,'LineWidth',1.5)
hold all
plot(projdqes_0p8keV.x,projdqes_0p8keV.dqe_density_4_bin(1,:),':x','Color',lineColor1,'LineWidth',1.5)

plot(projdqes_1p6keV.x,projdqes_1p6keV.dqe_density_unbinned(1,:),'-s','Color',lineColor2,'LineWidth',1.5)
plot(projdqes_1p6keV.x,projdqes_1p6keV.dqe_density_4_bin(1,:),':*','Color',lineColor2,'LineWidth',1.5)

plot(projdqes_3p2keV.x,projdqes_3p2keV.dqe_density_unbinned(1,:),'-v','Color',lineColor3,'LineWidth',1.5)
plot(projdqes_3p2keV.x,projdqes_3p2keV.dqe_density_4_bin(1,:),':^','Color',lineColor3,'LineWidth',1.5)

leg = legend(...
    '120 bins, \sigma_e=0.8 keV',...
    '4 bins,     \sigma_e=0.8 keV',...
    '120 bins, \sigma_e=1.6 keV',...
    '4 bins,     \sigma_e=1.6 keV',...
    '120 bins, \sigma_e=3.2 keV',...
    '4 bins,     \sigma_e=3.2 keV',...
    'Location','NorthEast');
leg.Position=leg.Position+[0.05 0 0 0];
text(19.1,0.45,{'30 cm water background','Density task'},'FontSize',12)

set(gca,'FontSize',font_size)

xlim([0 35])
ylim([0.3 0.7])

xlabel('Lowest threshold energy [keV]')
ylabel('DQE^{density}_{projection}')

if print_figs
    print('-depsc',[fig_dir 'Figure_noiselevels_a.eps'])
end
%%
pause(1) %Time to save plot
figure
grid on
plot(projdqes_0p8keV.x,projdqes_0p8keV.dqe_spectral_unbinned(1,:),'-o','Color',lineColor1,'LineWidth',1.5)
hold all
plot(projdqes_0p8keV.x,projdqes_0p8keV.dqe_spectral_4_bin(1,:),':x','Color',lineColor1,'LineWidth',1.5)

plot(projdqes_1p6keV.x,projdqes_1p6keV.dqe_spectral_unbinned(1,:),'-s','Color',lineColor2,'LineWidth',1.5)
plot(projdqes_1p6keV.x,projdqes_1p6keV.dqe_spectral_4_bin(1,:),':*','Color',lineColor2,'LineWidth',1.5)

plot(projdqes_3p2keV.x,projdqes_3p2keV.dqe_spectral_unbinned(1,:),'-v','Color',lineColor3,'LineWidth',1.5)
plot(projdqes_3p2keV.x,projdqes_3p2keV.dqe_spectral_4_bin(1,:),':^','Color',lineColor3,'LineWidth',1.5)

leg = legend(...
    '120 bins, \sigma_e=0.8 keV',...
    '4 bins,     \sigma_e=0.8 keV',...
    '120 bins, \sigma_e=1.6 keV',...
    '4 bins,     \sigma_e=1.6 keV',...
    '120 bins, \sigma_e=3.2 keV',...
    '4 bins,     \sigma_e=3.2 keV',...
    'Location','NorthEast');

text(2.5,0.35,{'30 cm water background','Spectral task'},'FontSize',12)

set(gca,'FontSize',font_size)

xlim([0 35])
ylim([0 0.45])

xlabel('Lowest threshold energy [keV]')
ylabel('DQE^{spectral}_{projection}')
set(leg, 'Position',[0.58   0.58 0.412499990420682   0.415476178555262])
if print_figs
    print('-depsc',[fig_dir 'Figure_noiselevels_b.eps'])
end
