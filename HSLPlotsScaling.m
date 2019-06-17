function HSLPlotsScaling(outSpec,outLevels,inputs)

%% Default formatting

% meComposites colours
meCblack=[0,0,0];
meCred=[0.5412,0,0];
meCblue=[0,0.4745,0.8314];
meCgreen=[0.4431,0.6824,0.1804];
meCyellow=[1,0.6431,0.1216];

meCcolours=[meCblack;meCred;meCblue;meCgreen;meCyellow];

% Defining size of Figures (all units in cm)
widthPlot=11;
heightPlot=11;
margin=2.0;

widthFig=widthPlot+2*margin;
heightFig=heightPlot+2*margin;

figPos=[0 0 widthFig heightFig];
axPos=[margin margin widthPlot heightPlot];

% Setting default formatting
set(groot,...
    'defaultLineLineWidth', 2,...
    'defaultTextUnits', 'centimeters',...
    'defaultRootUnits', 'centimeters',...
    'defaultFigureUnits', 'centimeters',...
    'defaultFigurePosition', figPos,...
    'defaultAxesUnits', 'centimeters',...
    'defaultAxesPosition', axPos,...
    'defaultAxesFontName', 'Cambria',...
    'defaultAxesFontSize', 16,...
    'defaultAxesColorOrder',meCcolours);

transparentFill=0.25;

%% Plot 1: Size effects on cross-section (constant length)
% (XavgC and Xavg vs nf & i-level)

plotName='1_xSection';

% Coordinated axes limits
iLim=[0 length(outLevels.Xavg)-1];
nfLim=2.^iLim;

XinfCLim=[1000 3500]/1000;
XinfLim=XinfCLim./inputs.Vf;

% Create figure & axes for plot
figNow=figure;
axNow=gca;

% Formatting plot 1 (nf vs XavgC)
% axes formatting
xtick=10.^(0:floor(log10(nfLim(end))));
xticklabels=cellstr(num2str(round(log10(xtick(:))), '10^%d'));

set(axNow,'XLim',nfLim,'XScale','log','XTick',xtick,...
    'YLim',XinfCLim);
set(axNow.XAxis,'TickLabels',xticklabels,'TickLabelInterpreter','tex');
set(axNow.YAxis,'TickLabelFormat','%.1f');

% label formatting
set(axNow.XLabel,'String','No. fibres in cross-section');
set(axNow.YLabel,'String','Composite stresses (GPa)');

% Plotting expected value
line(outLevels.nf,outLevels.XavgC/1000,...
    'DisplayName','Expected value');

% Plotting probability lines
probColours=jet(length(inputs.pValues));
for p=1:length(inputs.pValues)
    line(outLevels.nf,outLevels.XProbC(:,p)./1000,...
        'Color',probColours(p,:),'LineWidth',1,...
        'DisplayName',[num2str(inputs.pValues(p)*100), 'th percentile']);
end

% Create legend
legend('Location','northeast');

% Create second set of axes (with fibre stresses vs. i-level)
axNow=axes();
set(axNow,'XLim',iLim,'YLim',XinfLim,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'box','off');
set(axNow.XAxis,'TickLabelFormat','%.0f');
set(axNow.YAxis,'TickLabelFormat','%.1f');
% label formatting
set(axNow.XLabel,'String','Hierarchical level of bundle');
set(axNow.YLabel,'String','Fibre stresses (GPa)');

% Saving figure to PDF
set(figNow,'renderer','openGL',...
    'PaperPositionMode','Auto',...
    'PaperUnits','Centimeters',...
    'PaperSize',[widthFig, heightFig]);
print([inputs.matName '_' plotName '.pdf'],'-dpdf','-r600');

%% Plot 2: Size effects on cross-section (constant length)
% (XavgC and Xavg vs ls)

plotName='2_length';

% Calculations of strength distributions
outSpec.LSize=[0.5 5.0 50 500 5000];
lnSusSize=outSpec.lnSus.*outSpec.LSize./inputs.Ls;
FuSize=1-exp(lnSusSize);

%Average bundle strength
intFuDXSize=trapz(FuSize)*inputs.DX;
outSpec.XavgCSize=(outSpec.Xinf(end)-intFuDXSize).*inputs.Vf;

% Pre-allocating (for all i, all probabilties in pValues)
outSpec.XProbCSize=zeros(length(outSpec.LSize),length(inputs.pValues));

for L=1:length(outSpec.LSize)
    outSpec.XProbCSize(L,:)=...
        interp1q(FuSize(:,L),outSpec.XinfC,inputs.pValues')';
end

% Coordinated axes limits
XinfCLim=[1000 2500]/1000;
XinfLim=XinfCLim./inputs.Vf;

% Create figure & axes for plot
figNow=figure;
axNow=gca;

% Formatting plot 1 (LSize vs XavgC)
% axes formatting
set(axNow,'XLim',[min(outSpec.LSize) max(outSpec.LSize)],...
    'XScale','log','XTick',outSpec.LSize,...
    'YLim',XinfCLim);
set(axNow.XAxis,'TickLabels',outSpec.LSize,'TickLabelFormat','%.1e');
set(axNow.YAxis,'TickLabelFormat','%.1f');

% label formatting
set(axNow.XLabel,'String','Specimen length (mm)');
set(axNow.YLabel,'String','Composite stresses (GPa)');

% Plotting expected value
line(outSpec.LSize,outSpec.XavgCSize/1000,...
    'DisplayName','Expected value');

% Plotting probability lines
probColours=jet(length(inputs.pValues));
for p=1:length(inputs.pValues)
    line(outSpec.LSize,outSpec.XProbCSize(:,p)./1000,...
        'Color',probColours(p,:),'LineWidth',1,...
        'DisplayName',[num2str(inputs.pValues(p)*100), 'th percentile']);
end

% Create legend
legend('Location','northeast');

% Create second set of axes (with fibre stresses vs. i-level)
axNow=axes();
set(axNow,'XLim',[min(outSpec.LSize) max(outSpec.LSize)],...
    'XScale','log','YLim',XinfLim,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'box','off');
set(axNow.XAxis,'TickLabels',[]);
set(axNow.YAxis,'TickLabelFormat','%.1f');
% label formatting
set(axNow.XLabel,'String','');
set(axNow.YLabel,'String','Fibre stresses (GPa)');

% Saving figure to PDF
set(figNow,'renderer','openGL',...
    'PaperPositionMode','Auto',...
    'PaperUnits','Centimeters',...
    'PaperSize',[widthFig, heightFig]);
print([inputs.matName '_' plotName '.pdf'],'-dpdf','-r600');

%% Saving results file
save([inputs.matName '_Results'],'outLevels','outSpec','inputs');

end

