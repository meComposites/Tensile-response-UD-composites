function HSLPlotsParam()

%% Generating data
% Defining names of materials to be analysed
paramMat(1)=cellstr('parametricNom');
paramMat(2)=cellstr('parametricHm');
paramMat(3)=cellstr('parametricLm');
paramMat(4)=cellstr('parametricHX0');
paramMat(5)=cellstr('parametricLX0');

% Defining labels of materials to be analysed
paramLabel(1)=cellstr('nominal');
paramLabel(2)=cellstr('high Weibull modulus');
paramLabel(3)=cellstr('low Weibull modulus');
paramLabel(4)=cellstr('high Weibull strength');
paramLabel(5)=cellstr('low Weibull strength');

% Number of materials to be analysed
nMat=length(paramMat);

for mat=nMat:-1:1
    [outParam(mat),~,inParam(mat)]=HSL(paramMat{mat});
end

inputs.matName='parametric';

%% Default formatting

% meComposites colours
meCblack=[0,0,0];
meCred=[0.5412,0,0];
meCblue=[0,0.4745,0.8314];
meCgreen=[0.4431,0.6824,0.1804];
meCyellow=[1,0.6431,0.1216];

meCcoloursParam=[meCblack;...
    meCred; 1-(1-meCred)*0.5;...
    meCblue; 1-(1-meCblue)*0.5;];

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
    'defaultAxesColorOrder',meCcoloursParam);

transparentFill=0.25;

%% Plot 1: Expected stress vs strain curve
% (eExp vs XExpC)

plotName='1_eExp_XExpC';

% Create figure & axes
figNow=figure;
axNow=gca;

% Formatting plot
% axes formatting
set(axNow,'XLim',[0 2.0],'YLim',[0 2.5]);
set(axNow.XAxis,'TickLabelFormat','%.1f');
set(axNow.YAxis,'TickLabelFormat','%.1f');
% label formatting
set(axNow.XLabel,'String','Remote strain (%)');
set(axNow.YLabel,'String','Composite stresses (GPa)');

% Plotting data
for mat=nMat:-1:1
    line(outParam(mat).eExp*100,outParam(mat).XExpC/1000,...
        'Color',meCcoloursParam(mat,:),'DisplayName',paramLabel{mat});
end

% Create legend
legend('Location','northwest');

% Saving figure to PDF
set(figNow,'renderer','openGL',...
    'PaperPositionMode','Auto',...
    'PaperUnits','Centimeters',...
    'PaperSize',[widthFig, heightFig]);
print([inputs.matName '_' plotName '.pdf'],'-dpdf','-r600');

%% Plot 2: Strength and strain distributions
% (einfC & XinfC vs CDF)

plotName='2_einfC_XinfC';

% Coordinated axes limits
einfCLim=[1.2 2.2];

% Create figure & axes for plot 1
figNow=figure;
axNow=gca;

% Formatting plot 1 (einf vs CDF)
% axes formatting
set(axNow,'XLim',einfCLim,'YLim',[0 100]);
set(axNow.XAxis,'TickLabelFormat','%.1f');
set(axNow.YAxis,'TickLabelFormat','%.0f');
% label formatting
set(axNow.XLabel,'String','Remote strain (%)');
set(axNow.YLabel,'String','CDF (%)');

% Plotting data
for mat=nMat:-1:1
    line(outParam(mat).einfC*100,outParam(mat).Fu*100,...
        'Color',meCcoloursParam(mat,:),'DisplayName',paramLabel{mat});
end

% Create legend
legend('Location','southeast');

% Saving figure to PDF
set(figNow,'renderer','openGL',...
    'PaperPositionMode','Auto',...
    'PaperUnits','Centimeters',...
    'PaperSize',[widthFig, heightFig]);
print([inputs.matName '_' plotName '.pdf'],'-dpdf','-r600');

%% Plot 3: Fibre-breaks number and density
% (einfC & XinfC vs nBreaks & pBreaksC)

plotName='3_einfC_nBreaks';

% Coordinated axes limits
einfCLim=[0 2.0];
nBreaksLim=[0 700];
pBreaksCLim=nBreaksLim/outParam(1).Vs/1000;

% Create figure & axes for plot
figNow=figure;
axNow=gca;

% Formatting plot 1 (einfC vs nBreaks)
% axes formatting
set(axNow,'XLim',einfCLim,'YLim',nBreaksLim);
set(axNow.XAxis,'TickLabelFormat','%.1f');
set(axNow.YAxis,'TickLabelFormat','%.0f');
% label formatting
set(axNow.XLabel,'String','Remote strain (%)');
set(axNow.YLabel,'String','Number of fibre-breaks');

% Plotting data
for mat=nMat:-1:1
    line(outParam(mat).eExp*100,outParam(mat).nBreaksExp,...
        'Color',meCcoloursParam(mat,:),'DisplayName',paramLabel{mat});
end

% Create legend
legend('Location','northwest');

% Create second set of axes (with composite stresses vs. density of breaks)
axNow=axes();
set(axNow,'XLim',einfCLim,'YLim',pBreaksCLim,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'box','off');
set(axNow.XAxis,'TickLabels',[]);
set(axNow.YAxis,'TickLabelFormat','%.1f');
% label formatting
set(axNow.YLabel,'String','Density of fibre-breaks (1000 breaks/mm^3)');

% Saving figure to PDF
set(figNow,'renderer','openGL',...
    'PaperPositionMode','Auto',...
    'PaperUnits','Centimeters',...
    'PaperSize',[widthFig, heightFig]);
print([inputs.matName '_' plotName '.pdf'],'-dpdf','-r600');

%% Plot 4: Cluster number and density
% (einfC & XinfC vs nClusters & pClustersC)

plotName='4_einfC_nClusters';

% Coordinated axes limits
einfCLim=[0 2.0];
nClustersLim=[0 700; 0 50; 0 6; 0 1.0];
pClustersCLim=nClustersLim/outParam(1).Vs;

% Looping through the cluster sizes to be plotted
jMaxPlot=3;

% Labels for plots
j=0;
nClustersLabel(j+1)=cellstr('Number of single-fibre breaks');
pClustersLabel(j+1)=cellstr('Density of single-fibre breaks (breaks/mm^3)');

for j=1:jMaxPlot
    nClustersLabel(j+1)=cellstr(['Number of clusters with ', num2str(2^j), ...
        ' to ', num2str(2^(j+1)-1), ' fibres']);
    pClustersLabel(j+1)=cellstr(['Density of clusters with ', num2str(2^j), ...
        ' to ', num2str(2^(j+1)-1), ' fibres (clusters/mm^3)']);
end

% Format for nClusters tick labels
nClustersTickFormat=['%.0f'; '%.0f'; '%.0f'; '%.1f'];
% FontSize of nClusters tick labels
nClustersFontSize=[12; 14; 16; 16];

for j=0:jMaxPlot
    
    % Create figure & axes for plot
    figNow=figure;
    axNow=gca;

    % Formatting plot 1 (einfC vs nClusters)
    % axes formatting
    set(axNow,'XLim',einfCLim,'YLim',nClustersLim(j+1,:));
    set(axNow.XAxis,'TickLabelFormat','%.1f');
    set(axNow.YAxis,'TickLabelFormat',nClustersTickFormat(j+1,:));
    % label formatting
    set(axNow.XLabel,'String','Remote strain (%)');
    set(axNow.YLabel,'String',nClustersLabel(j+1));

    % Plotting data
    for mat=nMat:-1:1
        line(outParam(mat).eExp*100,outParam(mat).nClustersExp(:,j+1),...
            'Color',meCcoloursParam(mat,:),'DisplayName',paramLabel{mat});
    end
    
    % Create legend
    legend('Location','northwest');
    
    % Create second set of axes (with composite stresses vs. density of clusters)
    axNow=axes();
    set(axNow,'XLim',einfCLim,'YLim',pClustersCLim(j+1,:),...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'box','off');
    set(axNow.XAxis,'TickLabels',[]);
    set(axNow.YAxis,'TickLabelFormat',nClustersTickFormat(j+1,:),...
        'FontSize',nClustersFontSize(j+1));
    % label formatting
    set(axNow.YLabel,'String',pClustersLabel(j+1));

    % Saving figure to PDF
    set(figNow,'renderer','openGL',...
        'PaperPositionMode','Auto',...
        'PaperUnits','Centimeters',...
        'PaperSize',[widthFig, heightFig]);
    print([inputs.matName '_'  plotName ...
        '_j' num2str(j) '.pdf'],'-dpdf','-r600');

end

%% Plot 5: Size of largest cluster
% (einfC & XinfC vs clusterMaxNf & clusterMaxLevel)

plotName='5_einfC_clusterMaxNf';

% Coordinated axes limits
einfCLim=[0 2.0];
clusterMaxLevelLim=[-1 5];

% Tick labels for clusterMaxNf
clusterMaxNfTickLabel(1)=cellstr('0');
clusterMaxNfTickLabel(2)=cellstr('1');

for j=1:clusterMaxLevelLim(2)
    clusterMaxNfTickLabel(j+2)=cellstr(...
        [num2str(2^j),'-', num2str(2^(j+1)-1)]);
end

% Create figure & axes for plot
figNow=figure;
axNow=gca;

% Formatting plot 1 (einfC vs nBreaks)
% axes formatting
set(axNow,'XLim',einfCLim,'YLim',clusterMaxLevelLim);
set(axNow.XAxis,'TickLabelFormat','%.1f');
set(axNow.YAxis,'TickLabelFormat','%.0f');
% label formatting
set(axNow.XLabel,'String','Remote strain (%)');
set(axNow.YLabel,'String','Number of fibres in largest cluster');
set(axNow.YAxis,'TickLabels',clusterMaxNfTickLabel,...
    'FontSize',12);
set(axNow.YLabel,'String','Number of fibres in largest cluster');

% Plotting data
for mat=nMat:-1:1
    line(outParam(mat).eExp*100,outParam(mat).clusterMaxLevelExp,...
        'Color',meCcoloursParam(mat,:),'DisplayName',paramLabel{mat});
end

% Create legend
legend('Location','northwest');

% Create second set of axes (with composite stresses vs. density of breaks)
axNow=axes();
set(axNow,'XLim',einfCLim,'YLim',clusterMaxLevelLim,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'box','off');
set(axNow.XAxis,'TickLabels',[]);
set(axNow.YAxis,'TickLabelFormat','%.0f');
axNow.YAxis.TickLabels(1)=cellstr(' ');
% label formatting
set(axNow.YLabel,'String','Level of largest cluster');

% Saving figure to PDF
set(figNow,'renderer','openGL',...
    'PaperPositionMode','Auto',...
    'PaperUnits','Centimeters',...
    'PaperSize',[widthFig, heightFig]);
print([inputs.matName '_' plotName '.pdf'],'-dpdf','-r600');

end

