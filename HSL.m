function [outSpec,outLevels,inputs] = HSL(matFile)
% Runs Hierarchical Scaling Law, based on the following papers:
%     Pimenta (2017), DOI:
%     Pimenta and Pinho (2013), DOI:
%     Pimenta et al (2018), DOI:
%     
% Inputs: '<matFile>.mat' input mat-file.
% This file must have a structure "inputs", with the following fields:
% Numerical (discrete remote stress vector)
%     Xmax:	Maximum stress considered in hierarchical calculations
%     DX:     Stress increment
% 	XmaxOut:Maximum remote stress for (non-hierarchical) output calculations
% Weibull strength distribution of fibres: 
%     L0f:    Gauge length
%     X0f:    Weibull scale parameter
%     m:      Weibull shape parameter / Weibull modulus
% Other mechanical properties: 
%     Ef:     Longitudinal Young's modulus of fibres 
%     Tsl:    Shear-lag strength
% Geometry of composite and specimen: 
%     Df:     Fibre diameter
%     Vf:     Fibre volume fraction in composite
%     nfs:    Number of fibres in specimen
%     Ls:     Gauge length of specimen
% Model-specific inputs: 
%     k:      Stress concentration factor. Set value=2
%     Nthresh:Threshold for considering a broken cluster. Set value=0.5
%     pValues:Values of percentiles for plotting
    
%% Reading inputs & preliminary calculations
load(matFile,'inputs');
inputs.matName=matFile;

%% Preliminary calculations
%is:    Specimen bundle level
%iMax:  Maximum (integer) bundle level for model
%nf:    Number of fibres in sub-bundles for levels [0:iMax]
%A:     Sub-bundle cross-section area for levels [0:iMax-1]
%eta:   Non-linear strain factor
%Ck:    Stress concentration survival probability factor
%nX:    Number of increments in remote stress vector
%nK:    Number of increments in stress-concentration vector
%Xinf:  Remote stress vector (from 0 to Xmax)
%Le:    Effective recovery length for levels [0:iMax-1], all Xinf
%Vc:    Control volume of level for levels [1:iMax], all Xinf
[is,iMax,nf,Vs,VsI,eta,Ck,nX,nXOut,nK,Xinf,Le,Vc] = preCalc(inputs);

%Preallocating matrices
%Log-survivals under uniform stresses at specimen length
%for levels [0:iMax], all Xinf (Eq.5)
lnSus=zeros(nX,iMax+1);
%Relative non-linear strains for levels [0:iMax], all Xinf (Eq.9)
%(eRel(:,1)=1, i.e. single-fibre level is linear-elastic upt to failure)
eRel=ones(nX,iMax+1);
%Log-survival with one failed sub-bundle 
%for levels [1:iMax], all Xinf (Eq.5b)
lnSucK=zeros(nX,iMax);

%% Single-fibre to 2-fibre sub-bundle level (j=0->1) (Fig.3, Block II.3)
j=1;
%Single fibre log-survivals at specimen length, 
%considering Weibull distribution:
% Under uniform stresses (Eq.2)
lnSus(:,j)=-inputs.Ls/inputs.L0f.*((Xinf./inputs.X0f).^inputs.m);
%lnSksJm1:  under stress concentrations (for current j-1=0 level) (App.B)
lnSksJm1=Ck*lnSus(:,j);

%Level j=1 outputs:
%eRel:      relative non-linear strains (for all levels [0:iMax], all Xinf)
%lnSucK:    Log-survival with one failed sub-bundle 
%           (for all levels [1:iMax], all Xinf)
[lnSus(:,j+1),eRel(:,j+1),lnSucK(:,j)] =...
    subBundleLevel(inputs.Ls,inputs.k,Ck,eta,...
    Le(:,j),lnSus(:,j),lnSksJm1);

%% Sub-bundle level j>1
for j=2:iMax   
    lnSksJm1=...
        subBundleLnSks(inputs.k,Ck,lnSus(:,j),inputs.DX,nX,nK,Xinf);
    
    [lnSus(:,j+1),eRel(:,j+1),lnSucK(:,j)] =...
        subBundleLevel(inputs.Ls,inputs.k,Ck,eta,...
        Le(:,j),lnSus(:,j),lnSksJm1);
end

%% Post-processing for all bundle levels [0:iMax]

[outLevels,outSpec] = postProcess(lnSus,eRel,lnSucK,Vc,...
    inputs.Vf,inputs.Ef,inputs.Nthresh,inputs.DX,nXOut,inputs.pValues,...
    nf,Vs,VsI,iMax,is,Xinf);

% Saving input and output variables
save([inputs.matName '_Results'],'outLevels','outSpec','inputs');

end