function [is,iMax,nf,Vs,VsI,eta,Ck,nX,nXOut,nK,Xinf,Le,Vc] = preCalc(inputs)
% Preliminary calculations for HSL (Blocks II.1, II.2, III.1 and III.2 in
% flowchart,Fig.3)

%% j-levels
%Bundle-level of specimen
is=log2(inputs.nfs);    
%Maximum (integer) bundle level for model
if isempty(inputs.iMax)
    iMax=ceil(is);
else
    iMax=max(inputs.iMax,ceil(is));
end
%Vector of sub-bundle levels
levelsJ=(0:iMax);       

%% Cross-section
%Fig.3, Block II.1
%Fibre cross-section (App.A)
Af=pi*inputs.Df^2/4;
%Fibre perimeter (App.A)
Cf=pi*inputs.Df;        
%Inter-fibre distance (square packing) (App.A)
s=(sqrt(pi)/(2*sqrt(inputs.Vf))-1)*inputs.Df;   

%Fig.3, Block III.1
%Number of fibres in sub-bundles for levels [0:iMax] (Eq.1)
nf=2.^levelsJ; 
%Sub-bundle cross-section area for levels [0:iMax] (App.A)
A=nf*Af;
%Sub-bundle perimeter for levels [0:iMax] (App.A)
C=3*Cf+4*((sqrt(nf)-1)*s+(sqrt(nf)-2)*Cf/2);
%Sub-bundle volume for levels [0:iMax] (Eq.12) 
%(based on the cross-section of the fibres only)
Vs=A*inputs.Ls;
%Specimen volume for the exact nfs (Eq.12)
%(based on the cross-section of the fibres only)
VsI=inputs.nfs*Af*inputs.Ls;

%% Global parameters (Fig.3, Block II.1)
%Non-linear strain factor (Eq.9)
eta=1+(inputs.k-1)/4;
%Stress concentration survival probability factor (App.B)
if inputs.k==1
    Ck=1;
else
    Ck=(inputs.k^(inputs.m+1)-1)/((inputs.k-1)*(inputs.m+1));
end

%% Stress vector (Fig.3, Block II.2)
%Number of increments in remote stress vector, for hierarchical
%calculations
nX=ceil(inputs.Xmax/inputs.DX)+1;
%Number of increments in remote stress vector, for output
nXOut=ceil(inputs.XmaxOut/inputs.DX)+1;
%Number of increments in stress-concentration vector
nK=floor((nX-1)/inputs.k+1);
%Generates the remote stress vector
%Based on cross-section of fibres only:
Xinf=(linspace(0,inputs.DX*(nX-1),nX))';

%% Local lengths and volumes (Fig.3, Block III.2)
%Effective recovery length for levels [0:iMax-1], all Xinf
%(NOT normalised by any other length) (Eq.4)
Le=2*A(1:iMax)./(C(1:iMax).*inputs.Tsl).*Xinf;
%Control volume for levels [1:iMax], all Xinf (Eq.4,Eq.10)
Vc=4*A(1:iMax).*Le;

end