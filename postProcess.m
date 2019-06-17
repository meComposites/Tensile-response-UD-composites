function [outLevels,outSpec] = postProcess(lnSus,eRel,lnSucK,Vc,...
    Vf,Ef,Nthresh,DX,nXOut,pValues,...
    nf,Vs,VsI,iMax,is,Xinf)

%% Preliminary calculations
% Reducing the maximum stress considered from Xmax to XmaxOut, 
% or the number of stress increments from nX to nXOut
Xinf=Xinf(1:nXOut);
lnSus=lnSus(1:nXOut,:);
lnSucK=lnSucK(1:nXOut,:);
eRel=eRel(1:nXOut,:);
Vc=Vc(1:nXOut,:);

%% Calculating density of broken clusters (Fig.3, Block III.7; Eq.10)
% Probability of one broken cluster of level [0:iMax-1]
% in a level [1:iMax] bundle
SucK=real(exp(lnSucK));
% Density of broken clusters of levels [0:iMax-1] (Eq.10)
pCluster=SucK./Vc; 
% Zero-ing density of clusters for zero-stresses
pCluster(1,:)=0;

%% Post-processing for all bundle levels [0:iMax] (Fig.3, Block IV)

outLevels = postCalc(lnSus,pCluster,...
    Vf,Nthresh,DX,nf,Vs,nXOut,Xinf);

% Expected strains for stress-strain curve (Fig.3, Block IV.2)
% (for all Xinf)
outLevels.eExp=cumprod([Xinf./Ef,eRel(:,2:iMax+1)],2);
% Strain vector for plotting (same as expected strain)
outLevels.einfC=outLevels.eExp;
% Density of broken clusters (as calculated above)
outLevels.pCluster=pCluster;

%% Size effects on bundle strength (for all levels [0:iMax])

% Hierarchical level and number of fibres in all levels [0:iMax]
outLevels.i=0:iMax;
outLevels.nf=2.^outLevels.i;

% Calculating strength percentiles
% Pre-allocating (for all i, all probabilties in pValues)
outLevels.XProbC=zeros(iMax+1,length(pValues));

for j=1:iMax+1
    outLevels.XProbC(j,:)=...
        interp1q(outLevels.Fu(:,j),outLevels.XinfC,pValues')';
end

%% Post-processing for the specimen (no. of fibres: nfs; bundle level: is)
%  (Fig.3, Block IV)

% Interpolating lnSus for the exact level of the specimen is
lnSusI=interp1(0:iMax,lnSus',is)';

% Running post-processing function to generate outputs (Fig.3, Block IV)
[outSpec,mExpI] = postCalc(lnSusI,pCluster(:,1:ceil(is)),...
    Vf,Nthresh,DX,nf,VsI*ones(1,ceil(is)+1),nXOut,Xinf);

% Adding lnSusI to the output variable 
% (e.g. to calculate size effects along the length)
outSpec.lnSus = lnSusI;

% Interpolating eExp for the exact level of the specimen is
outSpec.eExp=interp1(0:iMax,outLevels.eExp',is)';
% Strain vector for plotting (same as expected strain)
outSpec.einfC=outSpec.eExp;

%% Broken clusters of all levels in specimen [0:floor(is)]

% Number of control volumes (for all Xinf, for all j-sublevels) 
nVcI=VsI./Vc(:,1:ceil(is));

% Number of clusters of all levels [0:ceil(is)], for all Xinf:
% Expected value (for all Xinf, and up to expected failure):
% (needed for the total number of fibre-breaks, see next section):
nClustersExpI=nVcI.*SucK(:,1:ceil(is));
% Zero-ing the number of clusters for zero stresses
nClustersExpI(1,:)=0;
% Expected number of clusters up to expected failure
outSpec.nClustersExp=nClustersExpI.*mExpI;
% Variance (needed for total number of fibre-breaks, see next section):
nClustersVarI=nVcI.*SucK(:,1:ceil(is)).*(1-SucK(:,1:ceil(is)));

% Pre-allocating matrix of probabilities 
% (for all Xinf, all cluster levels, all probabilties in pValues)
outSpec.nClustersProb=zeros(nXOut,ceil(is),length(pValues));

% Calculating distributions: at each Xinf and cluster level j=[0:ceil(is)], 
% the number of broken clusters follows a binomial distribution with 
% probability SucK(Xinf,j) and number of trials nVc(Xinf,j)
for p=1:length(pValues)
    outSpec.nClustersProb(:,:,p)=...
        binoinv(pValues(p),floor(nVcI),SucK(:,1:ceil(is)));
end

% Removing all NaN values (corresponding to SucK==0)
outSpec.nClustersProb(isnan(outSpec.nClustersProb))=0;

%% Statistical distribution of the total number of fibre breaks 

% Expected number of fibre-breaks (for all Xinf): 
nBreaksExpI=sum(nf(1:ceil(is)).*nClustersExpI,2);
% Variance: linear combination of variances of no. clusters:
nBreaksVarI=sum(nf(1:ceil(is)).^2.*nClustersVarI,2);

% Pre-allocating matrix of probability distribution of nBreaks (for all
% Xinf, and for the nProb number of probabilties to be calculated)
outSpec.nBreaksProb=zeros(nXOut,length(pValues));

% Distribution of total number of fibre-breaks. This assumes that the 
% linear combination of the number of broken clusters (each following a 
% binomial distribution) can be approximated by a normal distribution 
% (not 100% accurate, but it is the only analytical solution available)
for p=1:length(pValues)
    outSpec.nBreaksProb(:,p)=...
        norminv(pValues(p),nBreaksExpI,sqrt(nBreaksVarI));
end

%% Largest cluster for all probability levels

% For each Xinf, each probability: 
% 1. adds number of levels with nCluster>0 along level j (cumsum, dir 2);
% 2. find the maximum number of levels with nCluster>0 (max, direction 2);
% 3. re-arranges the Xinf*Levels*Probabilities matrix into a
% Xinf*Probabilitles matrix (permute);
% 4. sets j=0 for individual fibre breaks (column=1), j(j>0)=column-1,
% and j=-1 for no breaks (-1 at the end of calculation)
outSpec.clusterMaxLevelProb=...
    permute(max(cumsum(outSpec.nClustersProb>0,2),[],2),[1 3 2])-1;

% Calculating number of fibres in largest cluster (Eq.1)
outSpec.clusterMaxNfProb=2.^(outSpec.clusterMaxLevelProb);
% If there is no fibre break expected, sets number of fibre breaks to 0 
% (rather than 2^(-1)=1/2)
outSpec.clusterMaxNfProb(outSpec.clusterMaxNfProb<1)=0;

%% Specimen volume
outSpec.Vs=VsI;

end

function [out,mExp] = postCalc(lnSu,pClustersExp,...
    Vf,Nthresh,DX,nf,Vs,nX,Xinf)
% All calculations based on the cross-section of fibres only, 
% except for last section with suffix "C" 
% (which are based on the cross-section of the composite)

% Determining whether the calculations are 
% for all bundles (nLevels=jMax+1) or for a single specimen (nLevels=1)
jMax=size(pClustersExp,2);
nLevels=size(lnSu,2);

%% Bundle strength statistics
%  (Fig.3, Block IV.1; (see Pimenta & Pinho (2013) for further details)

%Statistical distribution
out.Xinf=Xinf;
out.Fu=1-exp(lnSu);

%Average bundle strength, based on cross-section of fibres only
intFuDX=trapz(out.Fu)*DX;
out.Xavg=Xinf(end)-intFuDX;

%CoV of bundle strength
intXFuDX=trapz(out.Fu.*Xinf)*DX;
out.CoVX=sqrt(Xinf(end).^2-out.Xavg.^2-2.*intXFuDX)./out.Xavg;

%% Expected stresses for stress-strain curve (up to expected failure)
%  (Fig.3, Block IV.2)

%Number of expected stress increments (up to Xavg)
nXExp=round(out.Xavg/DX)+1;

%Matrix of logical expected increments up to failure (all levels, all Xinf)
mExp=false(nX,nLevels);
for colOut=1:nLevels
    %Vector of expected stresses up to failure
    mExp(1:nXExp(colOut),colOut)=true();
end

out.XExp=repmat(Xinf,1,nLevels).*mExp;

%% Expected fibre breaks and clusters
%  (Fig.3, Block IV.3)

% Expected total density of fibre-breaks in level-j bundle (Eq.11)
% (for all levels [0:jMax]; column for j=0 is all zero)
out.pBreaksExp=[zeros(nX,1), cumsum(nf(1:jMax).*pClustersExp,2)].*mExp;
out.pBreaksExp=out.pBreaksExp(:,jMax+1-(nLevels-1):jMax+1);
% Expected total number of fibre-breaks in level-j bundle (Eq.12)
out.nBreaksExp=Vs(jMax+1-(nLevels-1):jMax+1).*out.pBreaksExp;

% Expected largest cluster in specimen, for all XExp
% (with density higher than Nthresh/Vs(j+1) )

% Pre-allocating matrix of largest cluster
out.clusterMaxLevelExp=-ones(nX,nLevels);

for colOut=1:nLevels
    % Calculating largest cluster level (Eq.14), 
    % starting with largest-level bundle (for Spec, only does ceil(is))
    % (If there is no fibre break expected, this sets clusterMaxLevel=-1)
    colj=jMax+1-(colOut-1);
    out.clusterMaxLevelExp(:,colOut)=...
        sum(cumsum(pClustersExp(:,1:colj-1)>(Nthresh/Vs(colj)))>0,2)-1;
end
%Turning clusterMaxLevel around for bundles, so that i=0 is in 1st column
out.clusterMaxLevelExp(:,:)=out.clusterMaxLevelExp(:,end:-1:1);

% Calculating number of fibres in largest cluster (Eq.1)
out.clusterMaxNfExp=2.^(out.clusterMaxLevelExp);
% If there is no fibre break expected, sets number of fibre breaks to 0 
% (rather than 2^(-1)=1/2)
out.clusterMaxNfExp(out.clusterMaxNfExp<1)=0;

% Capping results for Xinf>Xavg
out.clusterMaxLevelExp(~mExp)=-1;
out.clusterMaxNfExp=out.clusterMaxNfExp.*mExp;

%% Normalising outputs based on cross-section of composite 
%  (rather than the fibres only, suffix C)

%Stress vector for strength distribution
out.XinfC=Xinf*Vf;

%Average bundle strength
out.XavgC=out.Xavg*Vf;

%Stresses for stress-strain curve (up to expected failure)
out.XExpC=out.XExp*Vf;

%Expected total density of clusters and fibre-breaks
out.pClustersExp=pClustersExp;
out.pClustersExpC=pClustersExp*Vf;
out.pBreaksExpC=out.pBreaksExp*Vf;

end
