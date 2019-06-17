function [lnSusJ,eRelJ,lnSucKJ] =...
    subBundleLevel(Ls,k,Ck,eta,LeO,lnSusO,lnSksO)
% Applying the hierarchical scaling law to a generic sub-bundle level j. 
% Suffix O: from previous level calculated (j-1)
% Suffix J: for current level being calculated (j)

%% Input survival probabilities (at the effective recovery length)
% (Fig.3, Block III.4; App.B)
lnSueO=LeO./Ls.*lnSusO;
lnSkeO=LeO./Ls.*lnSksO;

%% Hierarchical scaling law for bundle survival probabilities 
%  (Fig.3, Block III.5)

% Log-survival at the control length under uniform stresses
% (Eq.3;See Pimenta and Pinho (2013) for further details))
if k==1
    lnSucJ=2*lnSueO+log(2-exp(2*lnSueO));
else
    if Ck>3
        lnSucJ=4*lnSueO+log(1+2*exp(lnSkeO-3*lnSueO)-2*exp(lnSkeO-lnSueO));
        lnSucJ(lnSucJ==Inf)=lnSueO(lnSucJ==Inf)+lnSkeO(lnSucJ==Inf)+log(2+exp(3*lnSueO(lnSucJ==Inf)-lnSkeO(lnSucJ==Inf))-2*exp(2*lnSueO(lnSucJ==Inf)));
    else
        lnSucJ=lnSueO+lnSkeO+log(2+exp(3*lnSueO-lnSkeO)-2*exp(2*lnSueO)); 
        lnSucJ(lnSucJ==Inf)=4*lnSueO(lnSucJ==Inf)+log(1+2*exp(lnSkeO(lnSucJ==Inf)-3*lnSueO(lnSucJ==Inf))-2*exp(lnSkeO(lnSucJ==Inf)-lnSueO(lnSucJ==Inf)));
    end
end

% Log-survival at the control length under uniform stresses (App.B)
lnSusJ=Ls./(2*LeO).*lnSucJ;
lnSusJ(1)=0;

%% Hierarchical scaling law for strain non-linearity
% (Fig.3, Block III.6; Eq.9) 
% Relative non-linear strains
eRelJ=1+(eta-1)./((1+exp(3*lnSueO-lnSkeO)./(2*(1-exp(2*lnSueO)))));

%% Hierarchical scaling law for cluster density
%  (Fig.3, Block III.7; Eq. 10) 
% log-survival with one failed sub-bundle (Eq.5b)
lnSucKJ=log(2)+log(1-exp(2*lnSueO))+lnSueO+lnSkeO;

end

