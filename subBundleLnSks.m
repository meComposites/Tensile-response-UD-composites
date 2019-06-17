function lnSks=subBundleLnSks(k,Ck,lnSus,DX,nX,nK,Xinf)
% Survival probability under stress concentrations 
% (Fig.3, Block III.3; App.B)
% (For further details, see Pimenta & Pinho (2013))

%% Calculates log-survival under linear stresses
lnSls=cumtrapz(lnSus)*DX./Xinf;
lnSls(1)=0;

%% Calculates log-survival under stress concentrations
if mod(k,1)==0
    %If k is integer, takes every other k-value:
    lnSks(1:nK,1)=(k*lnSls(1:k:end)-lnSls(1:nK))./(k-1);
else
    %If k is not integer (interpolation required):   
    %Interpolates lnSlr for k*X (non-integer k):
    lnSlsK(1:nK,1)=interp1(X,lnSls,XK,'linear');
    %uses lnSlr interpolated @ k*X (lnSlrK)
    lnSks(1:nK,1)=(k*lnSlsK-lnSls(1:nK))./(k-1);
end

%% Asymptotic survival probability under stress concentrations
%(for Xinf>Xmax/k) 
lnSks(nK+1:nX,1)=Ck*lnSus(nK+1:nX);
lnSks(1)=0;

end