function inputs = T700K315()
% Defining inputs for strength model
% 2nd benchmarking exercise, Blind material set 1

%% Numerical variables
inputs.Xmax=50000; inputs.DX=10;
inputs.XmaxOut=10000;

%% Mechanical properties
% Input fibre strength distribution
inputs.L0f=12; inputs.X0f=4219; inputs.m=4.14; 
% Fibre modulus
inputs.Ef=230000;
% Shear-lag strength
inputs.Tsl=60.4;

%% Geometry of composite and specimen
inputs.Df=0.0068; inputs.Vf=0.579;  
inputs.nfs=4430; inputs.Ls=1.2;
% Leave iMax empty to run the model only for the levels needed to analyse
% the specimen. Otherwise, set iMax equal to the maximum bundle level to be
% analysed
inputs.iMax=[];

%% Model-specific inputs
inputs.k=2;
inputs.Nthresh=0.5;
inputs.pValues=[0.01 0.10 0.50 0.90 0.99];

%% Saving .mat file
save('T700K315.mat');

end