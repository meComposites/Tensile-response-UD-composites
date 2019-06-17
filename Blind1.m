function inputs = Blind1()
% Defining inputs for strength model
% 2nd benchmarking exercise, Blind material set 1

%% Numerical variables
inputs.Xmax=50000; inputs.DX=10;
inputs.XmaxOut=10000;

%% Mechanical properties
% Input fibre strength distribution
inputs.L0f=10; inputs.X0f=4000; inputs.m=4; 
% Fibre modulus
inputs.Ef=200000;
% Shear-lag strength
inputs.Tsl=20;

%% Geometry of composite and specimen
inputs.Df=0.012; inputs.Vf=0.50;  
inputs.nfs=2000; inputs.Ls=4;
% Leave iMax empty to run the model only for the levels needed to analyse
% the specimen. Otherwise, set iMax equal to the maximum bundle level to be
% analysed
inputs.iMax=[];

%% Model-specific inputs
inputs.k=2;
inputs.Nthresh=0.5;
inputs.pValues=[0.01 0.10 0.50 0.90 0.99];

%% Saving .mat file
save('Blind1.mat');

end