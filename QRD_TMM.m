%Preamble to set figure text interpreter to Latex----------------------------------------------------%
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',18)
set(0,'DefaultFigureWindowStyle','normal')
path = convertCharsToStrings(fileparts(matlab.desktop.editor.getActiveFilename));
cd(path)
addpath(genpath('C:\Users\User\Dropbox\PhD\Projects\Fig 3.6'))

%% Constants
rho = 1.204; %density of air
c = 340; %speed of sound

%% Set parameters
N = 5; %number of wells and prime number
xlabels = 1:N;
fmax = Freq.f_max; % Maximum Freq of interest
fmin = Freq.f_min; % Minimum Freq of interest
df = Freq.df; %frequency step
                                                   
a = 3e-1; %size of half of panel in metres
W = 9e-2; %width of well in m.
step = lcm(N,100);
x = linspace(0,a,step);
theta = linspace(-pi/2,pi/2,181);

%% QRD depths
n = (0:(N-1));  %vector of desired sequence length
s_n = mod(n.^2,N);
dfreq = 500; %Hz
lambda_0 = c./dfreq; %design wavelength in meters
L = (s_n.*lambda_0)./(2*N); %well depths
L(L==0) = 1e-6;

%% Frequency vectors and QWR geometry
f_v = (fmin:df:fmax)';%frequency vector.
w = f_v.*2.*pi;
k = w./c; %wavenumber vector

%% Impedance of QWR (no losses)
z = rho.*c; %specific acoustic impedance of air at 20degC
Z_0 = z./W; %characteristic impedance of slit? Careful. This goes to infinity when slit is infinitely thin.
Zw = -1i*Z_0.*cot(L.*k);  %array is calculated by implicit expansion

%% Matrices
R = ((Zw./Z_0)-1)./((Zw./Z_0)+1);

%augmentation of R matrix so that Fraunhofer Integral can be calculated

%%
augind = length(x)./N;
Rbig = repelem(R, 1, augind);

%% Fraunhofer Integral

% Pstot = double.empty;
for ifr = 1:length(k)
   Ps(ifr,:) = fftfraunhofer(theta,Rbig(ifr,:),k(ifr),x);
end
PsA = abs(Ps);

Rflat = ones(size(Rbig));

for ifr = 1:length(k)
   Psflat(ifr,:) = fftfraunhofer(theta,Rflat(ifr,:),k(ifr),x);
end

%% Diffusion Coefficient

SI = abs(Ps).^2; %SPL is converted to sound intensity.
n_d = length(theta);

SIsum = sum(SI,2);
SIsq = sum(SI.^2,2);

delta_QRD = (SIsum.^2 - SIsq)./((n_d-1)*(SIsq)); %unnormalised diffusion coefficient

%% Flat plane
SIf = abs(Psflat).^2; %SPL is converted to sound intensity.

SIfsum = sum(SIf,2);
SIfsq = sum(SIf.^2,2);

deltaf = (SIfsum.^2 - SIfsq)./((n_d-1)*(SIfsq)); %unnormalised diffusion coefficient

%% Normalised Diffusion Coefficient
deltan = (delta_QRD - deltaf)./(1-deltaf);