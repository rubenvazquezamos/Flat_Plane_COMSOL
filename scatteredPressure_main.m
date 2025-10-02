% 1 Run COMSOL to MATLAB LiveLink
% 2 Then run this script.

% This code runs a numerical simulation of an acoustic 
% diffuser using COMSOL® and compares the obtained
% diffusion coefficient to that obtained using the
% transfer matrix method (TMM).

%% SESSION START UP COMMANDS
%-------------------------------------------------------------------------%
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',18)
set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultFigureColor',[1 1 1])
path = convertCharsToStrings(fileparts(matlab.desktop.editor.getActiveFilename));
cd(path)
clear; clear global *; clc; warning off; close all;

%% COMSOL FILE INFORMATION
%-------------------------------------------------------------------------%
File.Path = [pwd,filesep,'Models'];
File.Tag = 'Comsol_QRD5';
File.Extension = '.mph';
%-------------------------------------------------------------------------%

%% FREQUENCY
%-------------------------------------------------------------------------%
Freq.f_min = 250;                                  % Minimum Freq of interest
Freq.f_max = 251;                               % Maximum Freq of interest
Freq.df = 1;                                    % Freq discretization
Freq.Vector = (Freq.f_min:Freq.df:Freq.f_max);    % Freq vector
Freq.Nf = numel(Freq.Vector);                   % Number of frequencies

%% COMSOL PROBE INFORMATION
%-------------------------------------------------------------------------%
Probe.radius = 3; %radius of arc in meters
Probe.theta_min = 0; %arc starting angle
Probe.theta_max = pi; %arc end angle
Probe.Resolution = 181;
Probe.theta_vector = linspace(Probe.theta_min,Probe.theta_max,Probe.Resolution);
Probe.Coordinates(1,:) = Probe.radius*cos(Probe.theta_vector); %Probe x coordinates
Probe.Coordinates(2,:) = Probe.radius*sin(Probe.theta_vector); %Probe y coordinates

%% FEM MODELLING
%-------------------------------------------------------------------------%
tStart = tic;
Ps_1 = QR_5(Freq,Probe,File); %call COMSOL model for QRD
Psflatnum = 
tEnd = toc(tStart);
fprintf('FEM. time: %d minutes and  %.f seconds\n', floor(tEnd/60), rem(tEnd,60));
%-------------------------------------------------------------------------%

%% CALCULATE DIFFUSION COEFFICIENT
%-------------------------------------------------------------------------%
SI_1 = abs(Ps_1).^2; %sound intensity
SIsum_1 = sum(SI_1,2);
SIsq_1 = sum(SI_1.^2,2);

n_d = length(Probe.theta_vector);
delta_COMSOL = (SIsum_1.^2 - SIsq_1)./((n_d-1)*(SIsq_1));

run("QRD_TMM.m")
plot(Freq.Vector,delta_COMSOL); %plot diffusion coefficient
title("diffusion coefficient of N=5 QRD")
hold on
plot(Freq.Vector,delta_QRD)
legend("numerical","TMM")
ylim([0, 1])
%-------------------------------------------------------------------------%
