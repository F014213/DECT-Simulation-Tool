%Dect main v2:

%% IMAGE BASED APPROACH  
close all; 
clear all;
clc;

%define two photon energies
E1 = 47e-3; %MeV (70 keV)
E2 = 100e-3; %MeV (150 keV)
acquisition_angle_res = 1; %degrees
phantom = 'acr';


DECT_benchmark = addition_benchmark_DECT(E1, E2, acquisition_angle_res, phantom);
DECT_air = subtraction_DECT(E1, E2, acquisition_angle_res, phantom, 'air', 'air');
DECT_bone = subtraction_DECT(E1, E2, acquisition_angle_res, phantom, 'bone', 'bone');
DECT_poly = subtraction_DECT(E1, E2, acquisition_angle_res, phantom, 'polyethylene', 'polyethylene');
DECT_acryl = subtraction_DECT(E1, E2, acquisition_angle_res, phantom, 'acrylic', 'acrylic');

% t1 = DECT_bone;
% t1(t1>0.3) = 0;
% %figure, imshow(t1)
% t2 = t1/max(t1,[],'all');
% %figure,imshow(t2)
% 
% %R = t1;


R = DECT_benchmark;
G = DECT_air;
B = DECT_bone;

RGB = cat(3, R, G, -(B-1) );

figure, imshow(RGB)
