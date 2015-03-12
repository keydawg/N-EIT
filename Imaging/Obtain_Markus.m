% HAVE TO RUN FROM THE SAME FOLDER
addpath( genpath('/home/kdawg/Desktop/PEITS/dune-peits/matlab'));
addpath( genpath('/home/kdawg/Desktop/Rat_forward_linux/'));
%filepath=pwd;


%load('Prt_electrodes.mat')
%load('Mesh_hex_CROSSVAL_100um.mat'); % Change here
%load('Rat_027-Prt_ALL_for_inversion');
%prt_electrodes=dlmread('Rat_031-Prt_ALL_for_inversion.txt'); %here
direc='Markus'; %here
load('Mesh.mat');


files=dir([direc '/electrodevoltages*']);
files=files.name;
BV0 = load_electrode_voltages_binary([direc '/' files]);
disp('Jacobian cut shit...');

files=dir([direc '/sigmavector*']);
sigma_file=[direc '/' files.name];

files=dir([direc '/jacobian*']);
jac_file=[direc '/' files.name];

J = load_jacobian_binary_ID(jac_file, sigma_file, 1:length(Mesh.Tetra));
 

save('Jacobian_Tet_GSK_P9.mat','J','BV0','-v7.3');
sens=sum(J.^2,1).^0.5;
save('sens_GSK.mat','sens','-v7.3');
