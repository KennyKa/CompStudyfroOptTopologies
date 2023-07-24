%------------------------------- PolyStress ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
clear; clc; close all
restoredefaultpath; addpath(genpath('./')); %Use all folders and subfolders
set(0,'defaulttextinterpreter','latex')
%% ------------------------------------------------------------ CREATE Mesh

sajatvagynem = 0

fileName = "Mesh_Eyebar2";
dataName = str2func(fileName)

LoadMultiplierX = 0.9
LoadMultiplierY = 0.8

strLim = 450

haloszam = 10000

if sajatvagynem == 1
    haloszam = "NA"
    NUMBER = 2
    
    nana = load(fileName + "_nodes.txt");
    %elel = readcell("PolyProb1_elements.txt",Delimiter=";")
    elel = load(fileName + "_elements.txt")
    num2cell(elel, 2)
    x = nana(:, 1);
    y = nana(:, 2);
    
    Node = [x, y];
    Element = num2cell(elel, 2);
    class(Element);
    length(Element);
    
    Supp = load(fileName + "_supps.txt");
    
    Load = load(fileName + "_loads.txt")
    
    Load(:,2) = Load(:,2) * LoadMultiplierX;
    Load(:,3) = Load(:,3) * LoadMultiplierY;
    Nx = sum(Load(:,2))
    Ny = sum(Load(:,3))
    figure(2)
    PolyMshr_PlotMsh(Node, Element, length(Element), Supp, Load)
    hold on
else
%%

    [Node,Element,Supp,Load] = dataName(haloszam); 
    Load(:,2) = Load(:,2) * LoadMultiplierX;
    Load(:,3) = Load(:,3) * LoadMultiplierY;
    Nx = sum(Load(:,2))
    Ny = sum(Load(:,3))
end

NElem = size(Element,1); % Number of elements
%% ---------------------------------------------------- CREATE 'fem' STRUCT
%E0 = 70e3; % E0 in MPa
E0 = 200e3; % E0 in MPa
G = E0/2.5; Et = E0; Ec = E0;  % 0<=(Et,Ec)<=3*G; %Material props. (linear)
fem = struct(...
  'NNode',size(Node,1),...      % Number of nodes
  'NElem',size(Element,1),...   % Number of elements
  'Node',Node,...               % [NNode x 2] array of nodes
  'Element',{Element},...       % [NElement x Var] cell array of elements
  'Supp',Supp,...               % Array of supports
  'Load',Load,...               % Array of loads
  'Passive',[],...              % Passive elements  
  'Thickness',0.3,...             % Element thickness
  'MatModel','Bilinear',...     % Material model ('Bilinear','Polynomial')
  'MatParam',[Et,Ec,G],...      % Material parameters for MatModel
  'SLim',strLim,...                % Stress limit
  'TolR', 1e-8, ...             % Tolerance for norm of force residual
  'MaxIter', 355, ...            % Max NR iterations per load step
  'MEX', 'No');                 % Tag to use MEX functions in NLFEM routine
%% ---------------------------------------------------- CREATE 'opt' STRUCT
R = 0.04; q = 3; % Filter radius and filter exponent
p = 4; eta0 = 0.5;
m = @(y,B)MatIntFnc(y,'SIMP-H1',[p,B,eta0]);
%[P, ElCenter] = PolyFilter(fem,R,q);
%[P, ElCenter] = PolyFilter(fem,R,q, 'X');
[P, ElCenter] = PolyFilter(fem,R,q, 'Y');
zIni = 0.5*ones(size(P,2),1);
opt = struct(...               
  'zMin',0.0,...              % Lower bound for design variables
  'zMax',1.0,...              % Upper bound for design variables
  'zIni',zIni,...             % Initial design variables
  'MatIntFnc',m,...           % Handle to material interpolation fnc.
  'contB',[5,1,1,10],...      % Threshold projection continuation params.  
  'P',P,...                   % Matrix that maps design to element vars.
  'Tol',0.002,...             % Convergence tolerance on design vars.
  'TolS',0.003,...            % Convergence tolerance on stress constraints
  'MaxIter',150,...           % Maximum number of AL steps
  'MMA_Iter',5,...            % Number of MMA iterations per AL step
  'lambda0',zeros(NElem,1),...% Initial Lagrange multiplier estimators
  'mu0',10,...                % Initial penalty factor for AL function
  'mu_max',10000,...          % Maximum penalty factor for AL function
  'alpha',1.1,...             % Penalty factor update parameter  
  'Move',0.15,...             % Allowable move step in MMA update scheme
  'Osc',0.2,...               % Osc parameter in MMA update scheme
  'AsymInit',0.2,...          % Initial asymptote in MMA update shecme
  'AsymInc',1.2,...           % Asymptote increment in MMA update scheme  
  'AsymDecr',0.7...           % Asymptote decrement in MMA update scheme     
   );
%% ------------------------------------------------------- RUN 'PolyStress'
fem = preComputations(fem); % Run preComputations before running PolyStress
[z,V,fem, elapsedTime] = PolyStress(fem,opt);
% [z,V,fem, Frames] = PolyStressAvi(fem,opt);
%-------------------------------------------------------------------------%
%%

elapsedTime = round(elapsedTime);

kepnev = fileName + "_R=" + R + "_Nx=" + Nx + "_Ny=" + Ny +...
    "_sm=" + haloszam + "_t=" + elapsedTime + "s" + "_E0=" + E0 + "Et=" +...
    Et + "Ec=" + Ec + "G0=" + G + "SL=" + strLim

f = figure(3);
toprint = f.Children(1:2);
toprint2 = f.Children(3);

fig4 = figure(4);
set(gcf, 'Position', get(0, 'Screensize'));
copyobj(toprint, fig4);
f = gcf;
h = get(f,'Children');
%exportgraphics(f,'./DiplomaKepek/' + kepnev + '_vonM.pdf','Resolution',300);
exportgraphics(f,'./DiplomaKepek/' + kepnev + '_vonM.png','Resolution',300);

fig5 = figure(5);
set(gcf, 'Position', get(0, 'Screensize'));
copyobj(toprint2, fig5);
f = gcf;
h = get(f,'Children');
%exportgraphics(f,'./DiplomaKepek/' + kepnev + '_Dens.pdf','Resolution',300);
exportgraphics(f,'./DiplomaKepek/' + kepnev + '_Dens.png','Resolution',300);

close all;

return
Image = getframe(gcf);
imwrite(Image.cdata, 'mask_image.tiff');

return
video = VideoWriter('Proba.avi', 'Uncompressed AVI');
open(video);
writeVideo(video, Frames);
close(video);

return
kellindex = find(V>0.6);

ElCenterKell = ElCenter(kellindex,:);
figure(3);
plot(ElCenterKell(:,1), ElCenterKell(:,2), 'bd');
axis image;

fileID = fopen('./SajatFajlok/_OptimaltTextek/' + fileName + '_' + NUMBER + '.txt','w');
%fprintf(fileID,'%6s %12s\n','x', 'y');
fprintf(fileID,'%8.4f %8.4f\n',ElCenterKell');
fclose(fileID);