% Case launcher
clear; close all; 

%% Structure definition  

exercise = 2; %input('Elegir estructura: '); %Cada caso es una estructura propuesta
onlyBars = false;
fprintf('Exercise %i.\n',exercise)

switch exercise
        case 1
            structuralJointsArray=[0 6 0% in[mm]
                               0 4 0
                               0 6 5
                               0 4 5
                               3 5 2
                               6 5 2
                               11 5.3 2
                               0 6 -3
                               0 4 -3
                               6 6 5
                               6 6 0
                               6 4 5
                               6 4 0
                               6 4 2
                        2.55 2.55 .55]*1000; % for direction only
                    aux=length(structuralJointsArray);
            % Input of members connecting joints and which dof they connect
            % Begin | End | Cross Section Orientation
            structuralMembersArray.nodes=[1 5 aux
                                      2 5 aux
                                      3 5 aux
                                      4 5 aux
                                      1 6 aux
                                      2 6 aux
                                      3 6 aux
                                      4 6 aux
                                      5 6 aux
                                      1 11 aux
                                      2 13 aux
                                      3 10 aux
                                      4 12 aux
                                      5 10 aux
                                      5 11 aux
                                      6 7 aux
                                      11 7 aux
                                      10 7 aux
                                      8 11 aux
                                      9 13 aux
                                      12 7 aux
                                      13 7 aux
                                      13 11 aux
                                      5 12 aux
                                      5 13 aux
                                      10 12 aux];
            planeStructure = false;
            hollowBar = false;
        case 2
            structuralJointsArray=[ 0 10 0          % [mm]
                        0 15 0
                        5.8 10.3 0
                        11, 12, 1.3
                        11,12,1
                        0 10 2
                        0 15 2
                        5.8 10.3 2
                        2, 5.3458, 7.6589]*1000;
            aux = length(structuralJointsArray);          
            structuralMembersArray.nodes=[1 3 aux
                                          2 3 aux
                                          2 5 aux %thermal
                                          3 5 aux
                                          6 8 aux
                                          7 8 aux
                                          7 4 aux %thermal
                                          8 4 aux
                                          1 8 aux
                                          3 8 aux
                                          7 3 aux
                                          4 5 aux
                                          3 4 aux];
                                          
            planeStructure = false;
            hollowBar = false;
        
        case 3
            structuralJointsArray=[0 5 0
                                   0, 10.1, 0
                                   0, 15, 0
                                   5.5, 15, 0
                                   11, 12, 1.3
                                   5.5, 10.1, 0
                                   0 5 2
                                   0, 10.1, 2
                                   0, 15, 2
                                   5.5, 15, 2
                                   5.5, 10.1, 2 
                                   2, 5.3458, 7.6589]*1000;
            aux = length(structuralJointsArray);          
            structuralMembersArray.nodes=[3 4 aux %thermal
                                          4 5 aux %thermal
                                          5 6 aux
                                          5 6 aux
                                          6 1 aux
                                          2 6 aux
                                          3 6 aux
                                          4 6 aux
                                          9 10 aux %thermal
                                          10 5 aux %thermal
                                          10 6 aux
                                          5 11 aux
                                          11 7 aux
                                          8 11 aux
                                          9 11 aux
                                          10 11 aux
                                          9 4 aux %thermal
                                          2 10 aux
                                          7 6 aux];
            planeStructure = false;
            hollowBar = false;
        case 4
            structuralJointsArray=[ 0 10 0          % [mm]
                        0 15 0
                        5.8 10.3 0
                        11, 12, 1.3
                        0 10 2
                        0 15 2
                        5.8 10.3 2
                        2, 5.3458, 7.6589]*1000;
            aux = length(structuralJointsArray);          
            structuralMembersArray.nodes=[1 3 aux
                                          2 3 aux
                                          2 4 aux %thermal
                                          3 4 aux
                                          5 7 aux
                                          6 7 aux
                                          6 4 aux %thermal
                                          7 4 aux
                                          1 7 aux
                                          3 7 aux
                                          6 3 aux];
                                          
            planeStructure = false;
            hollowBar = false;
        

end

                          
% Connected Dof                          
structuralMembersArray.dof=true(size(structuralMembersArray.nodes,1),12);

% Number of elements in member
structuralMembersArray.refinement=ones(size(structuralMembersArray.nodes,1));

% Member cross section number
structuralMembersArray.crossSection=ones(size(structuralMembersArray.nodes,1));

% Member material number
structuralMembersArray.material=ones(size(structuralMembersArray.nodes,1));
                    
if hollowBar %circular bar in this case
    a = 60; %mm
    b = 40; %mm
    t = 5; %mm
    Izz1 = a*b^3/12-(a-2*t)*(b-2*t)^3/12;%mm4
    Iyy1 = a^3*b/12-(a-2*t)^3*(b-2*t)/12;%mm4
    A1 = b*a-(b-2*t)*(a-2*t);
    a = a/2; b = b/2; %para la formula de Tk
    Tk1 = a*b^3*(16/3-3.36*b/a*(1-b^4/(12*a^4)))-(a-t)*(b-t)^3*(16/3-3.36*(b-t)/(a-t)*(1-(b-t)^4/(12*(a-t)^4)));%mm4 (2*t^2*(a-t)^2*(b-t)^2)/(a*t+b*t-t^2-t^2)
else
    %circular bar
    Amax = 0.25*0.25*1000*1000; %Area maxima en mm2 
% %     r = sqrt(Amax/pi)*[sqrt(0.7) sqrt(0.13) sqrt(0.07) sqrt(0.05) sqrt(0.00001)]';%mm
    r = 60; %mm
    Izz1 = pi*r.^4/4;%mm4
    Iyy1 = Izz1;%mm4
    A1 = pi*r.^2;
    Tk1 = Izz1*2;
end

% Cross sections definition
% Area | Inertia Moment in P123 plane | Inertia Moment orthogonal to P123 plane | Torsional Stiffness
% membersCrossSection=[6980 16011520 108608400 278528]; % I Beam Section
% membersCrossSection=[6980 16011520 108608400 278528]; % I Beam Section
membersCrossSection=[A1  Izz1 Iyy1 Tk1]; % mm2 mm4

% Material definition
% Young Modulus | Transverse Modulus | Density 
% % membersMaterial=[68980 26531 2.7e-9]; %MPa kg/mm3/1000 Aluminium 6061
membersMaterial=[200000 80000 7800/1000^3/1000]; %Material del TP

% Structure plot
linearMeshPlot(structuralMembersArray.nodes(:,1:2),structuralJointsArray,'b','Yes');
             
%% Preprocess

% Mesh generation
[elementArray,nodesPositionArray]=trussFrameMeshGenerator(structuralMembersArray,structuralJointsArray);

% Problem parameters
nElements=size(elementArray.nodes,1);    %Number of elements
nNodes=size(nodesPositionArray,1);       %Number of nodes
nTotalDof=max(max(elementArray.dof));    %Number of total dofs

% Boundary conditions
switch exercise 
    case 1
        % Boundary conditions
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
        boundaryConditionsArray([1 2 3 4 8 9],:) = true;
        massNode = 7;

    case 2
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
% %         boundaryConditionsArray([1 2 5 6],[1 2 3]) = true; %Apoyos fijos
        boundaryConditionsArray([1 2 6 7],:) = true; %Empotramientos
        massNode = 4;
        
        
    case 3
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
        boundaryConditionsArray([1 2 3 7 8 9],[1 2 3]) = true; % Apoyos fijos
        massNode = 5;
    case 4
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
% %         boundaryConditionsArray([1 2 5 6],[1 2 3]) = true; %Apoyos fijos
        boundaryConditionsArray([1 2 5 6],:) = true; %Empotramientos
        massNode = 4;
        
end



        

%% Resonant modes determination

% Stiffness calculation and assembly
[stiffnessMatrix]=assemble1DStiffnessMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial);

% Stiffness calculation and assembly
[massMatrix]=assemble1DMassMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial);


sz = size(massMatrix);
puntualMassMatrix(sz(1),sz(2)) = 0;
puntualMass = 350/1000 ;%Kg/1000
puntualMassMatrix(6*(massNode-1)+1,6*(massNode-1)+1) = puntualMassMatrix(6*(massNode-1)+1,6*(massNode-1)+1)+ puntualMass;
puntualMassMatrix(6*(massNode-1)+2,6*(massNode-1)+2) = puntualMassMatrix(6*(massNode-1)+2,6*(massNode-1)+2)+ puntualMass;
puntualMassMatrix(6*(massNode-1)+3,6*(massNode-1)+3) = puntualMassMatrix(6*(massNode-1)+3,6*(massNode-1)+3)+ puntualMass;

massMatrix = massMatrix+puntualMassMatrix;

% Matrix reduction
isFixed = reshape(boundaryConditionsArray',1,[])';
isFree = ~isFixed;

% Eigenvalues determination
[eigenShapes, eigenAngularSquaredFrequencies] = eig(stiffnessMatrix(isFree,isFree),massMatrix(isFree,isFree));

% Eigenmode sorting and scaling
[eigenAngularSquaredFrequencies,eigenModesOrder]=sort(diag(eigenAngularSquaredFrequencies));
nModes=size(eigenAngularSquaredFrequencies,1)
naturalFrequencies=sqrt(eigenAngularSquaredFrequencies)/2/pi;

eigenShapes=eigenShapes(:,eigenModesOrder);

% Scaling

for iMode=1:nModes
    eigenShapes(:,iMode) = eigenShapes(:,iMode)/sqrt(eigenShapes(:,iMode)'*massMatrix(isFree,isFree)*eigenShapes(:,iMode));
end

% Recovery
displacementsEigenShapes = zeros(nTotalDof,nModes);
displacementsEigenShapes(isFree,:) = displacementsEigenShapes(isFree,:) + eigenShapes(:,:);

%% Postprocess
magnificationScale=100;

% First natural mode plot
nodalDisplacements=reshape(displacementsEigenShapes(:,1),6,size(nodesPositionArray,1))';
nodalPositions=nodesPositionArray+nodalDisplacements(:,1:3)*magnificationScale;
linearMeshPlot(elementArray.nodes(:,1:2),nodalPositions,'g','No');

% Second natural mode plot
nodalDisplacements=reshape(displacementsEigenShapes(:,2),6,size(nodesPositionArray,1))';
nodalPositions=nodesPositionArray+nodalDisplacements(:,1:3)*magnificationScale;
linearMeshPlot(elementArray.nodes(:,1:2),nodalPositions,'m','No');

% Third natural mode plot
nodalDisplacements=reshape(displacementsEigenShapes(:,3),6,size(nodesPositionArray,1))';
nodalPositions=nodesPositionArray+nodalDisplacements(:,1:3)*magnificationScale;
linearMeshPlot(elementArray.nodes(:,1:2),nodalPositions,'c','No');
title({['1st Natural Frequency   ' num2str(naturalFrequencies(1))], ['2nd Natural Frequency ' num2str(naturalFrequencies(2))], ['3rd Natural Frequency ' num2str(naturalFrequencies(3))]});
view ([1 1 1]);

% % break

%% Auxiliar calculations
E=68980; %[MPa];
rho=2700e-9; %[kg/mm3];
A=200*2*12+(294-2*12)*8; %[mm2]
Imax=294^3*200/12-(294-2*12)^3*(200-8)/12; %[mm4]
Imin=200^3*2*12/12+(294-2*12)*8^3/12; %[mm4]
L=2000; %[mm]
f1=1.875^2/pi/2/L^2*sqrt(E*Imin/rho/A*1e3);
f2=1.875^2/pi/2/L^2*sqrt(E*Imax/rho/A*1e3);