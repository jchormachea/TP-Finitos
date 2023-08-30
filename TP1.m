% Case launcher
clc; clear; close all; 

%% Structure definition
exercise = 4; %input('Elegir estructura: '); %Cada caso es una estructura propuesta
onlyBars = false;
fprintf('Exercise %i.\n',exercise)

switch exercise
        case 1
            
            structuralJointsArray=[0, 0, 0
                                   0, 10.1, 0
                                   0, 15, 0
                                   5.5, 15, 0
                                   10.8, 12, 0
                                   5.5, 10.1, 0                                   
                                   0 5 0
                                   2, 5.3458, 0]*1000; 
        
                        % Input of members connecting joints and which dof they connect
                        % Begin | End | Cross Section Orientation
            structuralMembersArray.nodes=[7 6 8
                                          3 4 8
                                          4 5 8
                                          5 6 8
                                          2 6 8
                                          3 6 8
                                          4 6 8
                                          6 1 8];
            planeStructure = true;
            hollowBar = false;
        case 2 
            structuralJointsArray=[ 0 10 0          % In [mm]
                        0 15 0
                        6 10 0
                        11 13 0
                        2 3.1415 0]*1000; 
        
                        % Input of members connecting joints and which dof they connect
                        % Begin | End | Cross Section Orientation
            structuralMembersArray.nodes=[1 3 5
                                          2 3 5
                                          2 4 5
                                          3 4 5];
            planeStructure = true;
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
structuralMembersArray.dof = true(size(structuralMembersArray.nodes,1),12);

% Number of elements in member
structuralMembersArray.refinement = ones(size(structuralMembersArray.nodes,1));

% Member cross section number
structuralMembersArray.crossSection = ones(size(structuralMembersArray.nodes,1),4);

switch exercise
    case 4
        structuralMembersArray.crossSection(:,:) = 1;
        structuralMembersArray.crossSection(9,:) = 1;
        structuralMembersArray.crossSection(11,:) = 1;
end


% Member material number
structuralMembersArray.material = ones(size(structuralMembersArray.nodes,1),3);



%% Cross sections definition
% exercise 1 and 2, retangular hollow section

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


% Area | Inertia Moment in P123 plane | Inertia Moment orthogonal to P123 plane | Torsional Stiffness
membersCrossSection=[A1  Izz1 Iyy1 Tk1]; % mm2 mm4

%% Material definition
% Young Modulus | Transverse Modulus | Density  %MPa kg/m3
membersMaterial=[200000 80000 7800]; %Material del TP
alpha = 1.26e-5; %1/°C
Sy = 200; %Mpa
%% Structure plot
linearMeshPlot(structuralMembersArray.nodes(:,1:2),structuralJointsArray,'b','Yes');
             
%% Preprocess ---------------------------------------------------------------------------------------------
    
% Mesh generation
[elementArray,nodesPositionArray]=trussFrameMeshGenerator(structuralMembersArray,structuralJointsArray);

% Problem parameters
nElements=size(elementArray.nodes,1);    %Number of elements
nNodes=size(nodesPositionArray,1);       %Number of nodes
nTotalDof=max(max(elementArray.dof));    %Number of total dofs

dT = 60; %Temperature difference - The same for all cases

switch exercise 
    case 1

        % Boundary conditions
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
        boundaryConditionsArray([1,2,3,7],:) = true;
        % Rotations elimination
        % Load definition
        pointLoadsArray = zeros(nNodes,6);     % Point load nodal value for each direction
        pointLoadsArray(5,6) = 20000/(8*pi/30)*1000; %Nmm
        pointLoadsArray(5,2) = -600000-350*9.81; %N
        pointLoadsArray(5,4) = -0.3*(+600000+350*9.81)*1000; %Nmm
        
        
        thermaLoads = zeros(nNodes,6); %thermal loads  
        thermalElements = 0; %Vector with the elements with thermal load
        pointLoadsArray = pointLoadsArray+thermaLoads;

    case 2
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
        boundaryConditionsArray([1 2],:) = true;
        
          % Load definition
        pointLoadsArray = zeros(nNodes,6);     % Point load nodal value for each direction
        pointLoadsArray(4,6) = 20000/(8*pi/30)*1000; %Nmm
        pointLoadsArray(4,2) = -600000-350*9.81; %Nmm
        pointLoadsArray(4,4) = 0.3*(+600000+350*9.81)*1000; %Nmm
        
        thermaLoads = zeros(nNodes,6); %thermal loads  
        thermalElements = 0; %Vector with the elements with thermal load
        
        pointLoadsArray = pointLoadsArray+thermaLoads;
        
    case 3
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
        boundaryConditionsArray([1 2 3 7 8 9],[1 2 3]) = true; % Apoyos fijos
        
          % Load definition
        pointLoadsArray = zeros(nNodes,6);     % Point load nodal value for each direction
        pointLoadsArray(5,6) = -20000/(8*pi/30)*1000; %Nmm
        pointLoadsArray(5,2) = -600000-350*9.81; %N
        
        thermaLoads = zeros(nNodes,6); %thermal loads  
        thermalElements = [9 17 1 10 2]; %Vector with the elements with thermal load
        
        
        for iElement = thermalElements
            ftLocal = alpha*membersMaterial(1)*membersCrossSection(elementArray.crossSection(iElement),1)*dT;
            
            node1 = nodesPositionArray(elementArray.nodes(iElement,1),:);
            node2 = nodesPositionArray(elementArray.nodes(iElement,2),:);
            anode = structuralJointsArray(elementArray.auxiliarPoint(iElement),:);
            lambda = RotationMatrix(node1,node2,anode);
            
            Ft = lambda'*[ftLocal;0;0];
            thermaLoads(elementArray.nodes(iElement,1),[1 2 3]) = thermaLoads(elementArray.nodes(iElement,1),[1 2 3])-Ft';
            thermaLoads(elementArray.nodes(iElement,2),[1 2 3]) = thermaLoads(elementArray.nodes(iElement,2),[1 2 3])+Ft';
            
        end
        
        pointLoadsArray = pointLoadsArray+thermaLoads;
        
    case 4
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
% %         boundaryConditionsArray([1 2 5 6],[1 2 3]) = true; %Apoyos fijos
        boundaryConditionsArray([1 2 5 6],:) = true; %Empotramientos
        
          % Load definition
        pointLoadsArray = zeros(nNodes,6);     % Point load nodal value for each direction
        pointLoadsArray(4,6) = -20000/(8*pi/30)*1000; %Nmm
        pointLoadsArray(4,2) = -600000-350*9.81; %N
        
        thermaLoads = zeros(nNodes,6); %thermal loads 
        thermalElements = [3 7]; %Vector with the elements with thermal load    
        
        
        for iElement = thermalElements
            ftLocal = alpha*membersMaterial(1)*membersCrossSection(elementArray.crossSection(iElement),1)*dT; %N
            
            node1 = nodesPositionArray(elementArray.nodes(iElement,1),:);
            node2 = nodesPositionArray(elementArray.nodes(iElement,2),:);
            anode = structuralJointsArray(elementArray.auxiliarPoint(iElement),:); %auxiliar node
            lambda = RotationMatrix(node1,node2,anode);
            
            Ft = lambda'*[ftLocal;0;0];
            thermaLoads(elementArray.nodes(iElement,1),[1 2 3]) = thermaLoads(elementArray.nodes(iElement,1),[1 2 3])-Ft';
            thermaLoads(elementArray.nodes(iElement,2),[1 2 3]) = thermaLoads(elementArray.nodes(iElement,2),[1 2 3])+Ft';
            
        end
        
        pointLoadsArray = pointLoadsArray+thermaLoads;
        
end



if onlyBars
	 boundaryConditionsArray(:,[ 4 5 6]) = true; 
end

if planeStructure
    boundaryConditionsArray(:,[3 6]) = true; 
end

%% Solver

% Stiffness calculation and assembly
[stiffnessMatrix]=assemble1DStiffnessMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial);

% Matrix reduction
isFixed = reshape(boundaryConditionsArray',1,[])';
isFree = ~isFixed;

% Loads vector rearrangement
loadsVector = reshape(pointLoadsArray',1,[])';

% Equation solving
displacementsReducedVector = stiffnessMatrix(isFree,isFree)\loadsVector(isFree);

% Reconstruction
displacementsVector = zeros(nTotalDof,1);
displacementsVector(isFree) = displacementsVector(isFree) + displacementsReducedVector;

%% Postprocess
magnificationScale=10;

% Nodal displacements rearranged
nodalDisplacements=reshape(displacementsVector,6,size(nodesPositionArray,1))';
nodalPositions=nodesPositionArray+nodalDisplacements(:,1:3)*magnificationScale;

% Deformed Structure plot
linearMeshPlot(elementArray.nodes(:,1:2),nodalPositions,'r','No');
axis([0 11400 0 15500])

% Reactions calculation

reactionsVector = stiffnessMatrix*displacementsVector;
nodalReactions = reshape(reactionsVector, 6, size(nodesPositionArray,1))';

%Stresses calculation

nElements=size(elementArray.nodes,1);      %Number of elements
for iElement = 1:nElements
    
    node1 = nodesPositionArray(elementArray.nodes(iElement,1),:);
    node2 = nodesPositionArray(elementArray.nodes(iElement,2),:);
    L = norm(node2-node1);
    anode = structuralJointsArray(elementArray.auxiliarPoint(iElement),:);
    lambdaProjectionMatrix = RotationMatrix(node1,node2,anode);
    
    projectionMatrix = blkdiag(lambdaProjectionMatrix,lambdaProjectionMatrix);
    localDisplacements(1,:) = projectionMatrix* nodalDisplacements(elementArray.nodes(iElement,1),:)';
    localDisplacements(2,:) = projectionMatrix* nodalDisplacements(elementArray.nodes(iElement,2),:)';
    axialStrain = 1/L*[-1 1]*localDisplacements([1 2],1);
    
    axialStress = membersMaterial(1)*axialStrain;
    thermal = 0;
    if sum(iElement==thermalElements) > 0
        thermal = 1;
    end
    ftLocal = alpha*membersMaterial(1)*membersCrossSection(elementArray.crossSection(iElement),1)*dT;
    soAxial = -ftLocal/membersCrossSection(elementArray.crossSection(iElement),1)*thermal; %0 si el elemento no tiene carga termica
    axialStress = membersMaterial(1)*axialStrain+soAxial;
    
   if onlyBars 
       bendingStress = 0;
       tau_circular = 0;
   else
        A = membersCrossSection(elementArray.crossSection(iElement),1);
        Izz= membersCrossSection(elementArray.crossSection(iElement),2);
        Iyy = membersCrossSection(elementArray.crossSection(iElement),3);
        Tk= membersCrossSection(elementArray.crossSection(iElement),4);
        
        x = 0;    
        Bbending = [-6/L^2+12*x/L^3, -4/L+6*x/L^2, 6/L^2-12*x/L^3, -2/L+6*x/L^2];
        
       
        Mz = membersMaterial(1)*Izz*Bbending*[localDisplacements(1,2), localDisplacements(1,6), localDisplacements(2,2), localDisplacements(2,6)]';
        My = membersMaterial(1)*Iyy*Bbending*[localDisplacements(1,3), -localDisplacements(1,5), localDisplacements(2,3), -localDisplacements(2,5)]';
        
        M = sqrt(My^2+Mz^2); %SOLO PARA CIRCULAR 
        bendingStress = M*r(elementArray.crossSection(iElement))/Izz;

        %Torsion
        %tau_rectangular = 3*T/(8*a*b^2)*(1+0.6095*(b/a)+0.8865*(b/a)^2-1.8023*(b/a)^3+0.91*(b/a)^4);
        T = membersMaterial(2)*Tk*(localDisplacements(2,4)-localDisplacements(1,4))/L;
        ct = r(elementArray.crossSection(iElement));
        tau_circular = T*ct/Tk; %Maximo en el radio
   end
    elementArray.stress(iElement,1) = axialStress+bendingStress*sign(axialStress); %max tension in modulo 
    elementArray.stress(iElement,2) = tau_circular;
    elementArray.stress(iElement,3) = sqrt(elementArray.stress(iElement,1)^2+3*elementArray.stress(iElement,2)^2); %Von Mises
    
    % Structure Weight
    volume = L*membersCrossSection(elementArray.crossSection(iElement),1)/(1000*1000*1000); %m3
    weight = volume*membersMaterial(elementArray.material(iElement),3);%kg
    elementArray.weight(iElement) = weight; %kg
    
end

%total structure weight

totalWeigth = sum(elementArray.weight)/1000; %tons

%information
fprintf('only bars = %s and plane structure = %s .\n',string(onlyBars),string(planeStructure));

staticSF = Sy*elementArray.stress(:,3).^-1;
fprintf('Von Mises stress and SF for element:\n');
disp([[1:nElements]', elementArray.stress(:,3), staticSF])
fprintf('Nodal Displacements (mm/rad):\n');
disp(nodalDisplacements);
fprintf('Total Structure Weight(tons):\n');
disp(totalWeigth);






