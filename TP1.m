% Case launcher
clc; clear; close all; 

%% Structure definition
exercise = 1; %input('Elegir estructura: '); %Cada caso es una estructura propuesta
%caso 1 es la estructura optimizando
%caso 2 es la estrucutra orginal que verifica
%caso 3 es una estructura vieja que no aguanta
%caso 4 es la estructura original pero con un solo nodo en la punta
onlyBars = false;
fprintf('Exercise %i.\n',exercise)

switch exercise
        case 1            
            structuralJointsArray=[ 0 10 0          % [mm]
                        0 15 0
                        5.8 10.3 0
                        11, 12, 0.7
                        11,12,1
                        0 10 2
                        0 15 2
                        5.8 10.3 2
                        0 8 1
                        2, 5.3458, 7.6589
                        0,10.3,1]*1000;
            aux = length(structuralJointsArray)-1; %anteultimo
            aux2 = length(structuralJointsArray); %ultimo
            structuralMembersArray.nodes=[1 3 2
                                          2 3 1
                                          2 4 aux %thermal
                                          3 4 aux
                                          6 8 7
                                          7 8 6
                                          7 5 2 %thermal
                                          8 5 aux
                                          1 8 aux
                                          3 8 aux2
                                          7 3 aux
                                          4 5 aux %thermal
                                          5 3 aux
                                          8 9 aux
                                          3 9 aux
                                          ];
                                          
            planeStructure = false;
        case 2
            structuralJointsArray=[ 0 10 0          % [mm]
                        0 15 0
                        5.8 10.3 0
                        11, 12, 0.7
                        11,12,1
                        0 10 2
                        0 15 2
                        5.8 10.3 2
                        2, 5.3458, 7.6589]*1000;
            aux = length(structuralJointsArray);          
            structuralMembersArray.nodes=[1 3 aux
                                          2 3 aux
                                          2 4 aux %thermal
                                          3 4 aux
                                          6 8 aux
                                          7 8 aux
                                          7 5 aux %thermal
                                          8 5 aux
                                          1 8 aux
                                          3 8 aux
                                          7 3 aux
                                          4 5 aux
                                          3 5 aux];
                                          
            planeStructure = false;
            
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
        

end

% Connected Dof                          
structuralMembersArray.dof = true(size(structuralMembersArray.nodes,1),12);

% Number of elements in member
structuralMembersArray.refinement = ones(size(structuralMembersArray.nodes,1));

% Member cross section number
structuralMembersArray.crossSection = ones(size(structuralMembersArray.nodes,1),4);



% Member material number
structuralMembersArray.material = ones(size(structuralMembersArray.nodes,1),3);



%% Cross sections definition
%defino diferentes secciones
Amax = 0.25*0.25*1000*1000; %Area maxima en mm2 --> rmax = 140mm

%circular bar
r = 30; %mm
Izr60 = pi*r.^4/4;%mm4
Iyr60 = Izr60;%mm4
Ar60 = pi*r.^2;
Tkr60 = Izr60*2;
%circular bar r80
r = 80; %mm
Izr80 = pi*r.^4/4;%mm4
Iyr80 = Izr80;%mm4
Ar80 = pi*r.^2;
Tkr80 = Izr80*2;

%circular bar r130 t10
ro = 130; %mm
t = 10;
ri = ro-t; %mm
Izr130t10 = pi*(ro^4-ri^4)/4;%mm4
Iyr130t10 = Izr130t10;%mm4
Ar130t10 = pi*(ro^2-ri^2);
Tkr130t10 = Izr130t10*2;

% Seccion cuadrada l20 t15
t = 15; %mm
a = 120; %mm
b = 120; %mm
IzC120t15 = a*b^3/12-(a-2*t)*(b-2*t)^3/12;%mm4
IyC120t15 = a^3*b/12-(a-2*t)^3*(b-2*t)/12;%mm4
AC120t15 = b*a-(b-2*t)*(a-2*t);
a = a/2; b = b/2; %para la formula de Tk
TkC120t15 = a*b^3*(16/3-3.36*b/a*(1-b^4/(12*a^4)))-(a-t)*(b-t)^3*(16/3-3.36*(b-t)/(a-t)*(1-(b-t)^4/(12*(a-t)^4)));%mm4 (2*t^2*(a-t)^2*(b-t)^2)/(a*t+b*t-t^2-t^2)
TkC120t15 = 18702997.496897; %mm4 según Nx

% Seccion cuadrada l20 t20
t = 20; %mm
a = 120; %mm
b = 120; %mm
IzC120t20 = a*b^3/12-(a-2*t)*(b-2*t)^3/12;%mm4
IyC120t20 = a^3*b/12-(a-2*t)^3*(b-2*t)/12;%mm4
AC120t20 = b*a-(b-2*t)*(a-2*t);
a = a/2; b = b/2; %para la formula de Tk
TkC120t20 = a*b^3*(16/3-3.36*b/a*(1-b^4/(12*a^4)))-(a-t)*(b-t)^3*(16/3-3.36*(b-t)/(a-t)*(1-(b-t)^4/(12*(a-t)^4)));%mm4 (2*t^2*(a-t)^2*(b-t)^2)/(a*t+b*t-t^2-t^2)
TkC110t20 = 22317010.477737; %mm4 según Nx


% Area | Inertia Moment in P123 plane | Inertia Moment orthogonal to P123 plane | Torsional Stiffness
membersCrossSection=[Ar60  Izr60 Iyr60 Tkr60
                     AC120t15 IzC120t15 IyC120t15 TkC120t15
                     Ar80  Izr80 Iyr80 Tkr80
                     Ar130t10  Izr130t10 Iyr130t10 Tkr130t10
                     AC120t20 IzC120t20 IyC120t20 TkC120t20
                     ]; % mm2 mm4


structuralMembersArray.circular = ones(size(structuralMembersArray.nodes,1),1); %para diferenciar secciones al momento de calcular tensiones
switch exercise
    case 1
        structuralMembersArray.crossSection(:,:) = 3; %circular r80        
        structuralMembersArray.crossSection(12,:) = 4; %circular r130 t10
        structuralMembersArray.crossSection([9 11],:) = 1; %circular r60
        
        structuralMembersArray.crossSection([2 6 10],:) = 2; %cuadrada 120t15
        structuralMembersArray.circular([2 6 10]) = 2; %cuadrada 120t15
                
        structuralMembersArray.crossSection([1 5],:) = 5; %cuadrada 120t20
        structuralMembersArray.circular([1 5]) = 2; %cuadrada 120t20
        
        
end   
                 
                 
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

%structure weight
for i = 1:nElements
    
    n1 = nodesPositionArray(elementArray.nodes(i,1),:);
    n2 = nodesPositionArray(elementArray.nodes(i,2),:);
    L = norm(n2-n1);
    % Structure Weight
    volume = L*membersCrossSection(elementArray.crossSection(i),1)/(1000*1000*1000); %m3
    weight = volume*membersMaterial(elementArray.material(i),3);%kg
    elementArray.weight(i) = weight; %kg
    
    
end

%Boundary conditions and Load Cases

dT = 0; %Temperature difference - The same for all cases
g = 9.81; %m/s^2
acceleration = [0,-1,0.2]*g;%m/s^2

switch exercise 
    case 1

        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
% %         boundaryConditionsArray([1 2 5 6],[1 2 3]) = true; %Apoyos fijos
        boundaryConditionsArray([1 2 6 7 9],:) = true; %Empotramientos
        
          % Load definition
        pointLoadsArray = zeros(nNodes,6);     % Point load nodal value for each direction
        pointLoadsArray(5,6) = -20000/(8*pi/30)*1000; %Nmm
        pointLoadsArray(5,2) = -600000; %N
        pointLoadsArray(4,2) = -350*g; %N
        massNode = 4; %donde está el motor
        
        
        for i = 1:nElements
            fAccelerationLocal = elementArray.weight(i)*acceleration/2;
            
            fAcceleration = fAccelerationLocal;
            pointLoadsArray(elementArray.nodes(i,1),[1 2 3]) = pointLoadsArray(elementArray.nodes(i,1),[1 2 3]) + fAcceleration ; %N
            pointLoadsArray(elementArray.nodes(i,2),[1 2 3]) = pointLoadsArray(elementArray.nodes(i,2),[1 2 3]) + fAcceleration; %N
        end
        
        thermaLoads = zeros(nNodes,6); %thermal loads 
        thermalElements = [3 7 12]; %Vector with the elements with thermal load    
        
        
        for iElement = thermalElements
            ftLocal = alpha*membersMaterial(1)*membersCrossSection(elementArray.crossSection(iElement),1)*dT; %N
            node1 = nodesPositionArray(elementArray.nodes(iElement,1),:);
            node2 = nodesPositionArray(elementArray.nodes(iElement,2),:);
            anode = structuralJointsArray(elementArray.auxiliarPoint(iElement),:);    
            lambda = RotationMatrix(node1,node2,anode);
            
            Ft = lambda'*[ftLocal;0;0];
            thermaLoads(elementArray.nodes(iElement,1),[1 2 3]) = thermaLoads(elementArray.nodes(iElement,1),[1 2 3])-Ft';
            thermaLoads(elementArray.nodes(iElement,2),[1 2 3]) = thermaLoads(elementArray.nodes(iElement,2),[1 2 3])+Ft';
            
        end
        
        pointLoadsArray = pointLoadsArray+thermaLoads;

    case 2
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
% %         boundaryConditionsArray([1 2 5 6],[1 2 3]) = true; %Apoyos fijos
        boundaryConditionsArray([1 2 6 7],:) = true; %Empotramientos
        
          % Load definition
        pointLoadsArray = zeros(nNodes,6);     % Point load nodal value for each direction
        pointLoadsArray(5,6) = -20000/(8*pi/30)*1000; %Nmm
        pointLoadsArray(5,2) = -600000; %N
        pointLoadsArray(4,2) = -350*g; %N
        massNode = 4; %donde está el motor
        
        
        for i = 1:nElements
            fAccelerationLocal = elementArray.weight(i)*acceleration/2;
            
            fAcceleration = fAccelerationLocal;
            pointLoadsArray(elementArray.nodes(i,1),[1 2 3]) = pointLoadsArray(elementArray.nodes(i,1),[1 2 3]) + fAcceleration ; %N
            pointLoadsArray(elementArray.nodes(i,2),[1 2 3]) = pointLoadsArray(elementArray.nodes(i,2),[1 2 3]) + fAcceleration; %N
        end
        
        thermaLoads = zeros(nNodes,6); %thermal loads 
        thermalElements = [3 7 12]; %Vector with the elements with thermal load    
        
        
        for iElement = thermalElements
            ftLocal = alpha*membersMaterial(1)*membersCrossSection(elementArray.crossSection(iElement),1)*dT; %N
            node1 = nodesPositionArray(elementArray.nodes(iElement,1),:);
            node2 = nodesPositionArray(elementArray.nodes(iElement,2),:);
            anode = structuralJointsArray(elementArray.auxiliarPoint(iElement),:);    
            lambda = RotationMatrix(node1,node2,anode);
            
            Ft = lambda'*[ftLocal;0;0];
            thermaLoads(elementArray.nodes(iElement,1),[1 2 3]) = thermaLoads(elementArray.nodes(iElement,1),[1 2 3])-Ft';
            thermaLoads(elementArray.nodes(iElement,2),[1 2 3]) = thermaLoads(elementArray.nodes(iElement,2),[1 2 3])+Ft';
            
        end
        
        pointLoadsArray = pointLoadsArray+thermaLoads;
     
        
    case 3
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
        boundaryConditionsArray([1 2 3 7 8 9],[1 2 3]) = true; % Apoyos fijos
        
          % Load definition
        pointLoadsArray = zeros(nNodes,6);     % Point load nodal value for each direction
        pointLoadsArray(5,6) = -20000/(8*pi/30)*1000; %Nmm
        pointLoadsArray(5,2) = -600000-350*9.81; %N
        massNode = 5;
        
        for i = 1:nElements
            fAccelerationLocal = elementArray.weight(i)*acceleration/2;
            
            fAcceleration = fAccelerationLocal;
            pointLoadsArray(elementArray.nodes(i,1),[1 2 3]) = pointLoadsArray(elementArray.nodes(i,1),[1 2 3]) + fAcceleration; %N
            pointLoadsArray(elementArray.nodes(i,2),[1 2 3]) = pointLoadsArray(elementArray.nodes(i,2),[1 2 3]) + fAcceleration; %N
        end
        
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
        massNode = 4;
        
        for i = 1:nElements
            fAccelerationLocal = elementArray.weight(i)*acceleration/2;
            
            fAcceleration = fAccelerationLocal;
            pointLoadsArray(elementArray.nodes(i,1),[1 2 3]) = pointLoadsArray(elementArray.nodes(i,1),[1 2 3]) + fAcceleration; %N
            pointLoadsArray(elementArray.nodes(i,2),[1 2 3]) = pointLoadsArray(elementArray.nodes(i,2),[1 2 3]) + fAcceleration; %N
        end
        
        thermaLoads = zeros(nNodes,6); %thermal loads 
        thermalElements = [3 7]; %Vector with the elements with thermal load    
        
        
        for iElement = thermalElements
            ftLocal = alpha*membersMaterial(1)*membersCrossSection(elementArray.crossSection(iElement),1)*dT; %N
            node1 = nodesPositionArray(elementArray.nodes(iElement,1),:);
            node2 = nodesPositionArray(elementArray.nodes(iElement,2),:);
            anode = structuralJointsArray(elementArray.auxiliarPoint(iElement),:);    
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
    boundaryConditionsArray(:,[3 5]) = true; 
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
    bendingStress = 0;
    bendingStressy = 0;
    bendingStressz = 0;
    
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
       tau = 0; %Bending ya lo defini como cero
   else
        A = membersCrossSection(elementArray.crossSection(iElement),1);
        Izz= membersCrossSection(elementArray.crossSection(iElement),2);
        Iyy = membersCrossSection(elementArray.crossSection(iElement),3);
        Tk= membersCrossSection(elementArray.crossSection(iElement),4);
        
        positions = [0 L];
        
        for i = 1:2 
                x = positions(i);
                Bbending = [-6/L^2+12*x/L^3, -4/L+6*x/L^2, 6/L^2-12*x/L^3, -2/L+6*x/L^2];
                
                if structuralMembersArray.circular(iElement) == 1 
                    r = sqrt(A/pi);
                    if elementArray.crossSection(iElement) == 4
                        r = 130;
                    end
                    
                    Mz = membersMaterial(1)*Izz*Bbending*[localDisplacements(1,2), localDisplacements(1,6), localDisplacements(2,2), localDisplacements(2,6)]';
                    My = membersMaterial(1)*Iyy*Bbending*[localDisplacements(1,3), -localDisplacements(1,5), localDisplacements(2,3), -localDisplacements(2,5)]';

                    M = sqrt(My^2+Mz^2); %SOLO PARA CIRCULAR 
                    bendingStress(i) = M*r/Izz;

                    %Torsion
                    %tau_rectangular = 3*T/(8*a*b^2)*(1+0.6095*(b/a)+0.8865*(b/a)^2-1.8023*(b/a)^3+0.91*(b/a)^4);
                    T = membersMaterial(2)*Tk*(localDisplacements(2,4)-localDisplacements(1,4))/L;
                    ct = r;
                    tau(i) = T*ct/Tk; %Maximo en el radio
                    
                elseif structuralMembersArray.circular(iElement) == 2 %cuadrada 120 lado
                    
                    y = 120/2; z = 120/2;                 
                    
                    % y axis bending                         
                    bendingStrainy = -y.*Bbending*[localDisplacements(1,2), localDisplacements(1,6), localDisplacements(2,2), localDisplacements(2,6)]';
                    bendingStressy(i) = membersMaterial(1)*bendingStrainy;            

                    % z axis bending
                    bendingStrainz = z.*Bbending*[localDisplacements(1,3), -localDisplacements(1,5), localDisplacements(2,3), -localDisplacements(2,5)]';
                    bendingStressz(i) = membersMaterial(1)*bendingStrainz;

                    %torsion
                    T = membersMaterial(2)*Tk*(localDisplacements(2,4)-localDisplacements(1,4))/L;
                    a = 120/2; b = 120/2;
                    tau(i) = 3*T/(8*a*b^2)*(1+0.6095*(b/a)+0.8865*(b/a)^2-1.8023*(b/a)^3+0.91*(b/a)^4);
                end
        end
   end
    [val, pos] = max(abs(bendingStressy));
    bendingStressy = val;
    [val, pos] = max(abs(bendingStressz));
    bendingStressz = val;
    [val, pos] = max(abs(bendingStress));
    bendingStress = val;
    [val, pos] = max(abs(tau));
    tau = val;
    %bendingStress se usa solo en circular, caso contrario vale cero. Los
    %otros dos son para la seccion cuadrada. 
    elementArray.stress(iElement,1) = axialStress+(bendingStress+bendingStressy+bendingStressz)*sign(axialStress); %max stress (absolute)
    elementArray.stress(iElement,2) = tau;
    elementArray.stress(iElement,3) = sqrt(elementArray.stress(iElement,1)^2+3*elementArray.stress(iElement,2)^2); %Von Mises
    elementArray.stress(iElement,4) = (bendingStressy+bendingStressz+bendingStress)*sign(axialStress); %bending Stress
        
end

%total structure weight

totalWeigth = sum(elementArray.weight)/1000; %tons

%print information
fprintf('only bars = %s and plane structure = %s .\n',string(onlyBars),string(planeStructure));

staticSF = Sy*elementArray.stress(:,3).^-1;
fprintf('ElementNumber | AxialStress | BendingStress | shear stress | Von Mises stress | SF for element:\n');
format bank
disp([[1:nElements]',elementArray.stress(:,1),elementArray.stress(:,4),elementArray.stress(:,2) ,elementArray.stress(:,3), staticSF])
fprintf('Nodal Displacements (mm/rad):\n');
format long
disp(nodalDisplacements);
fprintf('Total Structure Weight(tons): %d \n', totalWeigth);
x = [6 6 11 11]*1000;
y = [0 10 10 15]*1000;
z = [0 0 0 0]*1000;
line(x,y,z,'Color','#4DBEEE','LineWidth',1)
hold off
magnificationScale = 100;


[naturalFrequencies] = Vibrations(onlyBars,structuralJointsArray,structuralMembersArray,planeStructure,membersCrossSection,membersMaterial,boundaryConditionsArray,massNode,magnificationScale);
[filteredBucklingLoadFactor] = Buckling(structuralJointsArray,structuralMembersArray,planeStructure,membersCrossSection,membersMaterial,boundaryConditionsArray,pointLoadsArray,magnificationScale);





