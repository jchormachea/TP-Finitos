function [naturalFrequencies] = Vibrations(onlyBars,structuralJointsArray,structuralMembersArray,planeStructure,membersCrossSection,membersMaterial,boundaryConditionsArray,puntualMassNode,magnificationScale)
%Funcion que hace las vibraciones
hold on
membersMaterial(3)=membersMaterial(3)/(1000^4);                          
puntualMass = puntualMassNode(2)/1000;%Kg/1000
massNode = puntualMassNode(1);

% Connected Dof                          
% structuralMembersArray.dof=true(size(structuralMembersArray.nodes,1),12);

% Number of elements in member
% structuralMembersArray.refinement=ones(size(structuralMembersArray.nodes,1));

% Member cross section number
% structuralMembersArray.crossSection=ones(size(structuralMembersArray.nodes,1));

% Member material number
% structuralMembersArray.material=ones(size(structuralMembersArray.nodes,1));
                    

% Cross sections definition
% Area | Inertia Moment in P123 plane | Inertia Moment orthogonal to P123 plane | Torsional Stiffness
% membersCrossSection=[A1  Izz1 Iyy1 Tk1]; % mm2 mm4

% Material definition
% Young Modulus | Transverse Modulus | Density 
% % membersMaterial=[68980 26531 2.7e-9]; %MPa kg/mm3/1000 Aluminium 6061
% membersMaterial=[200000 80000 7800/1000^3/1000]; %Material del TP

% Structure plot
linearMeshPlot(structuralMembersArray.nodes(:,1:2),structuralJointsArray,'b','Yes');
             
%% Preprocess

% Mesh generation
[elementArray,nodesPositionArray]=trussFrameMeshGenerator(structuralMembersArray,structuralJointsArray);

% Problem parameters
nElements=size(elementArray.nodes,1);    %Number of elements
nNodes=size(nodesPositionArray,1);       %Number of nodes
nTotalDof=max(max(elementArray.dof));    %Number of total dofs
        

%% Resonant modes determination

% Stiffness calculation and assembly
[stiffnessMatrix]=assemble1DStiffnessMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial);

% Stiffness calculation and assembly
[massMatrix]=assemble1DMassMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial);


sz = size(massMatrix);
puntualMassMatrix(sz(1),sz(2)) = 0;
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
% magnificationScale=100;
figure

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
hold off
end

