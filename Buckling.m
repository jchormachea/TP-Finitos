function  Buckling(structuralJointsArray,structuralMembersArray,planeStructure,membersCrossSection,membersMaterial,boundaryConditionsArray,pointLoadsArray,magnificationScale)
%Funcion que hace el buckling
hold on
% Connected Dof                          
 structuralMembersArray.dof=true(size(structuralMembersArray.nodes,1),12);

% Number of elements in member
structuralMembersArray.refinement=ones(size(structuralMembersArray.nodes,1));

% Member cross section number
structuralMembersArray.crossSection=ones(size(structuralMembersArray.nodes,1));

% Member material number
structuralMembersArray.material=ones(size(structuralMembersArray.nodes,1));

                    
% Cross sections definition
% Area | Inertia Moment in P123 plane | Inertia Moment orthogonal to P123 plane | Torsional Stiffness
% membersCrossSection=[6980 16011520 108608400 278528]; % I Beam Section
% membersCrossSection=[A1  Izz1 Iyy1 Tk1]; % mm2 mm4

% Material definition
% Young Modulus | Transverse Modulus | Density 
% membersMaterial=[68980 26531 2700]; %MPa kg/m3 Aluminium 6061
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

%% Static Stress Solver

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
% magnificationScale=100;
figure
% Nodal displacements rearranged
nodalDisplacements=reshape(displacementsVector,6,size(nodesPositionArray,1))';
nodalPositions=nodesPositionArray+nodalDisplacements(:,1:3)*magnificationScale;

% Deformed Structure plot
linearMeshPlot(elementArray.nodes(:,1:2),nodalPositions,'r','No');

% Element resultant loads
[elementLocalNodalLoads]=getElementLocalNodalLoads(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial,displacementsVector);

%% Buckling Load Factor Solver

% Stiffness calculation and assembly
[geometricalStiffnessMatrix]=assemble1DGeometryStiffnessMatrix(elementArray,nodesPositionArray,structuralJointsArray,elementLocalNodalLoads(:,1));

% Eigenvalues determination
[bucklingEigenmode bucklingLoadFactor] = eig(stiffnessMatrix(isFree,isFree),geometricalStiffnessMatrix(isFree,isFree));

% Infinite eigenvalue elimination
activeBucklingModes=find(~isinf(diag(bucklingLoadFactor)));

% Ordering and filtering
[filteredBucklingLoadFactor blucklingLoadFactorOrder]=sort(diag(bucklingLoadFactor(activeBucklingModes,activeBucklingModes)));
displacementsBucklingModes = zeros(nTotalDof,size(activeBucklingModes,1));
displacementsBucklingModes(isFree,:) = displacementsBucklingModes(isFree,:) + bucklingEigenmode(:,activeBucklingModes(blucklingLoadFactorOrder));

% First buckling eigenmode plot
nodalDisplacements=reshape(displacementsBucklingModes(:,1),6,size(nodesPositionArray,1))';
nodalPositions=nodesPositionArray+nodalDisplacements(:,1:3)*magnificationScale;
linearMeshPlot(elementArray.nodes(:,1:2),nodalPositions,'g','No');

% Second buckling eigenmode plot
nodalDisplacements=reshape(displacementsBucklingModes(:,2),6,size(nodesPositionArray,1))';
nodalPositions=nodesPositionArray+nodalDisplacements(:,1:3)*magnificationScale;
linearMeshPlot(elementArray.nodes(:,1:2),nodalPositions,'m','No');

% Third buckling eigenmode plot
nodalDisplacements=reshape(displacementsBucklingModes(:,3),6,size(nodesPositionArray,1))';
nodalPositions=nodesPositionArray+nodalDisplacements(:,1:3)*magnificationScale;
linearMeshPlot(elementArray.nodes(:,1:2),nodalPositions,'c','No');

title({['1st Buckling Load Factor   ' num2str(filteredBucklingLoadFactor(1))], ['2nd Buckling Load Factor ' num2str(filteredBucklingLoadFactor(2))], ['3rd Buckling Load Factor ' num2str(filteredBucklingLoadFactor(3))]});
view ([1 1 1]);
hold off
end

