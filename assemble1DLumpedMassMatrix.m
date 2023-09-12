function [lumpedMassMatrix]=assemble1DLumpedMassMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial)
% 1D finite element matrix assembler
% 
% assemble1DGeometricalStiffnessMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersAxialLoad)
%
% geometryMassMatrix:      Assembled stiffness matrix
%
% elementArray.dof
% elementArray.nodes
% nodesPositionArray
% membersCrossSection
% membersMaterial
%

%% Definitions
nElements=size(elementArray.nodes,1);      %Number of elements
nTotalDof=max(max(elementArray.dof));    %Number of Dofs

%% Stiffness matrix assembly
lumpedMassMatrix = zeros(nTotalDof);

for iElement = 1:nElements
    % Rotation
    V1 = nodesPositionArray(elementArray.nodes(iElement,2),:) - nodesPositionArray(elementArray.nodes(iElement,1),:);
    elementLength = norm(V1);
    V1 = V1/elementLength;
    V2 = structuralJointsArray(elementArray.auxiliarPoint(iElement),:) - nodesPositionArray(elementArray.nodes(iElement,1),:);
    V2 = V2/norm(V2);
    V3 = cross(V1,V2);
    V3 = V3/norm(V3);
    V2 = cross(V3,V1);
    V2 = V2/norm(V2);
    lambdaProjectionMatrix = [V1
                              V2
                              V3];
    projectionMatrix = blkdiag(lambdaProjectionMatrix,lambdaProjectionMatrix,lambdaProjectionMatrix,lambdaProjectionMatrix);
    
    % Coefficients
    membersCrossSection(elementArray.crossSection(iElement),1)
    X = membersCrossSection(elementArray.crossSection(iElement),1)/2;
    Y = 12*membersCrossSection(elementArray.crossSection(iElement),1)/24;
    Z = Y;
    Xx = membersCrossSection(elementArray.crossSection(iElement),4)/2;
    Yy = membersCrossSection(elementArray.crossSection(iElement),1)*elementLength^2/24;
    Zz = Yy;
    
    
    
    elementalMassMatrix=[  X   0   0  0   0   0  0   0   0  0   0   0
                           0   Y   0  0   0   0  0   0   0  0   0   0
                           0   0   Z  0   0   0  0   0   0  0   0   0
                           0   0   0 Xx   0   0  0   0   0  0   0   0
                           0   0   0  0  Yy   0  0   0   0  0   0   0
                           0   0   0  0   0  Zz  0   0   0  0   0   0
                           0   0   0  0   0   0  X   0   0  0   0   0
                           0   0   0  0   0   0  0   Y   0  0   0   0
                           0   0   0  0   0   0  0   0   Z  0   0   0
                           0   0   0  0   0   0  0   0   0 Xx   0   0
                           0   0   0  0   0   0  0   0   0  0  Yy   0
                           0   0   0  0   0   0  0   0   0  0   0  Zz];
    
    factor = membersMaterial(elementArray.material(iElement),3)*elementLength;
    elementalMassMatrix = factor*projectionMatrix'*elementalMassMatrix*projectionMatrix;
     
    % Matrix assembly
    
    lumpedMassMatrix(elementArray.dof(iElement,:),elementArray.dof(iElement,:)) = lumpedMassMatrix(elementArray.dof(iElement,:),elementArray.dof(iElement,:)) + elementalMassMatrix;

end

