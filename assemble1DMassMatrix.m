function [massMatrix]=assemble1DMassMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial)
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
massMatrix = zeros(nTotalDof);

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
    X1  = membersCrossSection(elementArray.crossSection(iElement),1)/3;
    X2  = membersCrossSection(elementArray.crossSection(iElement),1)/6;
    S1  = (membersCrossSection(elementArray.crossSection(iElement),2)+membersCrossSection(elementArray.crossSection(iElement),3))/3;
    S2  = (membersCrossSection(elementArray.crossSection(iElement),2)+membersCrossSection(elementArray.crossSection(iElement),3))/6;
    Y1 = 156/420*membersCrossSection(elementArray.crossSection(iElement),1);
    Y2 =  22/420*membersCrossSection(elementArray.crossSection(iElement),1)*elementLength;
    Y3 =  54/420*membersCrossSection(elementArray.crossSection(iElement),1);
    Y4 =  13/420*membersCrossSection(elementArray.crossSection(iElement),1)*elementLength;
    Y5 =   4/420*membersCrossSection(elementArray.crossSection(iElement),1)*elementLength^2;
    Y6 =   3/420*membersCrossSection(elementArray.crossSection(iElement),1)*elementLength^2;
    Z1 = Y1;
    Z2 = Y2;
    Z3 = Y3;
    Z4 = Y4;
    Z5 = Y5;
    Z6 = Y6;
    
    elementalMassMatrix=[ X1   0   0  0   0   0 X2   0   0  0   0   0
                           0  Y1   0  0   0  Y2  0  Y3   0  0   0 -Y4
                           0   0  Z1  0 -Z2   0  0   0  Z3  0  Z4   0
                           0   0   0 S1   0   0  0   0   0 S2   0   0
                           0   0 -Z2  0  Z5   0  0   0 -Z4  0 -Z6   0
                           0  Y2   0  0   0  Y5  0  Y4   0  0   0 -Y6
                          X2   0   0  0   0   0 X1   0   0  0   0   0
                           0  Y3   0  0   0  Y4  0  Y1   0  0   0 -Y2
                           0   0  Z3  0 -Z4   0  0   0  Z1  0  Z2   0
                           0   0   0 S2   0   0  0   0   0 S1   0   0
                           0   0  Z4  0 -Z6   0  0   0  Z2  0  Z5   0
                           0 -Y4   0  0   0 -Y6  0 -Y2   0  0   0  Y5];
    
    elementalStiffnessMatrix = membersMaterial(elementArray.material(iElement),3)*elementLength*projectionMatrix'*elementalMassMatrix*projectionMatrix;
     
    % Matrix assembly
    
    massMatrix(elementArray.dof(iElement,:),elementArray.dof(iElement,:)) = massMatrix(elementArray.dof(iElement,:),elementArray.dof(iElement,:)) + elementalStiffnessMatrix;

end

