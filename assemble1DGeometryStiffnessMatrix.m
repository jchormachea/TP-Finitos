function [geometryStiffnessMatrix]=assemble1DGeometryStiffnessMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersAxialLoad)
% 1D finite element matrix assembler
% 
% assemble1DGeometricalStiffnessMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersAxialLoad)
%
% geometryStiffnessMatrix:      Assembled stiffness matrix
%
% elementArray.dof
% elementArray.nodes
% nodesPositionArray
% membersAxialLoad
%

%% Definitions
nElements=size(elementArray.nodes,1);      %Number of elements
nTotalDof=max(max(elementArray.dof));    %Number of Dofs

%% Stiffness matrix assembly
geometryStiffnessMatrix = zeros(nTotalDof);

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
    Y1 = 1.2/elementLength;
    Y2 = 0.1;
    Y3 = 2/15*elementLength;
    Y4 = 1/30*elementLength;
    Z1 = Y1;
    Z2 = Y2;
    Z3 = Y3;
    Z4 = Y4;
    
    elementalGeometryStiffnessMatrix=[  0   0   0  0   0   0  0   0   0  0   0   0
                                        0  Y1   0  0   0  Y2  0 -Y1   0  0   0  Y2
                                        0   0  Z1  0 -Z2   0  0   0 -Z1  0 -Z2   0
                                        0   0   0  0   0   0  0   0   0  0   0   0
                                        0   0 -Z2  0  Z3   0  0   0  Z2  0 -Z4   0
                                        0  Y2   0  0   0  Y3  0 -Y2   0  0   0 -Y4
                                        0   0   0  0   0   0  0   0   0  0   0   0
                                        0 -Y1   0  0   0 -Y2  0  Y1   0  0   0 -Y2
                                        0   0 -Z1  0  Z2   0  0   0  Z1  0  Z2   0
                                        0   0   0  0   0   0  0   0   0  0   0   0
                                        0   0 -Z2  0 -Z4   0  0   0  Z2  0  Z3   0
                                        0  Y2   0  0   0 -Y4  0 -Y2   0  0   0  Y3];
    
    elementalGeometryStiffnessMatrix = membersAxialLoad(iElement)*projectionMatrix'*elementalGeometryStiffnessMatrix*projectionMatrix;
     
    % Matrix assembly
    
    geometryStiffnessMatrix(elementArray.dof(iElement,:),elementArray.dof(iElement,:)) = geometryStiffnessMatrix(elementArray.dof(iElement,:),elementArray.dof(iElement,:)) + elementalGeometryStiffnessMatrix;

end

