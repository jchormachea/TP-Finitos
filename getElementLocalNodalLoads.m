function [elementLocalNodalLoads]=getElementLocalNodalLoads(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial,displacementsVector)
% Determines the nodal loads of a beam element in local coordinates
% 
% getElementLocalNodalLoads(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial,displacementsVector)
%
% elementLocalNodalLoads:      Nodal loads vector in local coordinates
%                              Axial' ShearY' ShearZ' Torsion BendingY' BendingZ'
%
% elementArray.dof
% elementArray.nodes
% nodesPositionArray
% membersCrossSection
% membersMaterial
% displacementsVector
%

%% Definitions
nElements=size(elementArray.nodes,1);    %Number of elements
nTotalDof=max(max(elementArray.dof));    %Number of Dofs

%% Stiffness matrix assembly
elementLocalNodalLoads = zeros(nElements,12);

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
    X  = membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),1)/elementLength;
    S  = membersMaterial(elementArray.material(iElement),2)*membersCrossSection(elementArray.crossSection(iElement),4)/elementLength;
    Y1 = 12*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),2)/elementLength^3;
    Y2 = 6*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),2)/elementLength^2;
    Y3 = 4*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),2)/elementLength^1;
    Y4 = 2*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),2)/elementLength^1;
    Z1 = 12*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),3)/elementLength^3;
    Z2 = 6*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),3)/elementLength^2;
    Z3 = 4*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),3)/elementLength^1;
    Z4 = 2*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),3)/elementLength^1;
    
    elementalStiffnessMatrix=[  X   0   0  0   0   0 -X   0   0  0   0   0
                                0  Y1   0  0   0  Y2  0 -Y1   0  0   0  Y2
                                0   0  Z1  0 -Z2   0  0   0 -Z1  0 -Z2   0
                                0   0   0  S   0   0  0   0   0 -S   0   0
                                0   0 -Z2  0  Z3   0  0   0  Z2  0  Z4   0
                                0  Y2   0  0   0  Y3  0 -Y2   0  0   0  Y4
                               -X   0   0  0   0   0  X   0   0  0   0   0
                                0 -Y1   0  0   0 -Y2  0  Y1   0  0   0 -Y2
                                0   0 -Z1  0  Z2   0  0   0  Z1  0  Z2   0
                                0   0   0 -S   0   0  0   0   0  S   0   0
                                0   0 -Z2  0  Z4   0  0   0  Z2  0  Z3   0
                                0  Y2   0  0   0  Y4  0 -Y2   0  0   0  Y3];
    
    elementLocalNodalLoads(iElement,:) = (elementalStiffnessMatrix*projectionMatrix*displacementsVector(elementArray.dof(iElement,:)))';
 
end

