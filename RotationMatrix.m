function [lambdaMatrix] = RotationMatrix(node1, node2,anode)
%RotationMatrix te devuelve la matriz de rotacion lambda

V1 = node2-node1;
V1 =  V1/norm(V1);
V2 = anode-node1;
V2 = V2/norm(V2);
    V3 = cross(V1,V2);
    V3 = V3/norm(V3);
    V2 = cross(V3,V1);
    V2 = V2/norm(V2);
    lambdaMatrix = [V1
                    V2
                    V3];

end

