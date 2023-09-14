function [NewJoints,NewStructural,NewThermal] = MeshRefinement(structuralJointsArray,structuralMembersArray,AuxiliarNodes, thermal)

Lold = length(structuralJointsArray);
nElements=size(structuralMembersArray.nodes,1);    %Number of elements
NewThermal = thermal;

for iElement = 1:nElements
        node1 = structuralMembersArray.nodes(iElement,1);
        node2 = structuralMembersArray.nodes(iElement,2);
        anode = structuralMembersArray.nodes(iElement,3);
        n1 = structuralJointsArray(node1,:);
        n2 = structuralJointsArray(node2,:);
    
    if structuralMembersArray.refinement(iElement) ==2
        naux = (n1+n2)/2;        
        structuralJointsArray = [structuralJointsArray(1:end-AuxiliarNodes,:);naux; structuralJointsArray(end-AuxiliarNodes+1:end,:)];
        L = length(structuralJointsArray);        
        if anode > Lold-AuxiliarNodes
            pos = Lold-anode;
            anode = L-pos;
        end
        structuralMembersArray.nodes(iElement,:) = [node1,L-AuxiliarNodes, anode];
        structuralMembersArray.nodes(end+1,:) = [L-AuxiliarNodes,node2, anode];
        last = size(structuralMembersArray.nodes,1);
        if sum(iElement==thermal) > 0
            NewThermal = [NewThermal, last]; 
        end
        
        structuralMembersArray.crossSection(last,:) = structuralMembersArray.crossSection(iElement,:);
        structuralMembersArray.circular(last) = structuralMembersArray.circular(iElement,1);
        
    elseif structuralMembersArray.refinement(iElement) ==4
        
        naux2 = (n1+n2)/2;
        naux1 = (n1+naux2)/2;
        naux3 = (naux2+n2)/2;
        structuralJointsArray = [structuralJointsArray(1:end-AuxiliarNodes,:);[naux1;naux2;naux3]; structuralJointsArray(end-AuxiliarNodes+1:end,:)];
        L = length(structuralJointsArray);
        
        if anode > Lold-AuxiliarNodes
            pos = Lold-anode;
            anode = L-pos;
        end
        
        structuralMembersArray.nodes(iElement,:) = [node1,L-AuxiliarNodes-2, anode];
        structuralMembersArray.nodes(end+1,:) = [L-AuxiliarNodes-2,L-AuxiliarNodes-1, anode];
        structuralMembersArray.nodes(end+1,:) = [L-AuxiliarNodes-1,L-AuxiliarNodes, anode];
        structuralMembersArray.nodes(end+1,:) = [L-AuxiliarNodes,node2, anode];
        last = size(structuralMembersArray.nodes,1);
        if sum(iElement==thermal) > 0
            NewThermal = [NewThermal, last-2:last]; 
        end
        structuralMembersArray.crossSection(last-2,:) = structuralMembersArray.crossSection(iElement,:);
        structuralMembersArray.crossSection(last-1,:) = structuralMembersArray.crossSection(iElement,:);
        structuralMembersArray.crossSection(last,:) = structuralMembersArray.crossSection(iElement,:);
        structuralMembersArray.circular(last-2:last,:) = structuralMembersArray.circular(iElement);
        
    elseif structuralMembersArray.refinement(iElement) ==8        
        
        naux4 = (n1+n2)/2;
        naux2 = (n1+naux4)/2;
        naux6 = (naux4+n2)/2;
        naux3 = (naux2+naux4)/2;
        naux1 = (n1+naux2)/2;
        naux5 = (naux4+naux6)/2;
        naux7 = (naux6+n2)/2;
        structuralJointsArray = [structuralJointsArray(1:end-AuxiliarNodes,:);[naux1;naux2;naux3;naux4;naux5;naux6;naux7]; structuralJointsArray(end-AuxiliarNodes+1:end,:)];
        L = length(structuralJointsArray);
        
        if anode > Lold-AuxiliarNodes
            pos = Lold-anode;
            anode = L-pos;
        end
        
        structuralMembersArray.nodes(iElement,:) = [node1,L-AuxiliarNodes-6, anode];
        structuralMembersArray.nodes(end+1,:) = [L-AuxiliarNodes-6,L-AuxiliarNodes-5, anode];
        structuralMembersArray.nodes(end+1,:) = [L-AuxiliarNodes-5,L-AuxiliarNodes-4, anode];
        structuralMembersArray.nodes(end+1,:) = [L-AuxiliarNodes-4,L-AuxiliarNodes-3, anode];
        structuralMembersArray.nodes(end+1,:) = [L-AuxiliarNodes-3,L-AuxiliarNodes-2, anode];
        structuralMembersArray.nodes(end+1,:) = [L-AuxiliarNodes-2,L-AuxiliarNodes-1, anode];
        structuralMembersArray.nodes(end+1,:) = [L-AuxiliarNodes-1,L-AuxiliarNodes, anode];
        structuralMembersArray.nodes(end+1,:) = [L-AuxiliarNodes,node2, anode];
        last = size(structuralMembersArray.nodes,1);
        
        if sum(iElement==thermal) > 0
            NewThermal = [NewThermal, last-6:last]; 
        end
        for i = 0:6
            structuralMembersArray.crossSection(last-i,:) = structuralMembersArray.crossSection(iElement,:);
            structuralMembersArray.circular(last-i,:) = structuralMembersArray.circular(iElement);
        end
        
    end
    
end

NewJoints = structuralJointsArray;
NewStructural = structuralMembersArray;

end

