files = dir('IM/*mat');
files={files.name}

 load('Mesh_16.mat');
 
  cnts = (Mesh.Nodes(Mesh.Tetra(:,1),:) +Mesh.Nodes(Mesh.Tetra(:,2),:) + Mesh.Nodes(Mesh.Tetra(:,3),:) + Mesh.Nodes(Mesh.Tetra(:,4),:) )/4;

    sel=cnts(:,3)<0.05 & cnts(:,3)>-0.05;

    Mesh1.Nodes=Mesh.Nodes;
    Mesh1.Tetra=Mesh.Tetra(sel,:);
     [Mesh1.Nodes, Mesh1.Tetra] = removeisolatednode(Mesh1.Nodes, Mesh1.Tetra);
    
for f=1:length(files)
    load(['IM/' files{f}]);
    
     for i= 1:size(X,2)
        writeVTKcell(['VTK/Raw' files{f}(1:end-4)  num2str(i)],Mesh1.Tetra,Mesh1.Nodes,X(sel,i));
    end
end