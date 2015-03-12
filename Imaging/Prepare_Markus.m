
addpath( genpath('/home/kdawg/Desktop/PEITS/dune-peits/matlab'))

mname='GSK_16_1';

load('Mesh_16.mat');
load pos.mat;
load('for_image_EIT3_P_6025Hz_10uA_20s_Prt1_200ms_..mat','Prt_0');
pos = [pos zeros(length(pos),1)];
electrodepositions = [pos; [0,0,1]] ;

z_gnd=min(Mesh.Nodes(:,3));
x_gnd=0;
y_gnd=0;
gnd_pos = [x_gnd,y_gnd,z_gnd];
sigma=0.3*ones(length(Mesh.Tetra),1);

srf = dubs3_2(Mesh.Tetra(:,1:4));
el_nds=unique([srf(:,1);srf(:,2);srf(:,3)],'rows');

% DisplayBoundaries(Mesh);
% hold on;

elecnodfile = fopen(['electrode_nodes_' mname],'w');

Prt=[Prt_0,repmat(length(electrodepositions),length(Prt_0),1)];

dlmwrite('Prt.txt',Prt);

for i = 1: length(pos)
    dist = ((Mesh.Nodes(el_nds,1)-pos(i,1)).^2 + (Mesh.Nodes(el_nds,2)-pos(i,2)).^2).^0.5;
    nodes_ = el_nds(dist<0.095 & abs(Mesh.Nodes(el_nds,3))<0.5);
    %scatter3(Mesh.Nodes(nodes_,1),Mesh.Nodes(nodes_,2),Mesh.Nodes(nodes_,3),50,'filled')
    fprintf(elecnodfile,'%d,',nodes_(1:end-1));
    fprintf(elecnodfile,'%d\n',nodes_(end));
end


nodes_ = el_nds(Mesh.Nodes(el_nds,3)>0.9 & Mesh.Nodes(el_nds,3)<0.999);
%scatter3(Mesh.Nodes(nodes_,1),Mesh.Nodes(nodes_,2),Mesh.Nodes(nodes_,3),50,'filled');
fprintf(elecnodfile,'%d,',nodes_(1:end-1));
fprintf(elecnodfile,'%d\n',nodes_(end));
fclose(elecnodfile);


dune_exporter(Mesh.Nodes(:,1:3)/1000,Mesh.Tetra(:,1:4),sigma,'',[mname '.dgf'],electrodepositions/1000,gnd_pos/1000);
