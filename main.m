%% Extract calfem notation from pdetool-mesh
enod=t(1:3,:)'; % nodes of elements
nelm=size(enod,1); % number of elements
nnod=size(coord,1); % number of nodes
dof=(1:nnod)'; % dof number is node number
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; % give each dof a number
edof_S = zeros(nelm, 3);
edof = zeros(nelm,2);

for ie=1:nelm
    edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
    edof(ie,:)=[ie,enod(ie,:)];
end

% Check which segments that should have convections
er = e([1 2 5],:); % Reduced e
conv_segments = [2 15 14 23 31,16, 1]; % Segments with convection.
edges_conv = [];
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv er(1:2,i)];
    end

%% a) Stationary temperature distribution
T_inf = 40; 
T_c = 20;

%


