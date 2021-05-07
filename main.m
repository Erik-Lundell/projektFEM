%% Extract calfem notation from pdetool-mesh
load('data.mat');

TI = 1;
GL = 0;

%element data
nelm = size(enod,1); % number of elements
enod=t(1:3,:)'; % nodes of elements
emat = (t(4,:)==TI)';  %  Material type of element

%node data
coord = p';
nnod=size(coord,1); % number of nodes
dof=(1:nnod)'; % dof (degrees of freedom) number is node number
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; % give each dof a number
edof_S = zeros(nelm, 7);   
edof = zeros(nelm,4);   %Element degrees of freedom

for ie=1:nelm
    edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
    edof(ie,:)=[ie,enod(ie,:)];
end

% separate edges conv and heat flow
er = e([1 2 5],:); % Reduced e [2x node number; edge label]
conv_segments = [2 15 14 23 31,16, 1]; % Segments with convection.
edges_conv = [];

for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv [er(1:2,i);er(5,i)==1]]; % [node labels of convection elements;convection type]
            %convection type = 0 -> T_inf outside
            %convection type = 1 -> T_c outside
    end
end

% find edges with given heat flow

%% a) Stationary temperature distribution
%Define constants
T_0 = 20; %base level as T_0 = zero stress.
T_outside = [40-T_0 20-T_0] % [T_inf T_c]
ep = 0.01; % m

D = containers.Map;
D(TI) = eye(2)*17; %W/(mK)
D(GL) = eye(2)*0.8; %W/(mK)

alpha_c = 100 %W/(m^2K)

%allocate memory
K = zeros(nnod,nnod);
F = zeros(nnod,1);

%Basic stiffness matrix
for ie = 1:nelm
    %find coords of element
    ex = coord(enod(ie,:),1);
    ey = coord(enod(ie,:),2);

    Ke = flw2te(ex,ey,ep,D(emat(ie)));
    
    indx = edof(ie,2:end);  % where to insert Ke
    K(indx,indx) = K(indx,indx)+Ke;  % insert
end

%Load vector from heat flow = 0, as q_e = 0 along all heat flow boundaries && Q=0
%Convection load vector + stiffness matrix
%qe = alpha(T_outside - T) ->
    %f_c = int_L N^T alpha_c*T_outside dL = alpha_c*T_outside * l/2
    %K_c = int_L N^T alpha_c N dL = alpha_c*l*[1/3 1/6; 1/6 1/3]
    % where L = boundary line segment (l = length of segment)
    % N = linear functions from

fce_const = [1; 1]*alpha_c/2;
Kce_const = [2 1; 1 2]*alpha_c/6;
    
for ib = 1:length(edges_conv)
    %find coordinates of boundary nodes
    indx = edges_conv(1:2,ib);
    ex = coord(indx,1);
    ey = coord(indx,2);
    l = sqrt((ex(1)-ex(2))^2+(ey(1)-ey(2))^2);  %calculate edge length
    
    fce = fce_const*l*T(edges_conv(3,ib));
    Kce = Kce_const*l;

    K(indx,indx) = K(indx,indx)+Kce;  % insert
    F(indx) = F(indx) + fce;
end

%Solve system, bc = [dof (=node number tror jag....) temp]
%Oklart vad vi har här men kanske T_inf vänstra boundary och T_c högra?
%eller bara högra?
bc = [???????????];
T = solveq(K,F,bc);

ex = coord(:,1);
ey = coord(:,2);
ed = extract(edof, T);
patch(ex',ey',ed');


