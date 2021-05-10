%% Constants
TI = 1;
GL = 0;

D = containers.Map({TI,GL},{eye(2)*17, eye(2)*0.8}); %W/(m K)
alpha_c = 100; %W/(m^2K)

rho = containers.Map({TI,GL},{4620,3860}); %Density, kg/m^3
E = containers.Map({TI,GL},{110, 67}); %Young's modulus, GPa
alpha = containers.Map({TI,GL},{9.4e-6,7e-6}); %Expansion coefficient, 1/K
c_p = containers.Map({TI,GL},{523, 670}); %Specific heat, J/(kg K)


%% Extract calfem notation from pdetool-mesh
load('data.mat');


%element data
enod=t(1:3,:)'; % nodes of elements
nelm = size(enod,1); % number of elements
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
        edges_conv = [edges_conv [er(1:2,i);er(3,i)==1]]; % [node labels of convection elements;convection type]
            %convection type = 0 -> T_inf outside
            %convection type = 1 -> T_c outside
    end
end

% find edges with given heat flow

%% a) Stationary temperature distribution
%Define constants
T_0 = 20; %base level as T_0 = zero stress.
T_outside = [40-T_0 20-T_0]; % [T_inf T_c]
ep = 0.01; % m

%allocate memory
K = zeros(nnod);
F = zeros(nnod,1);

%Basic stiffness matrix
for ie = 1:nelm
    %find coords of element
    ex = coord(enod(ie,:),1)';
    ey = coord(enod(ie,:),2)';

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
    
    convectionTemp = edges_conv(3,ib)+1;
    fce = fce_const*l*T_outside(convectionTemp);
    Kce = Kce_const*l;

    K(indx,indx) = K(indx,indx)+Kce;  % insert
    F(indx) = F(indx) + fce;
end

%Solve system
a0 = solveq(K,F)+T_0;
[ex, ey] = coordxtr(edof, coord, dof, 3);
ed = extract(edof, a0);

%Plot
patch(ex',ey',ed','EdgeColor','none');
hold on
patch(ex',-ey',ed','EdgeColor','none');

caxis([-100 50]);
axis([-0.1 1.2 -0.5 0.5]);
title("Temperature distribution, T_{\infty} = " + (T_outside(1)+T_0) + " [C]");
%colormap(cold);
colorbar;
xlabel('x-position [m]')
ylabel('y-postition [m]');

disp("MAX TEMP: " + max(max(ed)));
%Max T vid T_inf = -96: 19,9769
%Max T vid T_inf = 40: 40,0093

%% b) transient heat
%Run a) to generate a0 with desired initial condition.
delta_T = 1;
T_outside = [-96-T_0 20-T_0];
a = a0;

%Generate F with new boundary cond.
A = zeros(nnod);
F = zeros(nnod,1);

fce_const = [1; 1]*alpha_c/2;
for ib = 1:length(edges_conv)
    %find coordinates of boundary nodes
    indx = edges_conv(1:2,ib);
    ex = coord(indx,1);
    ey = coord(indx,2);
    l = sqrt((ex(1)-ex(2))^2+(ey(1)-ey(2))^2);  %calculate edge length
    
    convectionTemp = edges_conv(3,ib)+1;
    fce = fce_const*l*T_outside(convectionTemp);
    F(indx) = F(indx) + fce;
end

%Generate derivate matrix A = rho c int N^T N dV 
%MEN HUR?????????????????????

% Implicit Euler time step
time_step = @(a) (A + delta_T*K)\(F*delta_t+A*a);

% Step
a1 = time_step(a0);

[ex, ey] = coordxtr(edof, coord, dof, 3);
ed = extract(edof, a1);

%Plot
patch(ex',ey',ed','EdgeColor','none');
hold on
patch(ex',-ey',ed','EdgeColor','none');

caxis([0 50]);
axis([-0.1 1.2 -0.5 0.5]);
title("Temperature distribution, T_{\infty} = " + (T_outside(0)+T_0) + " [C]");
%colormap(cold);
colorbar;
xlabel('x-position [m]')
ylabel('y-postition [m]');

