%% Constants
TI = 1;
GL = 0;

D = containers.Map({TI,GL},{eye(2)*17, eye(2)*0.8}); %W/(m K)
alpha_c = 100; %W/(m^2K)
T_0 = 20; %Celcius, T_0 = zero stress from heat.

rho = containers.Map({TI,GL},{4620,3860}); %Density, kg/m^3
E = containers.Map({TI,GL},{110, 67}); %Young's modulus, GPa
alpha = containers.Map({TI,GL},{9.4e-6,7e-6}); %Expansion coefficient, 1/K
c_p = containers.Map({TI,GL},{523, 670}); %Specific heat, J/(kg K)
Poisson = containers.Map({TI,GL},{0.34, 0.2}); % Poission's ratio [-]

thickness = 0.01; % [m]

%% Extract calfem notation from pdetool-mesh

%element data
enod=t(1:3,:)'; % nodes of elements
nelm = size(enod,1); % number of elements
emat = (t(4,:)==TI)';  %  Material type of element

%node data
coord = p'/100;
nnod=size(coord,1); % number of nodes
dof=(1:nnod)'; % dof (degrees of freedom) number is node number
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; % give each dof a number
edof_S = zeros(nelm, 7);   
edof = zeros(nelm,4);   %Element degrees of freedom

for ie=1:nelm
    edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
    edof(ie,:)=[ie,enod(ie,:)];
end

% separate edges with convection and fixed edges
er = e([1 2 5],:); % Reduced e [2x node number; edge label]
conv_segments = [2 15 14 23 31, 1]; % Segments with convection.
edges_conv = [];

fixed_segments = [21 22 1 16 17 18 19 20];
edges_fixed = [];

for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv [er(1:2,i);er(3,i)==1]]; % [node labels of convection elements;convection type]
            %convection type = 0 -> T_inf outside
            %convection type = 1 -> T_c outside
    end
    if ismember(er(3,i),fixed_segments)
        edges_fixed = [edges_fixed [er(1:2,i);er(3,i)==1]]; % [node labels of convection elements; fixed type]
            %fixed type = 0 -> u_y = 0
            %fixed type = 1 -> u_x = 0
    end
end


% find edges with given heat flow

%% a) Stationary temperature distribution
%Define constants
T_outside = [-96 20]; % [T_inf T_c]

%allocate memory
K = zeros(nnod);
F = zeros(nnod,1);

%Basic stiffness matrix
for ie = 1:nelm
    %find coords of element
    ex = coord(enod(ie,:),1)';
    ey = coord(enod(ie,:),2)';

    Ke = flw2te(ex,ey,thickness,D(emat(ie)));
    
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
a0 = solveq(K,F);
[ex, ey] = coordxtr(edof, coord, dof, 3);
ed = extract(edof, a0);

%Plot
clf
patch(ex',ey',ed','EdgeColor','none');
hold on
patch(ex',-ey',ed','EdgeColor','none');

caxis([-100 50]);
axis([-0.1 1.2 -0.5 0.5]/100);
title("Temperature distribution, T_{\infty} = " + (T_outside(1)) + " [C]");
%colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-postition [m]');

disp("MAX TEMP: " + max(max(ed)));
%Max T vid T_inf = -96: 19,9769
%Max T vid T_inf = 40: 40,0093

%% b) transient heat
%Run a) to generate a0 with desired initial condition.
delta_t = 100;
T_outside = [-96 20]; %[T_inf T_c]
a = a0;

A = zeros(nnod);
F = zeros(nnod,1);

%Generate F with new boundary cond.
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
for ie = 1:nelm
    %find coords of element
    ex = coord(enod(ie,:),1);
    ey = coord(enod(ie,:),2);
   
    material = emat(ie);
    Ae = plantml(ex',ey',rho(material)*c_p(material));
    
    indx = edof(ie,2:end);  % where to insert Ae
    A(indx, indx) = A(indx,indx)+Ae;
end

% Implicit Euler time step
time_step = @(a) (A + delta_t*K)\(F*delta_t+A*a);

% Step
%Suggested steps: t= 200 1200 10,000 20,000
for i=1:200
    a = time_step(a);
end

[ex, ey] = coordxtr(edof, coord, dof, 3);

ed = extract(edof, a);

%Plot
patch(ex',ey',ed','EdgeColor','none');
hold on
patch(ex',-ey',ed','EdgeColor','none');

caxis([-100 50]);
axis([-0.1 1.2 -0.5 0.5]/100);
title("t = " + (i*delta_t) + " s");
colormap default;
colorbar;
xlabel('x-position [m]');
ylabel('y-postition [m]');

disp("MAX TEMP: " + max(max(ed)));
disp("MIN TEMP: " + min(min(ed)));


%% C Mises stress field
% Plain strain conditions hold -> epsilon = [exx eyy exy]'
% sigma_zz != 0 
% Ka = f_b + f_l + f_0
% K, Stiffness matrix, = int BtDtB dA = D*t*int BtB dA.
% f_0, initial vector from temperature distribution = t* int_ Bt epsilon_0 dA
% f_b from know traction at boundaries = 0.
% f_l = 0 because of no internal load. 

% Define isotropic D matrix for plane strain (13.35 in book)
calcD = @(E, v) E/(1+v)*(1-2*v)*[(1-v) v 0; v (1-v) 0; 0 0 (1-2*v)/2];
D_el = containers.Map({TI,GL},{calcD(E(TI),Poisson(TI)),calcD(E(TI),Poisson(GL))});

% Define constant part of D*epsilon_0 (13.34 in book)
calcConstEps_0 = @(mat) alpha(mat)*E(mat)/(1-2*Poisson(mat))*[1;1;0];
const_eps0 = containers.Map({TI,GL},{calcConstEps_0(TI),calcConstEps_0(GL)});

%Calculate K and f_0
K = zeros(nnod*2);
F = zeros(nnod*2,1);

for ie = 1:nelm
    ex = coord(enod(ie,:),1)';
    ey = coord(enod(ie,:),2)';
    material = emat(ie);
    
    Ke = plante(ex, ey, [2 thickness], D_el(material));
    
    %Absolut inte hundra på detta.
    dT = mean(ed(ie,:))-T_0; %Hämta delta_T jfrt stressfri i elementet beräknad i a) (medelvärdet av nodernas temperatur)
    es = const_eps0(material)*dT;
    f_0e = plantf(ex, ey, [2 thickness], es'); %plantf beräknar int Bt es t dA.
    
    %Insert
    indx = edof_S(ie,2:end);  % where to insert
    K(indx, indx) = K(indx,indx)+Ke;
    F(indx) = F(indx) + f_0e;
end

%Calculate bc, known extrusions (=0). 
%Koden blir trixig här då jag hämtat ut alla edges som har u = 0 vilket gör
%att det kan förekomma dubletter av noder.

already_added = zeros(nnod*2,1);
bc = [];
for ib = 1:length(edges_fixed)
    edge = edges_fixed(:,ib);
    x_or_y = 2-edge(3);
    
    node_id = dof_S(edge(1),x_or_y);
    if(already_added(node_id)==0)
        bc = [bc; node_id, 0];
        already_added(node_id) = 1;
    end
    
    node_id = dof_S(edge(2),x_or_y);
    if(already_added(node_id)==0)
        bc = [bc; node_id, 0];
        already_added(node_id) = 1;
    end
end

%solve system
a_S = solveq(K,F, bc);
[ex_S, ey_S] = coordxtr(edof_S, coord, dof_S, 3);
ed_S = extract(edof_S, a_S);

%Calculate von Mises stress per element
Seff_el = zeros(nelm,1);
for ie = 1: nelm
    ex = coord(enod(ie,:),1)';
    ey = coord(enod(ie,:),2)';
    material = emat(ie);
    
    a_index = [dof_S(enod(ie,:), 1) ; dof_S(enod(ie,:), 2)];
    
    % Boken har något om temperaturen här också men är osäker på om det
    % behövs/ hur det skulle fungera?
    sigma = plants(ex, ey, [2 thickness], D_el(material),a_S(a_index)'); % [sigma_xx sigma_yy sigma_xy]
    sigma_zz = Poisson(material)*(sigma(1) + sigma(2)); % 13.42 i boken
    
    vonMisesSquared = sigma*sigma' + sigma_zz^2 - sigma(1)*sigma(2)-sigma(1)*sigma_zz-sigma(2)*sigma_zz+2*sigma(3)^2;
    Seff_el(ie) = sqrt(vonMisesSquared);
end

%Calculate von Mieses stress per node (mean of connected elements)
Seff_nod = zeros(nnod,1);
for i=1:nnod
    [c0,c1]=find(edof(:,2:4)==i);
    Seff_nod(i,1)=sum(Seff_el(c0))/size(c0,1);
end
%Extract stress per node to patch format
eM = extract(edof, Seff_nod);

% Plot stress Field
patch(ex_S',ey_S',eM','EdgeColor','none');
hold on
patch(ex_S',-ey_S',eM','EdgeColor','none');

axis([-0.1 1.2 -0.5 0.5]/100);
title("Von Mises effective stress field.");
colormap default;
colorbar;
xlabel('x-position [m]');
ylabel('y-postition [m]');


%% plot displacements
mag = 5;
exd = ex_S + mag*ed_S(:,1:2:end);
eyd = ey_S + mag*ed_S(:,2:2:end);

figure();
patch(ex_S',ey_S',[0 0 0],"EdgeColor","none","FaceAlpha",0.3);
hold on
patch(ex_S',-ey_S',[0 0 0],"EdgeColor","none","FaceAlpha",0.3);
patch(exd',eyd',[0 0 0],"FaceAlpha",0.3);
patch(exd',-eyd',[0 0 0],"FaceAlpha",0.3);
axis equal



