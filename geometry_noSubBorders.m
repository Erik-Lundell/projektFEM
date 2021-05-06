% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
% 1) Export the required variables from pdetool and create a MATLAB script
%    to perform operations on these.
% 2) Define the problem completely using a MATLAB script. See
%    https://www.mathworks.com/help/pde/examples.html for examples
%    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[1 0.49828106450336368 1]);
set(ax,'PlotBoxAspectRatio',[1.743983725761773 1 2.4914053225168185]);
set(ax,'XLimMode','auto');
set(ax,'YLimMode','auto');
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');
pdetool('gridon','on');

% Geometry description:
pderect([0 1.1000000000000001 0 0.25],'R1');
pderect([0 1 0.14999999999999999 0],'R2');
pderect([0 0.10000000000000001 0.25 0.29999999999999999],'R3');
pderect([0.90000000000000002 1 0.050000000000000003 0],'R4');
pdeellip(0.10000000000000001,0,0.10000000000000001,0.25,...
0,'E1');
pdeellip(0.40000000000000002,0,0.10000000000000001,0.25,...
0,'E2');
pderect([0.59999999999999998 0.69999999999999996 0.14999999999999999 0],'R5');
pderect([0.10000000000000001 0.20000000000000001 0.14999999999999999 0],'R6');
pderect([0.40000000000000002 0.5 0.14999999999999999 0],'R7');
pderect([0 0.20000000000000001 0 -0.25],'R8');
pderect([0.29999999999999999 0.5 0 -0.25],'R9');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1-R2+R3+R4+E1-R8-R6+E2-R9-R7+R5')

% Boundary conditions:
pdetool('changemode',0)
pdetool('removeb',[8 ]);
pdetool('removeb',[33 ]);
pdetool('removeb',[32 ]);
pdetool('removeb',[34 ]);
pdetool('removeb',[33 ]);
pdetool('removeb',[8 ]);
pdesetbd(32,...
'dir',...
1,...
'1',...
'0')
pdesetbd(31,...
'dir',...
1,...
'1',...
'0')
pdesetbd(30,...
'dir',...
1,...
'1',...
'0')
pdesetbd(29,...
'dir',...
1,...
'1',...
'0')
pdesetbd(27,...
'dir',...
1,...
'1',...
'0')
pdesetbd(26,...
'dir',...
1,...
'1',...
'0')
pdesetbd(25,...
'dir',...
1,...
'1',...
'0')
pdesetbd(23,...
'dir',...
1,...
'1',...
'0')
pdesetbd(22,...
'dir',...
1,...
'1',...
'0')
pdesetbd(21,...
'dir',...
1,...
'1',...
'0')
pdesetbd(20,...
'dir',...
1,...
'1',...
'0')
pdesetbd(19,...
'dir',...
1,...
'1',...
'0')
pdesetbd(18,...
'dir',...
1,...
'1',...
'0')
pdesetbd(17,...
'dir',...
1,...
'1',...
'0')
pdesetbd(16,...
'dir',...
1,...
'1',...
'0')
pdesetbd(15,...
'dir',...
1,...
'1',...
'0')
pdesetbd(14,...
'dir',...
1,...
'1',...
'0')
pdesetbd(13,...
'dir',...
1,...
'1',...
'0')
pdesetbd(11,...
'dir',...
1,...
'1',...
'0')
pdesetbd(10,...
'dir',...
1,...
'1',...
'0')
pdesetbd(9,...
'dir',...
1,...
'1',...
'0')
pdesetbd(8,...
'dir',...
1,...
'1',...
'0')
pdesetbd(7,...
'dir',...
1,...
'1',...
'0')
pdesetbd(6,...
'dir',...
1,...
'1',...
'0')
pdesetbd(5,...
'dir',...
1,...
'1',...
'0')
pdesetbd(4,...
'dir',...
1,...
'1',...
'0')
pdesetbd(3,...
'dir',...
1,...
'1',...
'0')
pdesetbd(2,...
'dir',...
1,...
'1',...
'0')
pdesetbd(1,...
'dir',...
1,...
'1',...
'0')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'10.0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0 ';...
'0.0 ';...
'10.0';...
'1.0 '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','1000','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');
