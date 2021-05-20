function Te=NtN(ex,ey,t)
% T=(ex,ey,t)
%-------------------------------------------------------------
% PURPOSE
%  Compute the quantity: Te=int(N^T*N)t dA
%
% INPUT:  ex,ey;       Element coordinates
%         t;           Element thickness	
%	  
%
% OUTPUT: Te :      Matrix 3 x 3
%-------------------------------------------------------------

Area=1/2*det([ones(3,1) ex' ey']);

L1=[0.5 0 0.5];
L2=[0.5 0.5 0];
L3=[0 0.5 0.5];

NtN=zeros(6);


for i=1:3
	NtN=NtN+1/3*[L1(i)^2 0 L1(i)*L2(i) 0 L1(i)*L2(i) 0
		 		0 L1(i)^2 0 L1(i)*L2(i) 0 L1(i)*L2(i)
		  		L2(i)*L1(i) 0 L2(i)^2 0 L2(i)*L3(i) 0
		  		0 L2(i)*L1(i) 0 L2(i)^2 0 L2(i)*L3(i)
		 		L3(i)*L1(i) 0 L3(i)*L2(i) 0 L3(i)^2 0
		  		0 L3(i)*L1(i) 0 L3(i)*L2(i) 0 L3(i)^2];
end

Te=NtN*Area*t;
     

