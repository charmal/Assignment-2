function [p] = convAdvec(h,sigma,mu,q)

%   Input:
%           h: scalar. cellsize
%       sigma: vector. spatially varying coefficient, evaluated at cell centres.
%          mu: scalar. constant coefficient
%           q: vector. Source function, evaluated at nodes.
%  
%   Output: 
%           p: Approximate value of p on the nodes.
%           
%   You must write your discretization as a matrix equation A*p = q.
%   Once you have formed the matrix A, you can solve for p using the 
%   Matlab backslash operator. Make sure that you form A as a sparse
%   matrix. If A is dense the backslash operator will be very inefficient.

% Size of Vector
n = length(sigma);
h = 1/n;

%Centres and Nodes 
Xc = [h/2:h:(1-(h/2))]';
Xn = [h:h:1-h]';

%Sigma
sig = @(x) 1+(x).^2;
sigma = sig(Xc);
Sig = spdiags(sigma,0,n,n);

% Derivative Matrices

   % Nodes to Cell
    Dnc = (1/h)*spdiags([ones(n-1,1), -ones(n-1,1)],[0,-1],n,n-1);
    
    % Cell to nodes
    Dcn = (1/h)*spdiags([-ones(n-1,1), ones(n-1,1)],[0,1],n-1,n);
    
    % Nodes to Nodes
    Dnn = (1/(2*h))*spdiags([-ones(n-1,1), ones(n-1,1)],[-1,1],n-1,n-1);
    
% Define A
A = (Dcn)*(Sig)*(Dnc)+(mu)*(Dnn);

q= (2*pi*(cos(2*pi*Xn).*(2*Xn+mu)-2*pi*sin(2*pi*Xn).*(1+Xn.^2)));
 
p = A\q;

end