function y = RSreinit2D(nt,dt,uin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Russo-Smereka reinitialization.
% nt = Number of time steps to take.
% dt = Time step size. Explicit scheme.
% uin = Level set function to reinitialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

me = 1e-6; % Regularization.

[m,n] = size(uin);  % Size of the domain.

u = buffer2(uin);   % Buffer up the input array.
u0geq0 = double(u >= 0);
u0le0 = 1- u0geq0;

S = u./sqrt(u.^2+me);   % Regularized sign function.

% Allocation:
b = zeros(m+4,n+4);
a = b;
d = b;
c = b;



f = sign(u);
mb = m+4;
nb = n+4;
dfb = abs(centerx(f))+abs(centery(f));
dfb = double(dfb>0);

dfb(1,1:n+4)=0;
dfb(m+4,1:n+4)=0;
dfb(1:m+4,1)=0;
dfb(1:m+4,n+4)=0;



df = reshape(dfb,1,mb*nb);
ind = find(df);
% Convert the given distance function matrix into a long row vector:
R = reshape(u,1,mb*nb);
% Calculate the gradient at positions indicated by ind:
R_x = zeros(1,mb*nb);
R_x(ind) = (R(ind+nb)-R(ind-nb))*n/2;
R_y = zeros(1,mb*nb);
R_y(ind) = (R(ind+1)-R(ind-1))*m/2;
DR = sqrt( R_x.^2 + R_y.^2 );
% Compute the distance function at these positions:
D = zeros(1,mb*nb);
D(ind) = R(ind)./max(DR(ind),0.01);
D = reshape(D,mb,nb);



u = D;





for t=1:nt
    
    % Difference quotients:
    b = forwardx(u)*m;
    a(3:m+2,3:n+2) = b(3:m+2,2:n+1);
    a = upbuffer2(a);
    d = forwardy(u)*n;
    c(3:m+2,3:n+2) = d(2:m+1,3:n+2);
    c = upbuffer2(c);
    
    am = min(0,a);
    ap = max(0,a);
    bm = min(0,b);
    bp = max(0,b);
    cm = min(0,c);
    cp = max(0,c);
    dm = min(0,d);
    dp = max(0,d);
    
    G = sqrt( max(ap.^2,bm.^2) + max(cp.^2,dm.^2) ) - 1;
    G = G.*u0geq0 +...
        ( sqrt( max(am.^2,bp.^2) + max(cm.^2,dp.^2) ) - 1 ).*u0le0;
    
    % Update:
    %u = u - dt*S.*G;
    u = u - dt*(1-dfb).*S.*G;
    
end

y = u(3:m+2,3:n+2);

function [ub]=buffer2(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function buffers up a given matrix with a buffer layer
% of thickness 2. It uses periodic boundary conditions to fill
% the buffer layer.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses: upbuffer2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: If the input matrix u has size mxn, then the output 
%    matrix ub has size (m+4)x(n+4).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(u);

ub=zeros(m+4,n+4);  % Initialize ub.
ub(3:m+2,3:n+2)=u;  % Bulk part of ub is just u.

ub=upbuffer2(ub);   % Fills in buffer layer via periodicity.

function [ubo]=upbuffer2(ub)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function updates the buffer layer of a buffered matrix, 
% using periodicity, where the buffer layer thickness is 2. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: No error checking is built in. Make sure size of input
%    matrix is large enough: both dims. > 4 needed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=size(ub,1)-4;
n=size(ub,2)-4;

% Update the buffer:

ub(3:m+2,1:2)=ub(3:m+2,n+1:n+2);
ub(3:m+2,n+3:n+4)=ub(3:m+2,3:4);
ub(1:2,1:n+4)=ub(m+1:m+2,1:n+4);
ub(m+3:m+4,1:n+4)=ub(3:4,1:n+4);

ubo=ub;

function [dx]=forwardx(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward difference operator in the x-direction for the input 
% matrix u, which is assumed to be buffered (layer thickness 2). 
% Periodic boundary conditions used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses: upbuffer2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(u);  % Buffered size.
m=m-4;      % Find the unbuffered size.
n=n-4;      % Find the unbuffered size.

dx=u;       % Initialization for dx.

% Carry out operation on unbuffered part of u:
dx(3:m+2,3:n+2)=( u(3:m+2,4:n+3) - u(3:m+2,3:n+2) );

% Update the buffer region:
dx=upbuffer2(dx);

function [dy]=forwardy(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward difference operator in the y-direction for the input
% matrix u, which is assumed to be buffered (layer thickness 2).
% Periodic boundary conditions used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses: upbuffer2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(u);  % Buffered size.
m=m-4;      % Find the unbuffered size.
n=n-4;      % Find the unbuffered size.

dy=u;       % Initialization for dy.

% Carry out operation on unbuffered part of u:
dy(3:m+2,3:n+2)=( u(4:m+3,3:n+2) - u(3:m+2,3:n+2) );

% Update the buffer region:
dy=upbuffer2(dy);

function [dx]=centerx(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Centered difference operator in the x-direction for the input
% matrix u, which is assumed to be buffered (layer thickness 2).
% Periodic boundary conditions used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses: upbuffer2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(u);  % Buffered size.
m=m-4;      % Find the unbuffered size.
n=n-4;      % Find the unbuffered size.

dx=u;       % Initialization for dx.

% Carry out operation on unbuffered part of u:
dx(3:m+2,3:n+2)=( u(3:m+2,4:n+3) - u(3:m+2,2:n+1) )/2;

% Update the buffer region:
dx=upbuffer2(dx);

function [dy]=centery(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Centered difference operator in the y-direction for the input
% matrix u, which is assumed to be buffered (layer thickness 2).
% Periodic boundary conditions used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses: upbuffer2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(u);  % Buffered size.
m=m-4;      % Find the unbuffered size.
n=n-4;      % Find the unbuffered size.

dy=u;       % Initialization for dy.

% Carry out the operation on unbuffered part of u:
dy(3:m+2,3:n+2)=( u(4:m+3,3:n+2) - u(2:m+1,3:n+2) )/2;

% Update the buffer region:
dy=upbuffer2(dy);

function [dx]=backwardx(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backward difference operator in the x-direction for the input
% matrix u, which is assumed to be buffered (layer thickness 2).
% Periodic boundary conditions used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses: upbuffer2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(u);  % Buffered size.
m=m-4;      % Find the unbuffered size.
n=n-4;      % Find the unbuffered size.

dx=u;       % Initialization for dx.

% Carry out operation on unbuffered part of u:
dx(3:m+2,3:n+2)=( u(3:m+2,3:n+2) - u(3:m+2,2:n+1) );

% Update the buffer region:
dx=upbuffer2(dx);

function [dy]=backwardy(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backward difference operator in the y-direction for the input
% matrix u, which is assumed to be buffered (layer thickness 2).
% Periodic boundary conditions used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses: upbuffer2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(u);  % Buffered size.
m=m-4;      % Find the unbuffered size.
n=n-4;      % Find the unbuffered size.

dy=u;       % Initialization for dy.

% Carry out the operation on unbuffered part of u: 
dy(3:m+2,3:n+2)=( u(3:m+2,3:n+2) - u(2:m+1,3:n+2) );

% Update the buffer region:
dy=upbuffer2(dy);
