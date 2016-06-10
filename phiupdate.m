function [uo] = phiupdate(nt,dt,uin,f,fid,g1,g2,radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [uo] = phiupdate(nt,dt,uin,f,fid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Piecewise constant Mumford-Shah model based on Chan and Vese.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input:

%   nt = Number of time steps to take.

%   dt = Time step size.

%   uin = Initial level set function.

%   f = Given image to be segmented.

%   fid = Fidelity constant.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output:

%   uo = Final level set function.

%   g1 = Best unknown constant 1 for uo.

%   g2 = Best unknown constant 2 for uo.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[m,n] = size(uin); %size of level set function
[ms,ns,~] = size(f);
factor=m/ms;
me = 0.01;  % Regularization of denominators.
esp=.01;

u = buffer2(uin);   % Size is now (m+4)x(n+4). nice works around the boundaries
innerg1g2f=zeros(ms,ns);
g1g2=g1-g2;

for i = 1:ms,
    for j = 1:ns,
        innerg1g2f(i,j)=dot(g1g2,reshape(f(i,j,:),[1,4800]));
    end
end 
bigg1g2f=ebsdfilter(innerg1g2f,m,n,radius,factor,me);
g1g2g1=dot(g1g2,g1);
g1g2g2=dot(g1g2,g2);




for t = 1:nt % Main time loop
    % Determine the best choice of g1,g2 based on current u TO DO:

    %g1 = sum(sum( double(u>0).*fb ))/(sum(sum( double(u>0) )) + 0.0001);

    %g2 = sum(sum( double(u<0).*fb ))/(sum(sum( double(u<0) )) + 0.00001);


    % Difference quotients:  

    u_x = forwardx(u)*m;

    u_y = forwardy(u)*n;  

    Du = sqrt( u_x.^2 + u_y.^2 + me );
    
   
    %H(u,.01),m,n,radius,1,0)
    GGH=ebsdfilter(subres(H(u(3:m+2,3:n+2),esp),ms,ns,radius),m,n,radius,factor,me);
    %GGH2=ebsdfilter(subres(1-H(u(3:m+2,3:n+2),esp),ms,ns,radius),m,n,radius,factor);

    X=GGH*g1g2g1+(1-GGH)*g1g2g2-bigg1g2f;
    %X(i,j)=H1inner(g1-g2,g1*GGH(i,j)+g2*(GGH2(i,j))-vecf(Gfb(i,j),g1,g2,noise(i,j)),me);
    
    Xb=buffer2(X);
    rhs = backwardx(u_x./Du)*m + backwardy(u_y./Du)*n- fid*(Xb);
    
    u = u + dt*DH(u,esp).*rhs;
    

end



uo = u(3:m+2,3:n+2);


function [inner]=L2inner(vec1,vec2,me)

inner=dot(vec1,vec2)/(norm(vec1)*norm(vec2)+me);

function [value]=H1inner(vec1,vec2,me)
image1=reshape(vec1,[80,60]);
image2=reshape(vec2,[80,60]);
im1f=fftshift(fft2(image1));
im2f=fftshift(fft2(image2));
thecumsum=0;
for n=1:80,
    for m=1:60,
        if n~=41 && m~=31,
            thecumsum=thecumsum+im1f(n,m)*im2f(n,m)/((n-41)^2+(m-31)^2+me);
        end
    end
end
value=real(thecumsum);

function [newvec]=vecf(val,g1,g2,noise)

newvec=val*g1+g2*(1-val)+noise;

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



function [ub]=buffer2(u)



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



function [dy]=backwardy(u)



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



function [dx]=backwardx(u)



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





function [dy]=forwardy(u)



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





function [dx]=forwardx(u)



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



function y = DH(x,me)



d = 0.0001;



y = (H(x+d,me)-H(x-d,me))/(2*d);



function y = H(x,me)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function y = H(x,me)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Approximate Heaviside function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



y = (tanh(x/me)+1)/2;