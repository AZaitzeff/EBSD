% Read in the (built in) image:
%test = load('test.mat');
%f = test.f(1:100,1:100);
f=(repmat((1:400)-200,400,1).^2+repmat((1:400)'-200,1,400).^2)<100^2;
sizef=size(f);
sigma=20;
rows=10:20:400;
cols=10:20:400;
testf=ebsdfilter(f,rows,cols,sigma,1);

%dict = matfile('Dict_woprj_new333226nobg.mat');
%g1=dict.dict_woprj(1,:);
%g2=dict.dict_woprj(10,:);
%create f
%newf=zeros(20,20,4800);
%for i=1:20,
%    for j=1:20,
%       newf(i,j,:)=g1*testf(i,j)+g2*(1-testf(i,j))+normrnd(1,10,[1,4800]);
%    end
%end
% Look at it:
imagesc(f)
colormap gray

% Prepare initial level set function (a line):
u0 = zeros(sizef)-1;
u0(125:275,125:275)=1;
%for i=1:sizef(1)
%    for j=1:sizef(2)
%        if 250*i>j^2,
%            u0(i,j) =1;
%        end
%    end
%end
uin = RSreinit2D(1000,1/(5*500),u0);
%noise=zeros([20,20,1,4800]);
% Look at it on the image:
%imagesc(f)
%colormap gray(256)
%hold on
%contour(uin,[0 0],'r')

% Run the program with that initial condition:
%[u] = phiupdatesc(100,1/(5*100^2),uin,testf,100,sigma,rows,cols)
%[u] = phiupdate(100,1/(5*10^2),u0,Gfb,20,g1,g2,sigma,1);
Gfb=ebsdfilter(testf,rows,cols,sigma,0);
[u] = phiupdatesc(50,1/(5*100^2),uin,Gfb,10000,sigma,rows,cols);
save('u1.mat', 'u', '-v7.3');
[u] = phiupdatesc(50,1/(5*100^2),u,Gfb,10000,sigma,rows,cols);
save('u2.mat', 'u', '-v7.3');
[u] = phiupdatesc(50,1/(5*100^2),u,Gfb,10000,sigma,rows,cols);
save('u3.mat', 'u', '-v7.3');
[u] = phiupdatesc(50,1/(5*100^2),u,Gfb,10000,sigma,rows,cols);
save('u4.mat', 'u', '-v7.3');
[u] = phiupdatesc(50,1/(5*100^2),u,Gfb,10000,sigma,rows,cols);
save('u5.mat', 'u', '-v7.3');
[u] = phiupdatesc(50,1/(5*100^2),u,Gfb,10000,sigma,rows,cols);
save('u6.mat', 'u', '-v7.3');

% Look at the result:
%imagesc(f)
%colormap gray(256)
%hold on
%contour(uin,[0 0],'b')
%contour(u,[0 0],'g')
% Run some more:
%[u] = phiupdate(10,1/(5*20^2),u,testf,100,g1,g2,noise,sigma);

save('Dict_woprj_new333226nobg.mat', 'dict_woprj', '-v7.3');
% Look at the result:
%imagesc(f)
%colormap gray(256)
%hold on
%contour(round(Gfb/.0032)-.5,[0 0],'y')