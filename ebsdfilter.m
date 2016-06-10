function [newf]=ebsdfilter(f,m,n,radius,factor)

newf=zeros(m,n);

for i = 1:m,
    for j = 1:n,
        newf(i,j)=circleaver(f,i,j,radius,factor);
    end
end

function [averval]=circleaver(f,xpos,ypos,radius,factor)
sizef=size(f);
mask=(repmat((1:sizef(2))-(ypos/factor+1/2),sizef(1),1).^2+repmat((1:sizef(1))'-(xpos/factor+1/2),1,sizef(2)).^2)<(radius/factor)^2;
averval=sum(sum(mask.*f))/sum(sum(mask));