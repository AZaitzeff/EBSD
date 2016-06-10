function [smallf]=subres(f,m,n,radius)
sizef=size(f);
rowsize=sizef(1)/m;
colsize=sizef(2)/n;
smallf=zeros(m,n);
for i = 1:m,
    for j = 1:n,
        smallf(i,j)=circleaver(f,round(rowsize/2+(i-1)*rowsize),round(colsize/2+(j-1)*colsize),radius);
    end
end

function [averval]=circleaver(f,xpos,ypos,radius)
sizef=size(f);
mask=(repmat((1:sizef(2))-ypos,sizef(1),1).^2+repmat((1:sizef(1))'-xpos,1,sizef(2)).^2)<radius^2;
averval=sum(sum(mask.*f))/sum(sum(mask));