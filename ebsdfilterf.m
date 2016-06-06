function [newf]=ebsdfilterf(f,rows,cols,radius,sub)
sizef=size(f);
newf=zeros(sizef);
if sub==1,
    for i = rows,
        for j = cols,
            newf(i,j,:)=circleaver(f,i,j,radius);
        end
    end
else
    for i = 1:sizef(1),
        for j = 1:sizef(2),
            newf(i,j,:)=circleaver(f,i,j,radius);
        end
    end
end
function [averval]=circleaver(f,xpos,ypos,radius)
sizef=size(f);
mask=(repmat((1:sizef(2))-ypos,sizef(1),1).^2+repmat((1:sizef(1))'-xpos,1,sizef(2)).^2)<radius^2;
total=sum(sum(mask));
mask=repmat(mask,1,1,sizef(3));
averval=sum(sum(mask.*f))/total;