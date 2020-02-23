function [M,pos1,pos2,Xigual] =  proximidade(x1,y1,x2,y2)
% x1=xlateralEsqSup(indexEsq);
% y1=yLateralEsq(indexEsq);
% x2=xblob;
% y2=yblob_menor;

if(y1>0)
    d=zeros(size(x1,2),size(x2,2));
    for i=1:size(x1,2)
        for j=1:size(x2,2)
            d(i,j) = sqrt((x1(i)-x2(j)).^2.+(y1(i)-y2(j)).^2);
        end
    end
else
    d=zeros(size(x1,2),size(x2,2));
    for i=1:size(x1,2)
        for j=1:size(x2,2)
            d(i,j) = sqrt((x1(i)-x2(j)).^2);
        end
    end
end

[M,I]=min(d(:));
[pos1, pos2] = ind2sub(size(d),I);

Xbinario=x1(pos1)==x2;
[~,Xigual]=max(Xbinario);
%[x1(pos1),y1(pos1);[x2(pos2),y2(pos2)] ;[x2(Xigual),y2(Xigual)]]

