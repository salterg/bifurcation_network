function [ f, front ] = network_frontfinder( eta, sealevel, x, nodes, j)
f=zeros(1,nodes);
for i=1:nodes;
    s=max(find(eta(i,:)>=sealevel(i,j)));
    if isempty(s)==1;
        f(i)=NaN;
        front(i)=NaN;
    else
        f(i)=s;
        front(i)=x(s,i);
    end
end
    
end




