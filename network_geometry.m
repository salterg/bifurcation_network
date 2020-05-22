nodes=7;
Adjacency=zeros(nodes,nodes);
Adjacency([2,3],1)=1;%row i and column j, node j flows into nodes i
Adjacency([4,5],2)=1;
Adjacency([6,7],3)=1;
B=[100 70 70 50 50 50 50];%width of channel downstream of each node

alpha=3;
L=3*B.*[6 11 11 23 23 23 23];%length of each branch downstream of node

Ncells=1*[20 40 40 80 80 80 80];

dx_up=zeros(1,length(nodes));
dx_up(1)=L(1)/Ncells(1);
dx_up(2:3)=alpha*B(1);
dx_up(4:7)=alpha*B(2);
dx=(L-dx_up)./(Ncells-1);


%x vector
x=zeros(max(Ncells),nodes);%matrix of coordinate for cell centers
xface=zeros(max(Ncells)+1,nodes);%matrix of coordinate for cell faces
x(1,1)=dx_up(1)/2;
xface(1,1)=0;
x(2,1)=x(1,1)+(dx(1)+dx_up(1))/2;
xface(2,1)=dx_up(1);
for i=2:nodes;
    f_up=find(Adjacency(i,:)==1);
    xmin=L(f_up);
    f_up2=find(Adjacency(f_up,:)==1);
    if isempty(f_up2)==0;
        xmin=L(f_up)+L(f_up2); 
    end
    x(1,i)=dx_up(i)/2+xmin;
    xface(1,i)=0+xmin;
    x(2,i)=x(1,i)+(dx(i)+dx_up(i))/2;
    xface(2,i)=dx_up(i)+xmin;

end

for j=1:nodes
for k=3:Ncells(j);
    x(k,j)=x(k-1,j)+dx(j);
    xface(k,j)=xface(k-1,j)+dx(j);
end
xface(Ncells(j)+1,:)=xface(Ncells(j),:)+dx(j);
end

B1=B(1);