

if reset_topo==1;
    startslope=1*10^-3;
    eta=-x*startslope;
    eta=eta-min(min((eta)));
    for i=1:nodes;
        eta(Ncells(i)+1:max(Ncells),i)=NaN;
    end
end

Q=zeros(1,nodes);
for i=1:nodes;
    if sum(Adjacency(i,:))==0;
    Q(i)=Q1;
    else
        f1=find(Adjacency(i,:)==1);
        if length(f1)==1; %not confluence
            %Q(i)=Q(f1)/sum(Adjacency(:,f1));
            f2=find(Adjacency(:,f1)==1);
            f2=f2(find(f2~=i));
            Q(i)=Q(f1)*B(i)/(B(i)+B(f2));
        end
    end
end

%initialize model vectors

%initialize data vectors
H=zeros(max(Ncells)+1,nodes);
H(1,1)=H1;
taus=zeros(max(Ncells)+1,nodes);
taus(1,1)=taus1;
Qs=zeros(max(Ncells)+1,nodes);
Qs(1,1)=Qs1;

Q_cross=zeros(1,nodes);
Qs_cross=zeros(1,nodes);
%data collected every timestep
Q_save=zeros(nodes,timesteps/saveinterval);
Qs_save=zeros(nodes,timesteps/saveinterval);
eta_save=zeros(nodes,timesteps/saveinterval);
etad_save=zeros(nodes,timesteps/saveinterval);
if sealevel_mode==1;
front_save=zeros(nodes,timesteps/saveinterval);
end
if nonconst_chezy==1;
Cz_save=zeros(nodes,timesteps/saveinterval);
end
t_save=linspace(0,finaltime,timesteps/saveinterval);
mass_initial=(dx.*B_belt)*sum(eta(2:max(Ncells),:),1)'+(dx_up.*B_belt)*eta(1,:)'*(1-lamp);

eta_completesave=zeros(max(Ncells),nodes,timesteps/completesaveinterval);
Qs_completesave=zeros(max(Ncells)+1,nodes,timesteps/completesaveinterval);;


