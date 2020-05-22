%physical parameters

r=.5;%parameter in slope term of cross stream sed flux feedback
lamp=.4;%bed porosity
g=9.81*1;%gravity 
ps=1*2650;%particle density kg/m^3
p=1*1000;%density of water
Rr = (ps-p)/p;%submerged specific gravity of sediment
Cz1=15;%upstream dimensionless chezy - only for const Chezy, otherwise set by grain size
timeyear = 31557600;%seconds in a year

%fundamental dimensionless variables (imposed at upstream node)
taus1=.06;
asp_rat=100;
S1=.001;
%numerical parameters
dt_years=5*10^-5;
dt = dt_years * timeyear;%dt in seconds

timesteps=2000000;
finaltime=timesteps*dt_years;%final time in years
perturb_amp=2*10^-4;% adding a little noise every time step to the bed elevations
perturbtype=2;%0=no perturbation,1=each timestep,2=once per run at specific time
perturb_int=1;%every x timesteps apply perturbation (for pertubtype==1)
timestep_perturb=2;
years_perturb=timestep_perturb*dt_years;%year at which perturb (type 2) is applied
watlev_accuracy=1*10^-12;%how accurately you iteratively solve for the constancy of water surface elevation
recovery=.9999;%for recovering from complete avulsion within while loop
complete_threshold=1*10^-8;%lowest possible discharge before a complete avulsion occurs
t=linspace(0,finaltime,timesteps);%time in years
saveinterval=100;
completesaveinterval=1000;

%% switches

transport_form=1;%1=MPM,2=EH,3=other
down_bc=3;%set downstream b.c.
    %1=equal water surface elevation, 2=equal bed elevation, 3=sed flux
nonconst_chezy=0;%1 to have chezy coefficent which varies with depth
sealevel_mode=0;%1 to have time varying sea level in each branch
channelbelt=0;%0 =width is channel width, 1= specified width, 2= specified angle (radial)
reset_topo=1;%reset current topography? 1=yes, 0=no, 
H123_set=0;%if 0, H1 is substituded for H123
xi_down=0;%if 0 (default) xi computed only with eta(1), if 1, use (eta(1)+eta(2))/2
makeplots=0;
makemovie=0;


%% set switch specific parameters
%sediment transport formula

if nonconst_chezy==0;
    Cz=Cz1*ones(1,nodes);
elseif nonconst_chezy==1;
    Cz1=6+2.5*log(H1/(2.5*ds));  
end
    if transport_form==1;
        n=4;
        m=1.5;
        tausc=.047;
    elseif transport_form==2;
        n=.05*Cz1^2;
        m=2.5;
        tausc=0;
    else transport_form==3;
        n=1;
        m=.1;
        tausc=0;
    end 
%downstream boundary condition
    if down_bc==3;
    F_a=[NaN NaN NaN 0 0 0 0];
    F_i=[NaN NaN NaN .7 .7 .7 .7];
    end
%sea level, specified in m
if sealevel_mode==0;
    sealevel=-inf*ones(nodes,timesteps);%if sealevel is < bed elevation it has no effect
elseif sealevel_mode==1;
    sea_amp=[NaN NaN NaN 10 10 0 0];
    sea_period=[NaN NaN NaN 10 100 1 1];
    sea_initial=[NaN NaN NaN 5 5 5 5];
    sea_rate=[NaN NaN NaN .01 0 0 1];%in m/year
    sealevel=sea_initial'*ones(1,length(t))+sea_rate'*t + (sin(2*pi*(1./sea_period)'*t)).*(sea_amp'*ones(1,length(t)));
    sealevel(2,:)=max(sealevel(4,:),sealevel(5,:));
    sealevel(3,:)=max(sealevel(6,:),sealevel(7,:));
    sealevel(1,:)=max(sealevel(2,:),sealevel(3,:));
end 
%channel belt
if channelbelt==1;
    B_belt=B*2;
elseif channelbelt==0;
    B_belt=B;
else channelbelt==2;
%radial expansion - not yet implemented
end

%calculated physical parameters
H1=B1/asp_rat;
ds=S1*H1/(Rr*taus1);
Q1=sqrt(((H1^2)*((B1^2)*(Cz1^2)*taus1*(ps-p)*g*ds))/p);
Qs1=(n*(taus1-tausc).^m)*((Rr * g * ds) ^ 0.5) * ds*B1;


