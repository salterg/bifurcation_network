
run network_geometry
run network_params
run network_initialize

for j=1:timesteps
    

    if perturbtype==1&rem(j,perturb_int)==0||perturbtype==2&j==years_perturb/dt_years
            eta=eta+perturb_amp*randn(size(eta));
    end
    
    if down_bc==1
    A1=[1 -1 0 0; 0 0 1 -1; .5 .5 -.5 -.5; 1 1 1 1];
    b1=[eta(1,4)-eta(1,5);eta(1,6)-eta(1,7);eta(1,2)-eta(1,3);0];
    b2=A1\b1;
    for i=1:4
        eta(Ncells(i+3),i+3)=b2(i);
    end
    elseif down_bc==2
    eta(Ncells(4:7),4:7)=0;
    end
    
    
        
        f=zeros(1,nodes);
        for i=1:nodes;
        s=max(find(eta(:,i)>sealevel(i,j)));
        if isempty(s)==1;
        f(i)=0;
        front(i)=NaN;
        else
        f(i)=s;
        front(i)=x(s,i);
        end
        end
        
        [ Q, H ] = network_flowsplit(eta,Q,H,dx_up,dx,B,Cz,g,nonconst_chezy,xi_down,watlev_accuracy,recovery,complete_threshold,ds,Adjacency,Ncells);
        %fill in the rest of H vectors;
        H(2,1)=(Q(1)^(2/3))./((B(1)^(2/3))*(Cz(1)^(2/3))*(((eta(1,1)-max(sealevel(1,j),eta(2,1)))/(mean([dx(1),dx_up(1)]))).^(1/3))*(g^(1/3)));
        for i=1:nodes;
            if f(i)<Ncells(i)
            H(3:f(i)+1,i)=(Q(i)^(2/3))./((B(i)^(2/3))*(Cz(i)^(2/3))*(((eta(2:f(i),i)-max(sealevel(i,j),eta(3:f(i)+1,i)))/dx(i)).^(1/3))*(g^(1/3)));
            else
                H(3:f(i),i)=(Q(i)^(2/3))./((B(i)^(2/3))*(Cz(i)^(2/3))*(((eta(2:f(i)-1,i)-max(sealevel(i,j),eta(3:f(i),i)))/dx(i)).^(1/3))*(g^(1/3)));
            end
            H(f(i)+2:Ncells(i),i)=0;%depth=0 beyond shoreline
            fdown=find(Adjacency(:,i)==1);
            if isempty(fdown)==0;
            H(Ncells(i)+1,i)=(Q(i)^(2/3))./((B(i)^(2/3))*(Cz(i)^(2/3))*(((eta(Ncells(i),i)-max(sealevel(i,j),mean([eta(1,fdown(1)),eta(1,fdown(2))])))/(mean([mean([dx_up(fdown(1)),dx_up(fdown(2))]),dx(i)]))).^(1/3))*(g^(1/3)));
            end
        end
        
        %compute stresses
        for i=1:nodes;
                taus(2:f(i)+1,i)=p*((Q(i)).^2)./((ps-p)*g*ds*(Cz(i)^2)*(B(i)^2)*(H(2:f(i)+1,i).^2));
                taus(f(i)+2:Ncells(i),i)=0;
                fdown=find(Adjacency(:,i)==1);
                if isempty(fdown)==0;
                    taus(Ncells(i)+1,i)=p*((Q(i)).^2)./((ps-p)*g*ds*(Cz(i)^2)*(B(i)^2)*(H(Ncells(i)+1,i).^2));
                end
        end

        %compute sediment fluxes
        for i=1:nodes;
            for k=2:Ncells(i);%loop to compute sediment flux in downstream branches
               if taus(k,i)>tausc;
                   Qs(k,i) = (n*(taus(k,i)-tausc).^m)*((Rr * g * ds) ^ 0.5) * ds*B(i);
               else
                   Qs(k,i)=0;
               end
            end
            fdown=find(Adjacency(:,i)==1);
            if isempty(fdown)==0;
                if taus(Ncells(i)+1,i)>tausc;
                Qs(Ncells(i)+1,i)=(n*(taus(Ncells(i)+1,i)-tausc).^m)*((Rr * g * ds) ^ 0.5) * ds*B(i);
                else
                 Qs(Ncells(i)+1,i)=0;
                end
            end
            
        end
        for i=1:nodes;
            fup=find(Adjacency(i,:)==1);
            if isempty(fup)==0;
                if i==2|i==4|i==6
                Qs(1,i)=Qs(Ncells(fup)+1,fup)*B(i)./(B(i)+B(i+1)); 
                elseif i==3|i==5|i==7;
                Qs(1,i)=Qs(Ncells(fup)+1,fup)*B(i)./(B(i)+B(i-1)); 
                end
            end
        end
        
        if down_bc==3;
            for i=1:nodes;
                fdown=find(Adjacency(:,i)==1);
                fup=find(Adjacency(i,:)==1);
                if isempty(fdown)==1;
                    
                    Qs(Ncells(i)+1,i)=F_i(i)*Qs(Ncells(i),i)./(F_i(i)+(dx(i)*(sum(B(4:7)))/sum(L.*B))*((1-F_i(i)))); 
                end
            end
        end
        
        %cross-stream fluxes
        [Q_cross, Qs_cross]=network_cross(Q_cross, Qs_cross, eta,Q,Qs, H,dx_up,dx,B,Cz,g,ds,alpha,r,taus1, H123_set,Adjacency,Ncells);
        
        %update bed elevation
        for i=1:nodes
            eta(1,i)=eta(1,i)+(1/(1-lamp))*(Qs(1,i)-Qs(2,i))*dt/(dx_up(i)*B_belt(i))+dt*(1/(1-lamp))*(Qs_cross(i))/(dx_up(i)*B_belt(i));%alpha*B1 reach
            eta(2:Ncells(i),i)=eta(2:Ncells(i),i)+(1/(1-lamp))*dt*(Qs(2:Ncells(i),i)-Qs(3:Ncells(i)+1,i))./(dx(i).*B_belt(i));
        end
        
        
        if rem(j, saveinterval)==0;
        Q_save(:,j/saveinterval)=Q';
        eta_save(:,j/saveinterval)=eta(1,:)';
        Qs_save(:,j/saveinterval)=Qs(2,:)';
        etad_save(:,j/saveinterval)=nansum(eta(2:end,:),1);
        end
        if rem(j,completesaveinterval)==0;
            eta_completesave(:,:,j/completesaveinterval)=eta;
           
            Qs_completesave(:,:,j/completesaveinterval)=Qs;
        end
         if rem(j,500)==1;
            j
        end

        if makemovie==1;
            if rem(j,vid_interval)==0;
                fig=figure;
                subplot(2,1,1);
                plot(x,eta-551);
                axis([0 x(end, end) 0 50])
                subplot(2,1,2);
                plot(t,Q_save(:,4),t,Q_save(:,5),t,Q_save(:,6),t,Q_save(:,7))
                
            end
        end
    
end
t_save_star=t_save*(1/(1-lamp))*timeyear*(1-Fnew)*(Qs(1,1)/sum(L.*B))/H(1,1);
