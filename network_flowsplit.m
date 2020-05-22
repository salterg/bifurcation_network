function [ Q, H ] = network_flowsplit(eta,Q,H,dx_up,dx,B,Cz,g,nonconst_chezy,xi_down,watlev_accuracy,recovery,complete_threshold,ds,Adjacency,Ncells);
for i=[2,4,6];
    
    S2=(eta(1,i)-eta(2,i))/(.5*dx(i)+.5*dx_up(i));
    S3=(eta(1,i+1)-eta(2,i+1))/(.5*dx(i+1)+.5*dx_up(i+1));
    eta2=eta(1,i);
    eta3=eta(1,i+1);
    eta2d=eta(2,i);
    eta3d=eta(2,i+1);
    Q1=Q(find(Adjacency(i,:)==1));
    Q2=Q(i);
    Q3=Q(i+1);
    B2=B(i);
    B3=B(i+1);
    Cz1=Cz(Adjacency(i,:)==1);
    Cz2=Cz(i);
    Cz3=Cz(i+1);
    
    [Q2,Q3,H2,H3]=compute_Q(S2,S3,eta2,eta3,eta2d,eta3d,Q1,Q2,Q3,B2,B3,Cz1,Cz2,Cz3,g,nonconst_chezy,xi_down,watlev_accuracy,recovery,ds,complete_threshold);
    Q(i)=Q2;
    Q(i+1)=Q3;
    H(2,i)=H2;
    H(2,i+1)=H3;
    
    
end

end

function [Q2,Q3,H2,H3]=compute_Q(S2,S3,eta2,eta3,eta2d,eta3d,Q1,Q2,Q3,B2,B3,Cz1,Cz2,Cz3,g,nonconst_chezy,xi_down,watlev_accuracy,recovery,ds,complete_threshold);
    if isnan(Q2)==1;
        Q2=Q1/2;
    end
    if isnan(Q3)==1;
        Q3=Q1/2;
    end
    Qold1=Q2+Q3;
    if Qold1>0;
    Q2=Q2*Q1/Qold1;
    Q3=Q3*Q1/Qold1;
    else
        Q2=Q1/2;
        Q3=Q1/2;
    end
    if S2<=0;

        Q2=0;
        Q3=Q1;
        H2=0;
        H3=(Q3^(2/3))/((B3^(2/3))*(Cz1^(2/3))*(S3^(1/3))*(g^(1/3)));


    elseif S3<=0;

        Q2=Q1;
        Q3=0;
        H3=0;
        H2=(Q2^(2/3))/((B2^(2/3))*(Cz1^(2/3))*(S2^(1/3))*(g^(1/3)));


    else

H2=(Q2^(2/3))/((B2^(2/3))*(Cz2^(2/3))*(S2^(1/3))*(g^(1/3)));%normal flow
H3=(Q3^(2/3))/((B3^(2/3))*(Cz3^(2/3))*(S3^(1/3))*(g^(1/3)));
if nonconst_chezy==1;
Cz2=6+2.5*log(H2/(2.5*ds));
Cz3=6+2.5*log(H3/(2.5*ds));
end

xi2=H2+eta2(1);%water surface elevation at cell boundary, iterated until xi2=xi3;
xi3=H3+eta3(1);
    if xi_down==1;
            xi2=H2+.5*(eta2+eta2d);%water surface elevation at cell boundary, iterated until xi2=xi3;
            xi3=H3+.5*(eta3+eta3d);
    end
    


        while abs(xi2-xi3)>watlev_accuracy;%iteration to find fluxes Q2 and Q3 such that xi2=xi3

            if Q2==0;
                Q2=Q1-recovery*Q1;
                Q3=Q1-Q2;
            elseif Q3==0;
                Q3=Q1-recovery*Q1;
                Q2=Q1-Q3;
            end
            %improve estimate of Q2 and Q3
            %newton raphson iteration method
            Q2=Q2-(3/2)*(xi2-xi3)/((Q2^(-1/3))/((B2^(2/3))*(Cz2^(2/3))*(S2^(1/3))*(g^(1/3)))+((Q1-Q2)^(-1/3))/((B3^(2/3))*(Cz3^(2/3))*(S3^(1/3))*(g^(1/3))));%newton raphson iteration
            Q3=Q1-Q2;
            %if discharge is below a threshold, complete avulsion
            %occurs and exit while loop
            if Q2<complete_threshold;
                Q2=0;
                Q3=Q1;
                H2=0;
                H3=(Q3^(2/3))/((B3^(2/3))*(Cz3^(2/3))*(S3^(1/3))*(g^(1/3)));
                if nonconst_chezy==1;
                    error('complete avulsion')
                    exit
                end
                break
            end
            if Q3<complete_threshold;
                Q3=0;
                Q2=Q1;
                H3=0;
                H2=(Q2^(2/3))/((B2^(2/3))*(Cz2^(2/3))*(S2^(1/3))*(g^(1/3)));
                if nonconst_chezy==1;
                    error('complete avulsion')
                    exit
                end
                break
            end
            %after finding Q2 and Q3
            H2=(Q2^(2/3))/((B2^(2/3))*(Cz2^(2/3))*(S2^(1/3))*(g^(1/3)));%normal flow
            H3=(Q3^(2/3))/((B3^(2/3))*(Cz3^(2/3))*(S3^(1/3))*(g^(1/3)));

            if nonconst_chezy==1;
            Cz2=6+2.5*log(H2/(2.5*ds));
            Cz3=6+2.5*log(H3/(2.5*ds));
            end

            xi2=H2+eta2;%water surface elevation at cell boundary, iterated until xi2=xi3;
            xi3=H3+eta3;
            if xi_down==1;
                xi2=H2+.5*(eta2+eta2d);%water surface elevation at cell boundary, iterated until xi2=xi3;
                xi3=H3+.5*(eta3+eta3d);
            end


        end
    end
end