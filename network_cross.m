function [Q_cross, Qs_cross]=network_cross(Q_cross, Qs_cross, eta,Q,Qs, H,dx_up,dx,B,Cz,g,ds,alpha,r,taus1, H123_set,Adjacency,Ncells);

for i=[2,4,6];
    f=find(Adjacency(i,:)==1);
    B1=B(f);
    B2=B(i);
    B3=B(i+1);

    H2=H(2,i);
    H3=H(2,i+1);
    S32=(eta(1,i+1)-eta(1,i))/((B2+B3)/2);

    Qs1=Qs(Ncells(f),f);
    H1=H(Ncells(f)+1,f);
    

    Q2=Q(i);
    Q3=Q(i+1);
    Q1=Q(f);

    [Q23,Qs23]=cross(B1,B2,B3,S32,H2,H3,H1,alpha,Qs1,r,taus1,H123_set,Q1,Q2,Q3);
    Q_cross(i)=-Q23;%convention - flux in is positive
    Q_cross(i+1)=Q23;
    Qs_cross(1,i)=-Qs23;
    Qs_cross(1,i+1)=Qs23;
    
end

end

function [Q23,Qs23]=cross(B1,B2,B3,S32,H2,H3,H1,alpha,Qs1,r,taus1,H123_set,Q1,Q2,Q3)

          Q23=(B2/(B3+B2))*Q1-Q2;     
       
        
        if H123_set==1;
        H123=.5*((H2(2)+H3(2))/2+H1);
        else
            H123=H1;
        end
       
        Qs23=alpha*B1*(Qs1/B1)*(Q23*H1/(Q1*alpha*H123)-(r/sqrt(taus1))*S32);

        if isnan(Qs23)==1;
            Qs23=0;
        end
end