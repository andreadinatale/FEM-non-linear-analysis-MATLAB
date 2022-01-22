function [STEP,Fint_unc]=build_F_int(Gauss_number,type_SF,lambda_step,STEP,GEOMETRY)

%==============Vector of Gauss weigth==============
if Gauss_number==3
w=GEOMETRY.int_rule.three_point.w;
end
if Gauss_number==2
w=GEOMETRY.int_rule.two_point.w;
end
if Gauss_number==1
w=GEOMETRY.int_rule.one_point.w;
end
%===================================================

Fint_unc=zeros(2*GEOMETRY.N_nodes,1);

for i=1:GEOMETRY.N_elem

%=================Internal forces============================= 
    
     fint=zeros(2*type_SF,1);
     for row=1:Gauss_number
         for column=1:Gauss_number
             t=STEP(lambda_step).Neo_Hookean(i).t(row,column).t;
             B=STEP(lambda_step).K_F(i).B(row,column).B;
             Sigma=STEP(lambda_step).Neo_Hookean(i).sigma(row,column).sigma;     
             Sigma=[Sigma(1,1);Sigma(2,2);Sigma(1,2)];
             j=STEP(lambda_step).KINEMATICS(i).j(row,column).j;
             fint=fint+(B'*Sigma).*(t*det(j)*w(row)*w(column));
         end
     end
     STEP(lambda_step).K_F(i).fint=fint;

%============================================================
     
     fint_unc=GEOMETRY.fint;
     fint_unc(GEOMETRY.ptrs(i,:))=fint;
     Fint_unc=Fint_unc+fint_unc;

end % END STEP cycle     

end % END function