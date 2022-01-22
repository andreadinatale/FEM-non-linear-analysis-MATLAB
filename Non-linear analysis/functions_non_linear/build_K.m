function [STEP,Kt_unco]=build_K(Gauss_number,type_SF,lambda_step,STEP,GEOMETRY)

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
%======================================================================================================

Kt_unco=zeros(2*GEOMETRY.N_nodes,2*GEOMETRY.N_nodes);

for i=1:GEOMETRY.N_elem

%==================Calculate B matrix at Gauss points=================================    

     B=struct();
     for row=1:Gauss_number
         for column=1:Gauss_number
             B_vect=[];
             for j=1:type_SF
                 grad=STEP(lambda_step).KINEMATICS(i).gradN(row,column).gradN(:,j);  
                 B_alpha=[grad(1) 0;0 grad(2);grad(2) grad(1)];
                 B_vect(:,[2*j-1 2*j])=B_alpha;
             end
             B(row,column).B=B_vect;
         end
     end
     STEP(lambda_step).K_F(i).B=B;

%====================================================================================

%====================Material stiffness=======================================

     Kc=zeros(2*type_SF,2*type_SF);
     for row=1:Gauss_number
         for column=1:Gauss_number
             t=STEP(lambda_step).Neo_Hookean(i).t(row,column).t;
             B=STEP(lambda_step).K_F(i).B(row,column).B;
             C=STEP(lambda_step).Neo_Hookean(i).C_Neo(row,column).C_Neo;
             j=STEP(lambda_step).KINEMATICS(i).j(row,column).j;
             Kc=Kc+(B'*C*B).*(det(j)*w(row)*w(column)*t);
         end
     end
     STEP(lambda_step).K_F(i).Kc=Kc;

%=============================================================================     

%=================Geometrical stiffness=======================================     
  
     Kg=zeros(2*type_SF,2*type_SF);
     for row=1:Gauss_number
         for column=1:Gauss_number
             t=STEP(lambda_step).Neo_Hookean(i).t(row,column).t;
             Sigma=STEP(lambda_step).Neo_Hookean(i).sigma(row,column).sigma;     
             grad=STEP(lambda_step).KINEMATICS(i).gradN(row,column).gradN; 
             j=STEP(lambda_step).KINEMATICS(i).j(row,column).j;
             Kg_small=(grad'*Sigma*grad).*(det(j)*w(row)*w(column)*t);
             Kg=Kg+kron(Kg_small,eye(2));
         end
     end
     STEP(lambda_step).K_F(i).Kg=Kg;

%=============================================================================  

     Kt=Kg+Kc;
     K_unc_local=GEOMETRY.K_unc;
     K_unc_local(GEOMETRY.ptrs(i,:),GEOMETRY.ptrs(i,:))=Kt;
     Kt_unco=Kt_unco+K_unc_local;
 
end % END STEP cycle

end % END function 