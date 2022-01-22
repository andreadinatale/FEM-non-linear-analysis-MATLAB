function [STEP]=build_B(Gauss_number,type_SF,lambda_step,STEP,GEOMETRY)

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

end % END STEP cycle

end % END function