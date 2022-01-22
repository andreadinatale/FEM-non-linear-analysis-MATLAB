function [STEP]=build_Heo_Hookean_model(inc_com,lambda_step,Gauss_number,STEP,GEOMETRY,MATERIAL)

dir=@(i,j) 1.*(i==j);

for i=1:GEOMETRY.N_elem

%=============Calculate Neo_Hookean properties=====================================
     sigma=struct();
     for row=1:Gauss_number
         for column=1:Gauss_number 
             F=STEP(lambda_step).KINEMATICS(i).F(row,column).F;
             J=det(F);
             STEP(lambda_step).Neo_Hookean(i).mu_prime(row,column).mu_prime=MATERIAL.mu/J^2;
             STEP(lambda_step).Neo_Hookean(i).lambda_prime(row,column).lambda_prime=2*MATERIAL.mu/J^2;
             switch inc_com
                 case 'compressible'
                 STEP(lambda_step).Neo_Hookean(i).t(row,column).t=MATERIAL.T/J;  
                 case 'incompressible'
                 STEP(lambda_step).Neo_Hookean(i).t(row,column).t=MATERIAL.T;    
             end
             for ii=1:2
                  for j=1:2
                       for r=1:2
                            for s=1:2
                                 c(ii,j,r,s)=STEP(lambda_step).Neo_Hookean(i).lambda_prime(row,column).lambda_prime*dir(ii,j)*dir(r,s)+STEP(lambda_step).Neo_Hookean(i).mu_prime(row,column).mu_prime*(dir(ii,r)*dir(j,s)+dir(ii,s)*dir(j,r));
                            end
                       end
                  end
             end
             STEP(lambda_step).Neo_Hookean(i).C_Neo(row,column).C_Neo=[c(1,1,1,1) c(1,1,2,2) c(1,1,1,2)
                           c(2,2,1,1) c(2,2,2,2) c(2,2,1,2)
                           c(1,2,1,1) c(1,2,2,2) c(1,2,1,2)];
             sigma(row,column).sigma=MATERIAL.mu.*(F*F'-(1/J^2)*eye(2));
         end
     end
     STEP(lambda_step).Neo_Hookean(i).sigma=sigma;
%==================================================================================

end % END STEP cycle

end % END function