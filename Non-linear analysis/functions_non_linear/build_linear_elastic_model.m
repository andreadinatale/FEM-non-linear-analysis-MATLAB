function [STEP]=build_linear_elastic_model(type_SF,Gauss_number,lambda_step,inc_com,STEP,GEOMETRY,MATERIAL)


%=============Define MATERIAL linear-elastic struct====================================
for i=1:GEOMETRY.N_elem

sigma=struct();
Almansi=struct();
     for row=1:Gauss_number
         for column=1:Gauss_number 
             F=STEP(lambda_step).KINEMATICS(i).F(row,column).F;
             J=det(F);  
             almansi=0.5.*(eye(2)-F'\inv(F));
             almansi_voigt=[almansi(1,1) almansi(2,2) 2*almansi(1,2)]';
             switch inc_com
                 case 'compressible'
                 STEP(lambda_step).Neo_Hookean(i).t(row,column).t=MATERIAL.T/J;
                 case 'incompressible'
                 STEP(lambda_step).Neo_Hookean(i).t(row,column).t=MATERIAL.T;
             end

                 STEP(lambda_step).Neo_Hookean(i).C_Neo(row,column).C_Neo=(MATERIAL.E/(1-MATERIAL.nu^2)).*[1 MATERIAL.nu 0
                                    MATERIAL.nu 1 0
                                    0 0 (1-MATERIAL.nu)/2];
                 Sigma=STEP(lambda_step).Neo_Hookean(i).C_Neo(row,column).C_Neo*almansi_voigt;
                 sigma(row,column).sigma=[Sigma(1) Sigma(3)
                     Sigma(3) Sigma(2)];
                 Almansi(row,column).Almansi=almansi;

         end
     end
     STEP(lambda_step).Neo_Hookean(i).Almansi=Almansi;
     STEP(lambda_step).Neo_Hookean(i).sigma=sigma;
end



end % END function
