function [STEP]=build_KINEMATICS(Gauss_number,lambda_step,grad_xsi_eta,STEP,GEOMETRY)

for i=1:GEOMETRY.N_elem

     el=GEOMETRY.nodes(GEOMETRY.elements(i,:),:);
     STEP(lambda_step).KINEMATICS(i).XY_nodes=el;
     
%=============Define Jacobian J at each Gauss point ===================
     
     coord_nodes=el';
     J_struct=struct();
     for row=1:Gauss_number
         for column=1:Gauss_number
              
             J_struct(row,column).J=coord_nodes*grad_xsi_eta(row,column).grad_xsi_eta;
              
         end
     end
     STEP(lambda_step).KINEMATICS(i).J=J_struct;

%=============Define Gradient in MATERIAL framework===================

     GradN=struct();
     for row=1:Gauss_number
         for column=1:Gauss_number
             J=J_struct(row,column).J;
           
             GradN(row,column).GradN=J'\grad_xsi_eta(row,column).grad_xsi_eta';
             
         end
     end
     STEP(lambda_step).KINEMATICS(i).GradN=GradN;
     
%=====================================================================
%=====================================================================
%=============Calculate SPATIAL coordinates===============================

     displ=STEP(lambda_step).Displ(GEOMETRY.elements(i,:)',:);
     if lambda_step==1
     STEP(lambda_step).KINEMATICS(i).xy_nodes=STEP(lambda_step).KINEMATICS(i).XY_nodes+displ;
     else
     STEP(lambda_step).KINEMATICS(i).xy_nodes=STEP(lambda_step-1).KINEMATICS(i).xy_nodes+displ;
     end
     
%=========================================================================

%=============Calculate j at Gauss points=================================  
     coord_nodes=STEP(lambda_step).KINEMATICS(i).xy_nodes';
     j_struct=struct();
     for row=1:Gauss_number
         for column=1:Gauss_number
            
             j_struct(row,column).j=coord_nodes*grad_xsi_eta(row,column).grad_xsi_eta;
              
         end
     end
     STEP(lambda_step).KINEMATICS(i).j=j_struct;
%=========================================================================

%=============Calculate gradient in SPATIAL coordinates at Gauss points=================================
     gradN=struct();
     for row=1:Gauss_number
         for column=1:Gauss_number
             j=j_struct(row,column).j;
            
                 gradN(row,column).gradN=j'\grad_xsi_eta(row,column).grad_xsi_eta';
             
         end
     end
     STEP(lambda_step).KINEMATICS(i).gradN=gradN;
%=======================================================================================================

%=============Calculate F at Gauss points=================================
     F=struct();
     for row=1:Gauss_number
         for column=1:Gauss_number   
             F(row,column).F=STEP(lambda_step).KINEMATICS(i).xy_nodes'*(STEP(lambda_step).KINEMATICS(i).GradN(row,column).GradN)';
         end
     end
     STEP(lambda_step).KINEMATICS(i).F=F;
%=========================================================================

end % END STEP cycle

end % END function 