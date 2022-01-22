function [ITERATIONS,niter]=initial_stress(inc_com,type_material,Kt_unco_fixed,Fint_unc,Fext_unc,grad_xsi_eta,Gauss_number,type_SF,lambda_step,ITERATIONS,GEOMETRY,MATERIAL)

%======================Classic Newton-Raphson iterations========================

    fprintf('INITIAL STRESS NEWTON-RAPHSON ITERATIONS: LAMBDA = %f \n\n',GEOMETRY.lambda_vect(lambda_step))

    niter=1;

    while ITERATIONS(niter).normResidual>GEOMETRY.norm_res_max && niter<=GEOMETRY.nitermax

        fprintf('Newton-Raphson iteration # %d \n',niter);

%============Calculate constrained system===============================================
        fprintf('Build constrained system... \n')
        [K,Fext,Fint,Delete]=build_constrained_sys(GEOMETRY,Kt_unco_fixed,Fint_unc,Fext_unc);
%=======================================================================================

%============Solve the system===============================================
        fprintf('Solve the system... \n\n')
        [ITERATIONS]=FEM_solver(K,Fext,Fint,Delete,niter,GEOMETRY,ITERATIONS);
        clear Delete
%===========================================================================
%=============Update SPATIAL data=====================================
        fprintf('Updating SPATIAL framework... \n')
        [ITERATIONS]=build_KINEMATICS(Gauss_number,niter+1,grad_xsi_eta,ITERATIONS,GEOMETRY);
%========================================================================
%============Update Material model=================================
         switch type_material
            case 'Neo_Hookean'
        fprintf('Updating Neo-Hookean model... \n')
        [ITERATIONS]=build_Heo_Hookean_model(inc_com,niter+1,Gauss_number,ITERATIONS,GEOMETRY,MATERIAL);
            case 'linear_elastic'
         fprintf('Updating linear elastic model... \n')    
         [ITERATIONS]=build_linear_elastic_model(type_SF,Gauss_number,niter+1,inc_com,ITERATIONS,GEOMETRY,MATERIAL);       
        end
%========================================================================   
%============Calculate B matricies=================================
         [ITERATIONS]=build_B(Gauss_number,type_SF,niter+1,ITERATIONS,GEOMETRY);
%===========================================================================    
%============Calculate internal forces======================================
        fprintf('Updating internal forces... \n')
        [ITERATIONS,Fint_unc]=build_F_int(Gauss_number,type_SF,niter+1,ITERATIONS,GEOMETRY);
%=========================================================================== 
%============Check the residual=============================================
       [~,Fext,Fint,~]=build_constrained_sys(GEOMETRY,Kt_unco_fixed,Fint_unc,Fext_unc);
       fprintf('Updating residual... \n\n') 
       ITERATIONS(niter+1).normResidual=norm(Fext-Fint);
%===========================================================================  
       niter=niter+1;

    end % END NR iterations

end % END fucntion 