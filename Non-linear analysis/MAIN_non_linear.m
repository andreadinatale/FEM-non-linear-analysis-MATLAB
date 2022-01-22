%% NON LINEAR FE MODEL

clear all
close all
clc

%% IMPORT GEOMETRY WITH MESH 

addpath("functions_non_linear")
addpath("Matlab FEM Input files")

%==========CHOOSE ANALYSIS YOU WANT TO PERFORM==========================
GEOMETRY=inputNL_cantilever; 
GEOMETRY.dlambda=50;  
Gauss_number=3;
Ampl_factor=1;
type_NR='classicNR';
inc_com='compressible';
type_material='Neo_Hookean';
%=======================================================================
Size=size(GEOMETRY.elements);
type_SF=Size(2);
clear Size


subplot(3,1,1);
x_matrix=GEOMETRY.nodes(:,1);
y_matrix=GEOMETRY.nodes(:,2);
z_matrix=zeros(size(GEOMETRY.nodes(:,1)));
plot(x_matrix,y_matrix,'k*'), grid on, hold on;
num_ele=size(GEOMETRY.elements);
title('Amplification factor: ',num2str(Ampl_factor),'interpreter','latex');
xlabel('x (mm)','fontsize',15,'interpreter','latex');
ylabel('y (mm)','fontsize',15,'interpreter','latex');
axis equal

%=============Define MATERIAL linear-elastic struct====================================
MATERIAL=struct();
MATERIAL.E=GEOMETRY.E;
MATERIAL.nu=GEOMETRY.nu;
MATERIAL.G=MATERIAL.E/(2*(1+MATERIAL.nu));
MATERIAL.T=GEOMETRY.t;
MATERIAL.mu=MATERIAL.G;
%=======================================================================

%===================CREATE GAUSS POINTS=================================
GEOMETRY.int_rule=struct();
GEOMETRY.int_rule.one_point=struct();
GEOMETRY.int_rule.one_point.x=[0];
GEOMETRY.int_rule.one_point.w=[2];
GEOMETRY.int_rule.two_point=struct();
GEOMETRY.int_rule.two_point.x=[-1/sqrt(3),1/sqrt(3)];
GEOMETRY.int_rule.two_point.w=[1,1];
GEOMETRY.int_rule.three_point=struct();
GEOMETRY.int_rule.three_point.x=[-sqrt(0.6),0,sqrt(0.6)];
GEOMETRY.int_rule.three_point.w=[5/9,8/9,5/9];
%======================================================================

GEOMETRY.N_elem=num_ele(1);
num_nodes=size(GEOMETRY.nodes);
GEOMETRY.N_nodes=num_nodes(1);
clear num_nodes num_ele

%% RESHAPING THE ELEMENTS AND THE NODES

[GEOMETRY]=reshaping(GEOMETRY,type_SF);

%% IMPORT SHAPE FUNCTIONS AND THEIR DERIVATIVES

[vect_dN_xsi,vect_dN_eta]=SF(type_SF);

%% INIZIALIZE STRUCTURES

STEP=struct();
STEP.Neo_Hookean=struct();
STEP.KINEMATICS=struct();
STEP.K_F=struct();
STEP(1).Displ=zeros(GEOMETRY.N_nodes,2);
grad_xsi_eta=struct();

%% CALCULATE GRADIENT W.R.T XSI AND ETA

[grad_xsi_eta]=build_grad_xsi_eta(Gauss_number,vect_dN_xsi,vect_dN_eta,grad_xsi_eta,type_SF,GEOMETRY);
clear vect_dN_eta vect_dN_xsi 

%% ITERATE OVER LOAD STEPS

for lambda_step=1:length(GEOMETRY.lambda_vect)  % Perfomr an incremental procedure, increasing each time the loads applied
   
    fprintf('LAMBDA = %f \n\n',GEOMETRY.lambda_vect(lambda_step))
    
%=============Calculate KINEMATICS data=====================================
    fprintf('Calculating KINEMATICS... \n')
    [STEP]=build_KINEMATICS(Gauss_number,lambda_step,grad_xsi_eta,STEP,GEOMETRY);
%========================================================================
%============Calculate Neo-Hookean model=================================
switch type_material
    case 'Neo_Hookean'
    fprintf('Calculating Neo-Hookean model... \n')
    [STEP]=build_Heo_Hookean_model(inc_com,lambda_step,Gauss_number,STEP,GEOMETRY,MATERIAL);
    case 'linear_elastic'
    fprintf('Calculating linear elastic model... \n')
    [STEP]=build_linear_elastic_model(type_SF,Gauss_number,lambda_step,inc_com,STEP,GEOMETRY,MATERIAL);
end
%========================================================================    
%============Calculate B matricies=================================   
    [STEP]=build_B(Gauss_number,type_SF,lambda_step,STEP,GEOMETRY);
%============Calculate internal forces======================================
    fprintf('Calculating internal forces... \n')
    [STEP,Fint_unc]=build_F_int(Gauss_number,type_SF,lambda_step,STEP,GEOMETRY);
%===========================================================================    
%============Calculate external forces======================================
    fprintf('Calculating external forces... \n')
    [Fext_unc]=build_F_ext(lambda_step,GEOMETRY);
%===========================================================================   
    fprintf('Calculating the residual before NR iterations... \n\n')
    STEP(lambda_step).normResidual=norm(Fext_unc-Fint_unc);

    ITERATIONS=struct();
    ITERATIONS(1).normResidual=norm(Fext_unc-Fint_unc);
    ITERATIONS(1).KINEMATICS=STEP(lambda_step).KINEMATICS;
    ITERATIONS(1).Neo_Hookean=STEP(lambda_step).Neo_Hookean;
    ITERATIONS(1).K_F=STEP(lambda_step).K_F;

%% NEWTON-RAPHSON

switch type_NR

    case 'classicNR'

%============Calculate stiffneess matricies=================================   
    fprintf('Calculating stiffness matricies... \n')
    [STEP,Kt_unco]=build_K(Gauss_number,type_SF,lambda_step,STEP,GEOMETRY);
%=========================================================================== 
[ITERATIONS,niter]=classicNR(inc_com,type_material,Kt_unco,Fint_unc,Fext_unc,grad_xsi_eta,Gauss_number,type_SF,lambda_step,ITERATIONS,GEOMETRY,MATERIAL);

    case 'modifiedNR'

%============Calculate stiffneess matricies=================================   
    fprintf('Calculating stiffness matricies... \n')
    [STEP,Kt_unco]=build_K(Gauss_number,type_SF,lambda_step,STEP,GEOMETRY);
%=========================================================================== 
[ITERATIONS,niter]=modifiedNR(inc_com,type_material,Kt_unco,Fint_unc,Fext_unc,grad_xsi_eta,Gauss_number,type_SF,lambda_step,ITERATIONS,GEOMETRY,MATERIAL);

    case 'initial_stress'
        if lambda_step==1
            %============Calculate stiffneess matricies=================================   
                   fprintf('Calculating stiffness matricies... \n')
                   [STEP,Kt_unco_fixed]=build_K(Gauss_number,type_SF,lambda_step,STEP,GEOMETRY);
            %=========================================================================== 
        end
[ITERATIONS,niter]=initial_stress(inc_com,type_material,Kt_unco_fixed,Fint_unc,Fext_unc,grad_xsi_eta,Gauss_number,type_SF,lambda_step,ITERATIONS,GEOMETRY,MATERIAL);

end % END switch type NR

Total_displ=zeros(GEOMETRY.N_nodes,2);
for i=2:niter
    Total_displ=Total_displ+ITERATIONS(i).Displ;
end 



STEP(lambda_step+1).Displ=Total_displ;
Total_displ=zeros(GEOMETRY.N_nodes,2);
for i=1:lambda_step+1
    Total_displ=Total_displ+STEP(i).Displ;
end
subplot(3,1,2)
plot(max(abs(Total_displ(:,1))),GEOMETRY.lambda_vect(lambda_step),'r*'), grid on, hold on;
title('Displacement-Load','interpreter','latex');
xlabel('Max U (mm)','fontsize',15,'interpreter','latex');
ylabel('\lambda','fontsize',15);
axis([0 inf 0 GEOMETRY.lambda_vect(end)])
subplot(3,1,3)
plot(max(abs(Total_displ(:,2))),GEOMETRY.lambda_vect(lambda_step),'r*'), grid on, hold on;
title('Displacement-Load','interpreter','latex');
xlabel('Max V (mm)','fontsize',15,'interpreter','latex');
ylabel('\lambda','fontsize',15);
axis([0 inf 0 GEOMETRY.lambda_vect(end)])
clear Total_displ
STEP(lambda_step).ITERATIONS_NR=ITERATIONS;
clear ITERATIONS

end  % END load step

%% POST PROCESSING

fprintf('POST-PROCESSING \n\n')

[STEP]=build_KINEMATICS(Gauss_number,lambda_step+1,grad_xsi_eta,STEP,GEOMETRY);

switch type_material
    case 'Neo_Hookean'
[STEP]=build_Heo_Hookean_model(inc_com,lambda_step+1,Gauss_number,STEP,GEOMETRY,MATERIAL);
    case 'linear_elastic'
[STEP]=build_linear_elastic_model(type_SF,Gauss_number,lambda_step+1,inc_com,STEP,GEOMETRY,MATERIAL);
sigma=struct();
if Gauss_number==3
x=GEOMETRY.int_rule.three_point.x;
end
if Gauss_number==2
x=GEOMETRY.int_rule.two_point.x;
end
if Gauss_number==1
x=GEOMETRY.int_rule.one_point.x;
end
for i=1:GEOMETRY.N_elem
     for row=1:Gauss_number
         for column=1:Gauss_number 
             F=STEP(lambda_step+1).KINEMATICS(i).F(row,column).F;
             J=det(F);  
             almansi=0.5.*(eye(2)-F'\inv(F));
             almansi_voigt=[almansi(1,1) almansi(2,2) 2*almansi(1,2)]';
             t=STEP(lambda_step+1).Neo_Hookean(i).t(row,column).t;
             Sigma=STEP(lambda_step+1).Neo_Hookean(i).C_Neo(row,column).C_Neo.*t*almansi_voigt;
             sigma(row,column).sigma=[Sigma(1) Sigma(3)
                     Sigma(3) Sigma(2)];
         end
     end
     STEP(lambda_step+1).Neo_Hookean(i).sigma=sigma;
     fprintf('Calculating stress distributions on the %d element... \n',i);
    sigma_xx_vect=[];
    sigma_yy_vect=[];
    sigma_xy_vect=[];
    Matrix=[];
    for row=1:Gauss_number
        for column=1:Gauss_number
                 stress_Gauss=[sigma(row,column).sigma(1,1), sigma(row,column).sigma(2,2), sigma(row,column).sigma(1,2)];
                 sigma_xx_vect=[sigma_xx_vect; stress_Gauss(1)];
                 sigma_yy_vect=[sigma_yy_vect; stress_Gauss(2)];
                 sigma_xy_vect=[sigma_xy_vect; stress_Gauss(3)];
                 if Gauss_number==2
                     Matrix=[Matrix;1, x(column), x(row), x(column)*x(row)];
                 end
                 if Gauss_number==3
                     Matrix=[Matrix;1, x(column), x(row), x(column)*x(row), x(column)^2, x(row)^2, x(column)^2*x(row), x(column)*x(row)^2, x(column)^2*x(row)^2];
                 end
        end
    end
    a_xx=Matrix\sigma_xx_vect;
    a_yy=Matrix\sigma_yy_vect;
    a_xy=Matrix\sigma_xy_vect;
    if Gauss_number==2
        STEP(lambda_step+1).Neo_Hookean(i).sigma_xx=@(xsi,eta) a_xx(1)+a_xx(2).*xsi+a_xx(3).*eta+a_xx(4).*xsi.*eta;
        STEP(lambda_step+1).Neo_Hookean(i).sigma_yy=@(xsi,eta) a_yy(1)+a_yy(2).*xsi+a_yy(3).*eta+a_yy(4).*xsi.*eta;
        STEP(lambda_step+1).Neo_Hookean(i).sigma_xy=@(xsi,eta) a_xy(1)+a_xy(2).*xsi+a_xy(3).*eta+a_xy(4).*xsi.*eta;
    end
    if Gauss_number==3
        STEP(lambda_step+1).Neo_Hookean(i).sigma_xx=@(xsi,eta) a_xx(1)+a_xx(2).*xsi+a_xx(3).*eta+a_xx(4).*xsi.*eta+a_xx(5).*xsi.^2+a_xx(6).*eta.^2+a_xx(7).*xsi.^2.*eta+a_xx(8).*xsi.*eta.^2+a_xx(9).*xsi.^2.*eta.^2;
        STEP(lambda_step+1).Neo_Hookean(i).sigma_yy=@(xsi,eta) a_yy(1)+a_yy(2).*xsi+a_yy(3).*eta+a_yy(4).*xsi.*eta+a_yy(5).*xsi.^2+a_yy(6).*eta.^2+a_yy(7).*xsi.^2.*eta+a_yy(8).*xsi.*eta.^2+a_yy(9).*xsi.^2.*eta.^2;
        STEP(lambda_step+1).Neo_Hookean(i).sigma_xy=@(xsi,eta) a_xy(1)+a_xy(2).*xsi+a_xy(3).*eta+a_xy(4).*xsi.*eta+a_xy(5).*xsi.^2+a_xy(6).*eta.^2+a_xy(7).*xsi.^2.*eta+a_xy(8).*xsi.*eta.^2+a_xy(9).*xsi.^2.*eta.^2;
    end
end
end

Total_displ=zeros(GEOMETRY.N_nodes,2);
for i=1:length(GEOMETRY.lambda_vect)+1
    Total_displ=Total_displ+STEP(i).Displ;
end

subplot(3,1,1)

plot(x_matrix+Ampl_factor*Total_displ(:,1),y_matrix+Ampl_factor*Total_displ(:,2),'bo','LineWidth',3)
legend('Nodes (Initial)','Nodes (Current)','Location','northeast')











