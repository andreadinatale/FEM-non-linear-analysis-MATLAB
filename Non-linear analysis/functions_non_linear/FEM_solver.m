function [STEP]=FEM_solver(K,Fext,Fint,Delete,lambda_step,GEOMETRY,STEP)

F=Fext-Fint;
q=K\F;
q_uncon=zeros(2*GEOMETRY.N_nodes,1);
non_zero_pos=[1:2*GEOMETRY.N_nodes];
non_zero_pos(Delete)=[];
q_uncon(non_zero_pos)=q;
for i=1:GEOMETRY.N_nodes
     q_uncon_matrix(i,:)=[q_uncon(2*i-1) q_uncon(2*i)];
end
STEP(lambda_step+1).Displ=q_uncon_matrix;
end