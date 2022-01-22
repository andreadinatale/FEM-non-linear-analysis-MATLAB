function [Fext_unc]=build_F_ext(lambda_step,GEOMETRY)

Fext_unc=zeros(2*GEOMETRY.N_nodes,1);
for i=1:length(GEOMETRY.load(:,1))
Fext_unc((GEOMETRY.load(i,1)-1)*2+GEOMETRY.load(1,2))=GEOMETRY.lambda_vect(lambda_step)*GEOMETRY.load(i,3);
end