function [K,Fext,Fint,Delete]=build_constrained_sys(GEOMETRY,Kt_unco,Fint_unc,Fext_unc)

%=========Evaluate the constrained dof's================
Delete=[];    % vector that contains the constrained dof's
for i=1:length(GEOMETRY.spc(:,1))
    
    if GEOMETRY.spc(i,2)==1
        
        Delete(i)=(GEOMETRY.spc(i,1)-1)*2+1;
    end
    if GEOMETRY.spc(i,2)==2
       
        Delete(i)=(GEOMETRY.spc(i,1)-1)*2+2;
    end

end
%=======================================================

%==========Cancel out the rows and coloums related to constrained dof's====
Kt_unco(Delete,:)=[];
Kt_unco(:,Delete)=[];
Fint_unc(Delete)=[];
Fext_unc(Delete)=[];

K=Kt_unco;
Fext=Fext_unc;
Fint=Fint_unc;

end  % END function 
