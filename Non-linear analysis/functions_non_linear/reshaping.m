function [GEOMETRY]=reshaping(GEOMETRY,type_SF)

fprintf('Reshaping... \n');

% if type_SF==4
%     for i=1:GEOMETRY.N_elem
%         el=[GEOMETRY.elements(i,4) GEOMETRY.elements(i,3) GEOMETRY.elements(i,1) GEOMETRY.elements(i,2)];
%         GEOMETRY.elements(i,:)=el;
%     end
% end
% if type_SF==8
%     for i=1:GEOMETRY.N_elem
%         el=[GEOMETRY.elements(i,4) GEOMETRY.elements(i,7) GEOMETRY.elements(i,3) GEOMETRY.elements(i,8) GEOMETRY.elements(i,6) GEOMETRY.elements(i,1) GEOMETRY.elements(i,5) GEOMETRY.elements(i,2)];
%         GEOMETRY.elements(i,:)=el;
%     end
% end

GEOMETRY.K_unc=zeros(GEOMETRY.N_nodes*2,GEOMETRY.N_nodes*2);
GEOMETRY.F_unc=zeros(GEOMETRY.N_nodes*2,1);
GEOMETRY.Fhat=zeros(GEOMETRY.N_nodes*2,1);
GEOMETRY.fint=zeros(GEOMETRY.N_nodes*2,1);
GEOMETRY.res=zeros(GEOMETRY.N_nodes*2,1);
GEOMETRY.lambda_vect=[GEOMETRY.dlambda:GEOMETRY.dlambda:GEOMETRY.lambda_max];


GEOMETRY.pos=[[1:2:GEOMETRY.N_nodes*2-1]',[2:2:GEOMETRY.N_nodes*2]'];   % Position of dof's of each node in the global vector of unkonwns 

for i=1:GEOMETRY.N_elem
    el=GEOMETRY.pos([GEOMETRY.elements(i,:)],:);
    el_n=[];
    for j=1:length(el)
        el_n=[el_n el(j,:)];
    end
    GEOMETRY.ptrs(i,:)=el_n;     % Position of dof's of each element
end 

fprintf('Reshaping is over \n\n');

end % END function

