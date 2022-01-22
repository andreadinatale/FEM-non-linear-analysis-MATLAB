function [grad_xsi_eta]=build_grad_xsi_eta(Gauss_number,vect_dN_xsi,vect_dN_eta,grad_xsi_eta,type_SF,GEOMETRY)


%==============Vector of Gauss point==============

if Gauss_number==3
x=GEOMETRY.int_rule.three_point.x;
end
if Gauss_number==2
x=GEOMETRY.int_rule.two_point.x;
end
if Gauss_number==1
x=GEOMETRY.int_rule.one_point.x;
end

%===============================================================
for row=1:Gauss_number
         for column=1:Gauss_number
              if type_SF==8
                  grad_xsi_eta(row,column).grad_xsi_eta=[[vect_dN_xsi{1}(x(column),x(row)) vect_dN_xsi{2}(x(column),x(row)) vect_dN_xsi{3}(x(column),x(row)) vect_dN_xsi{4}(x(column),x(row)) vect_dN_xsi{5}(x(column),x(row)) vect_dN_xsi{6}(x(column),x(row)) vect_dN_xsi{7}(x(column),x(row)) vect_dN_xsi{8}(x(column),x(row))]',[vect_dN_eta{1}(x(column),x(row)) vect_dN_eta{2}(x(column),x(row)) vect_dN_eta{3}(x(column),x(row)) vect_dN_eta{4}(x(column),x(row)) vect_dN_eta{5}(x(column),x(row)) vect_dN_eta{6}(x(column),x(row)) vect_dN_eta{7}(x(column),x(row)) vect_dN_eta{8}(x(column),x(row))]'];
              end
              if type_SF==4
                  grad_xsi_eta(row,column).grad_xsi_eta=[[vect_dN_xsi{1}(x(column),x(row)) vect_dN_xsi{2}(x(column),x(row)) vect_dN_xsi{3}(x(column),x(row)) vect_dN_xsi{4}(x(column),x(row))]',[vect_dN_eta{1}(x(column),x(row)) vect_dN_eta{2}(x(column),x(row)) vect_dN_eta{3}(x(column),x(row)) vect_dN_eta{4}(x(column),x(row))]'];
              end
         end
end



end % END function 