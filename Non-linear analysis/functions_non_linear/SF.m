function [vect_dN_xsi,vect_dN_eta]=SF(type_SF)

if type_SF==4    % S4 elementss
%==========Shape functions================

N1=@(xsi,eta) 0.25.*(1-xsi).*(1-eta);
N2=@(xsi,eta) 0.25.*(1+xsi).*(1-eta);
N3=@(xsi,eta) 0.25.*(1+xsi).*(1+eta);
N4=@(xsi,eta) 0.25.*(1-xsi).*(1+eta);

%==========Derivatives of shape functions w.r.t xsi=================

dN1_xsi=@(xsi,eta)  eta./4 - 1/4; 
dN2_xsi=@(xsi,eta)     1/4 - eta./4; 
dN3_xsi=@(xsi,eta)    eta./4 + 1/4;
dN4_xsi=@(xsi,eta)    - eta./4 - 1/4; 

%==========Derivatives of shape functions w.r.t eta=================
 
dN1_eta=@(xsi,eta)  xsi./4 - 1/4;
dN2_eta=@(xsi,eta)   - xsi./4 - 1/4;
dN3_eta=@(xsi,eta)   xsi./4 + 1/4; 
dN4_eta=@(xsi,eta)  1/4 - xsi./4;

SF=struct();
dSF_xsi=struct();
dSF_eta=struct();

%=============Structure of shape functions===============

SF.N1=N1;
SF.N2=N2;
SF.N3=N3;
SF.N4=N4;

%=============Structure of derivatives of shape functions w.r.t xsi===============

dSF_xsi.dN1_xsi=dN1_xsi;
dSF_xsi.dN2_xsi=dN2_xsi;
dSF_xsi.dN3_xsi=dN3_xsi;
dSF_xsi.dN4_xsi=dN4_xsi;

%=============Structure of derivatives of shape functions w.r.t eta===============

dSF_eta.dN1_eta=dN1_eta;
dSF_eta.dN2_eta=dN2_eta;
dSF_eta.dN3_eta=dN3_eta;
dSF_eta.dN4_eta=dN4_eta;

% xsi=[];
% eta=[];
% vect=[-1:0.1:1];
% xsi(1,:)=vect;
% eta(:,1)=flipud(vect');
% for i=2:length(vect)
%      xsi(i,:)=vect;
%      eta(:,i)=flipud(vect');
% end

% figure;
% surf(xsi,eta,N1(xsi,eta));
% title('N1');
% xlabel('xsi');
% ylabel('eta');
% 
% figure;
% surf(xsi,eta,N2(xsi,eta));
% title('N2');
% xlabel('xsi');
% ylabel('eta');
% 
% figure;
% surf(xsi,eta,N3(xsi,eta));
% title('N3');
% xlabel('xsi');
% ylabel('eta');
% 
% figure;
% surf(xsi,eta,N4(xsi,eta));
% title('N4');
% xlabel('xsi');
% ylabel('eta');

end %END if

if type_SF==8   % S8 elements

%==========Shape functions================

N1=@(xsi,eta)   0.25.*(1-xsi).*(1-eta).*(-xsi-eta-1);
N2=@(xsi,eta)   0.25.*(1+xsi).*(1-eta).*(xsi-eta-1);
N3=@(xsi,eta)   0.25.*(1+xsi).*(1+eta).*(xsi+eta-1);
N4=@(xsi,eta)   0.25.*(1-xsi).*(1+eta).*(-xsi+eta-1);
N5=@(xsi,eta)   0.5.*(1-xsi.^2).*(1-eta);      
N6=@(xsi,eta)   0.5.*(1+xsi).*(1-eta.^2);
N7=@(xsi,eta)   0.5.*(1-xsi.^2).*(1+eta);  
N8=@(xsi,eta)   0.5.*(1-xsi).*(1-eta.^2);

%==========Derivatives of shape functions w.r.t xsi=================

dN1_xsi=@(xsi,eta)   - (xsi./4 - 1/4).*(eta - 1) - ((eta - 1).*(eta + xsi + 1))./4;
dN2_xsi=@(xsi,eta)   ((eta - 1).*(eta - xsi + 1))./4 - (xsi./4 + 1/4).*(eta - 1);
dN3_xsi=@(xsi,eta)   (xsi./4 + 1/4).*(eta + 1) + ((eta + 1).*(eta + xsi - 1))./4;
dN4_xsi=@(xsi,eta)   (xsi./4 - 1/4).*(eta + 1) + ((eta + 1).*(xsi - eta + 1))./4;
dN5_xsi=@(xsi,eta)  xsi.*(eta - 1);
dN6_xsi=@(xsi,eta)  1/2 - eta.^2./2;
dN7_xsi=@(xsi,eta)   -xsi.*(eta + 1);
dN8_xsi=@(xsi,eta)   eta.^2./2 - 1/2; 

%==========Derivatives of shape functions w.r.t eta=================

dN1_eta=@(xsi,eta)   - (xsi./4 - 1/4).*(eta - 1) - (xsi./4 - 1/4).*(eta + xsi + 1);
dN2_eta=@(xsi,eta)   (xsi./4 + 1/4).*(eta - xsi + 1) + (xsi./4 + 1/4).*(eta - 1);
dN3_eta=@(xsi,eta)   (xsi./4 + 1/4).*(eta + 1) + (xsi./4 + 1/4).*(eta + xsi - 1);
dN4_eta=@(xsi,eta)   (xsi./4 - 1/4).*(xsi - eta + 1) - (xsi./4 - 1/4).*(eta + 1);
dN5_eta=@(xsi,eta)    xsi.^2./2 - 1/2;
dN6_eta=@(xsi,eta)    -2.*eta.*(xsi./2 + 1/2);
dN7_eta=@(xsi,eta)    1/2 - xsi.^2./2; 
dN8_eta=@(xsi,eta)   2.*eta.*(xsi./2 - 1/2);

SF=struct();
dSF_xsi=struct();
dSF_eta=struct();

%=============Structure of shape functions===============

SF.N1=N1;
SF.N2=N2;
SF.N3=N3;
SF.N4=N4;
SF.N5=N5;
SF.N6=N6;
SF.N7=N7;
SF.N8=N8;

%=============Structure of derivatives of shape functions w.r.t xsi===============

dSF_xsi.dN1_xsi=dN1_xsi;
dSF_xsi.dN2_xsi=dN2_xsi;
dSF_xsi.dN3_xsi=dN3_xsi;
dSF_xsi.dN4_xsi=dN4_xsi;
dSF_xsi.dN5_xsi=dN5_xsi;
dSF_xsi.dN6_xsi=dN6_xsi;
dSF_xsi.dN7_xsi=dN7_xsi;
dSF_xsi.dN8_xsi=dN8_xsi;

%=============Structure of derivatives of shape functions w.r.t eta===============

dSF_eta.dN1_eta=dN1_eta;
dSF_eta.dN2_eta=dN2_eta;
dSF_eta.dN3_eta=dN3_eta;
dSF_eta.dN4_eta=dN4_eta;
dSF_eta.dN5_eta=dN5_eta;
dSF_eta.dN6_eta=dN6_eta;
dSF_eta.dN7_eta=dN7_eta;
dSF_eta.dN8_eta=dN8_eta;

% xsi=[];
% eta=[];
% vect=[-1:0.1:1];
% xsi(1,:)=vect;
% eta(:,1)=flipud(vect');
% for i=2:length(vect)
%      xsi(i,:)=vect;
%      eta(:,i)=flipud(vect');
% end
% 
% figure;
% surf(xsi,eta,N1(xsi,eta));
% title('N1');
% xlabel('xsi');
% ylabel('eta');
% 
% figure;
% surf(xsi,eta,N2(xsi,eta));
% title('N2');
% xlabel('xsi');
% ylabel('eta');
% 
% figure;
% surf(xsi,eta,N3(xsi,eta));
% title('N3');
% xlabel('xsi');
% ylabel('eta');
% 
% figure;
% surf(xsi,eta,N4(xsi,eta));
% title('N4');
% xlabel('xsi');
% ylabel('eta');
% 
% figure;
% surf(xsi,eta,N5(xsi,eta));
% title('N5');
% xlabel('xsi');
% ylabel('eta');
% 
% figure;
% surf(xsi,eta,N6(xsi,eta));
% title('N6');
% xlabel('xsi');
% ylabel('eta');
% 
% figure;
% surf(xsi,eta,N7(xsi,eta));
% title('N7');
% xlabel('xsi');
% ylabel('eta');
% 
% figure;
% surf(xsi,eta,N8(xsi,eta));
% title('N8');
% xlabel('xsi');
% ylabel('eta');

end  %END if
%==============Cell array of shape functions and their derivatives=========
if type_SF==8
% vect_N={SF.N1 SF.N2 SF.N3 SF.N4 SF.N5 SF.N6 SF.N7 SF.N8};
vect_dN_xsi={dSF_xsi.dN1_xsi dSF_xsi.dN2_xsi dSF_xsi.dN3_xsi dSF_xsi.dN4_xsi dSF_xsi.dN5_xsi dSF_xsi.dN6_xsi dSF_xsi.dN7_xsi dSF_xsi.dN8_xsi};
vect_dN_eta={dSF_eta.dN1_eta dSF_eta.dN2_eta dSF_eta.dN3_eta dSF_eta.dN4_eta dSF_eta.dN5_eta dSF_eta.dN6_eta dSF_eta.dN7_eta dSF_eta.dN8_eta};
end
if type_SF==4
% vect_N={SF.N1 SF.N2 SF.N3 SF.N4};
vect_dN_xsi={dSF_xsi.dN1_xsi dSF_xsi.dN2_xsi dSF_xsi.dN3_xsi dSF_xsi.dN4_xsi};
vect_dN_eta={dSF_eta.dN1_eta dSF_eta.dN2_eta dSF_eta.dN3_eta dSF_eta.dN4_eta}; 
end
end %END function



