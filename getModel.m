function [MODEL] = getModel(g,In,M,m,l,b,xi,km)

% 9.81,0.099,2.4,0.23,0.4,0.05,0.005,8

global MODEL
MODEL = {};


syms s k1 k2 k3 k4
% g = 9.81;
% In = 0.099;
% M = 2.4;
% m = 0.23;
% l= 0.4;
% b= 0.05;
% xi = 0.005;
% km = 8;

delta = ((In*(m+M)) + (m*(l^2)*M));

%Terminos matrices:
a11 = -((xi*(m+M))/(delta));
a12 = ((b*m*l)/(delta));
a13 = ((m*g*l*(m+M))/(delta));
a14 = 0;

a21 = -((xi*m*l)/(delta));
a22 = -((b*(In+m*(l^2)))/((delta)));
a23 = -(((g*(m^2)*(l^2)))/(delta));
a24 = 0;

a31 = 1;
a32 = 0;
a33 = 0;
a34 = 0;

a41 =  0;
a42 = -1;
a43 =  0;
a44 =  0;




b11 = -((km*m*l)/(delta));
b21 = ((km*(In+(m*(l^2))))/(delta));
b31 = 0;
b41 = 0;

c11 = 1;
c12 = 0;
c13 = 0;
c14 = 0;

c21 =  0;
c22 = -1;
c23 =  0;
c24 =  0;

Ka = [k1,k2,k3,k4];

Aa= [a11 a12 a13 a14;
    a21 a22 a23 a24;
    a31 a32 a33 a34;
    a41 a42 a43 a44];

Ba = [b11;b21;b31;b41];

Ca = [c11 c12 c13 c14;
      c21 c22 c23 c24];
  
Da = [0;0];

MODEL.A = Aa;
MODEL.B = Ba;
MODEL.C = Ca;
MODEL.D = Da;

I = eye(4);  
  
  
At= (s*I) - (Aa)+(Ba*Ka);
    
poli = det(At);

% pretty(poli)


end

