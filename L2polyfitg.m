function [fpol, rd] = L2polyfitg(dm,m,My,m1,m2)
%L2POLYFITG Finds the L2 least-squares polynomial to match given data
%   dm  data to fit the polynomial to
%   number of 2x2 subsystems
%   My  order of the polynomial approximation
%   g1m start point of summation
%   g2m end point of summation

syms y
% Legendre basis
bpy = cell(1,My+1);
for k = 1:My+1
  bpy{k} = legendreP(k-1,y);
end

% initialize and compute coefficient matrices and vectors
BV = zeros(My+1); DV = zeros(My+1,1); Aeq = zeros(1,My+1);
g1 = 2*(m1-1)/m-1; g2 = 2*m2/m-1;
for k = 0:My
  for l = 0:My
    BV(k+1,l+1) = int(bpy{k+1}*bpy{l+1},y,g1,g2);
  end
  Aeq(k+1) = int(bpy{k+1},y,g1,g2);
  DVtmp = 0; 
  for l = m1:m2
    DVtmp = DVtmp + dm(l)*int(bpy{k+1},y,2*(l-1)/m-1,2*l/m-1);
  end
  DV(k+1) = DVtmp;
end

% solve unknown coefficients in the polynoial fit
beq = 0;
for k = 1:m
  beq = beq + dm(k);
end
beq = 2*beq/m;
opts = optimset('Display','none');
CF = quadprog(BV,-DV,[],[],Aeq,beq,[],[],[],opts);
% compute the polynomial
fpol = 0;
for k = 0:My
  fpol = fpol + CF(k+1)*bpy{k+1};
end
% compute L2 residual
rd = 0;
for l = m1:m2
  rd = rd + int((dm(l) - fpol)^2,y,2*(l-1)/m-1,2*l/m-1);
end
rd = double(rd);
end