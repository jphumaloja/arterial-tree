function [fpol, rd] = L2polyfit(dm,m,My)
%L2POLYFIT Finds the L2 least-squares polynomial to match given data
%   dm  data to fit the polynomial to
%   m   number of functions in dm
%   My  order of the polynomial approximation

syms y
% Legendre basis
bpy = cell(1,My+1);
for k = 1:My+1
  bpy{k} = legendreP(k-1,y);
end

% initialize and compute coefficient matrices and vectors
BV = zeros(My+1); DV = zeros(My+1,1);

for k = 0:My
  BV(k+1,k+1) = int(bpy{k+1}^2,y,-1,1);
  DVtmp = 0; 
  for l = 1:m
    DVtmp = DVtmp + dm(l)*int(bpy{k+1},y,2*(l-1)/m-1,2*l/m-1);
  end
  DV(k+1) = DVtmp;
end

% solve for the unknown coefficients in the polynomial fit
CF = BV\DV;

% compute the polynomial
fpol = 0;
for k = 0:My
  fpol = fpol + CF(k+1)*bpy{k+1};
end
% compute L2 residual
rd = 0;
for l = 1:m
  rd = rd + int((dm(l) - fpol)^2,y,2*(l-1)/m-1,2*l/m-1);
end
rd = double(rd);
end