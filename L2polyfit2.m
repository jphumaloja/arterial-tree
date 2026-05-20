function [fpol, rd] = L2polyfit2(fm,m,M,My,dx)
%L2POLYFIT2 Finds the L2 least-squares polynomial to match given functions
%   fm  functions to fit the polynomial to
%   m   number of functions in fm
%   M   total order of the polynomial approximation
%   My  order of the polynomial approximation in y (My <= M)
%   dx  fit with respect to fm/dx as well (true or false)

syms x y
kind = @(i,j) i + M*j - (j-1)*(j-2)/2 + 2; % indexing function
Mk = M*(My+1) - My*(My-1)/2 + 1; % number of unknowns

% Legendre basis
bpx = cell(1,M+1); bpy = cell(1,My+1);
for k = 1:M+1
  bpx{k} = legendreP(k-1,x);
end
for k = 1:My+1
  bpy{k} = legendreP(k-1,y);
end

% initialize and compute coefficient matrices and vectors
BV = zeros(Mk); FV = zeros(Mk,1);
if dx
  BD = zeros(Mk); FD = zeros(Mk,1);
end
for j = 0:My
  wpy2 = int(bpy{j+1}^2,y,-1,1);
  for i = 0:M-j
    k = kind(i,j);
    BV(k,k) = int(bpx{i+1}^2,x,-1,1)*wpy2;
    if dx
      BD(k,k) = int(diff(bpx{i+1},x)^2,x,-1,1)*wpy2;
     if i > 2
       BD(k,k-2) = int(diff(bpx{i+1},x)*diff(bpx{i-1},x),x,-1,1)*wpy2;
      end
      if i < (M-j-1)
        BD(k,k+2) = int(diff(bpx{i+1},x)*diff(bpx{i+3},x),x,-1,1)*wpy2;
      end
      FDtmp = 0;
    end
    FVtmp = 0; 
    for l = 1:m
      wpym = int(bpy{j+1},y,2*(l-1)/m-1,2*l/m-1);
      FVtmp = FVtmp + int(fm(l)*bpx{i+1},x,-1,1)*wpym;
      if dx
        FDtmp = FDtmp + int(diff(fm(l),x)*diff(bpx{i+1},x),x,-1,1)*wpym;
      end
    end
    FV(k) = FVtmp;
    if dx
      FD(k) = FDtmp;
    end
  end
end
% solve for the unknown coefficients in the polynomial fit
if ~dx
  CF = BV\FV;
else
  CF = (BV+BD)\(FV+FD);
end
% compute the polynomial
fpol = 0;
for j = 0:My
  for i = 0:M-j
    k = kind(i,j);
    fpol = fpol + CF(k)*bpx{i+1}*bpy{j+1};
  end
end
% compute L2 residual
rd = 0;
for l = 1:m
  rd = rd + int(int((fm(l) - fpol)^2,x,-1,1),y,2*(l-1)/m-1,2*l/m-1);
  if dx
    rd = rd + int(int(diff(fm(l) - fpol,x)^2,x,-1,1),y,2*(l-1)/m-1,2*l/m-1);
  end
end
rd = double(rd);
end