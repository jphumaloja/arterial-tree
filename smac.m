function [MM,KK,BB,CC,PP] = smac(Nx,Ny,lam,mu,W,TH,Q,R,F,g,g1,g2,L,G1,G2)
%SMAC Construct spectral modal approximation MM*\dot{a} = KK*a + BB*U for a
%continuum observer with output y = CC*a and output injection gain L
%   Nx: approximation order in x
%   Ny: approximation order in y
%   lam, mum W, TH, Q, R, F,g, m1, m2: continuum system parameters
%   L, G1, G2: for computing output injection gains

syms x y
% Legendre polynomials up to orders Nx-1 and Ny-1 in x and y, respectively
bx = cell(1,Nx); by = cell(1,Ny);
for k = 1:Nx
  bx{k} = legendreP(k-1,x);
end
for k = 1:Ny
  by{k} = legendreP(k-1,y);
end

% initializations of auximilary matrices
MMu = zeros(Nx*Ny); Kuu = zeros(Nx*Ny); Kvv = zeros(Nx*Ny);
Kvu = zeros(Nx*Ny); Kuv = zeros(Nx*Ny);
BBv = zeros(Nx*Ny,1); CCv = zeros(1,Nx*Ny);
PP1 = zeros(Nx*Ny,1); PP2 = zeros(Nx*Ny,1);

% fill auxiliary matrices
for k = 1:(Nx*Ny)
  [xk, yk] = ixy(k,Nx);
  MMu(k,k) = int(int(bx{xk}^2*by{yk}^2,x,-1,1),y,-1,1);
  BBv(k) = int(F*subs(mu,x,1)*subs(bx{xk},x,1)*by{yk},y,-1,1);
  % division by 2 comes from extending the domain from [0,1] to [-1,1]
  CCv(k) = int(g*subs(bx{xk},x,-1)*by{yk},y,g1,g2)/2;
  % take real parts of G1s,G2s in case complex residuals due to round-off
  PP1(k) = real(int(int(G1*L*bx{xk}*by{yk},x,-1,1),y,-1,1));
  PP2(k) = real(int(int(G2*L*bx{xk}*by{yk},x,-1,1),y,-1,1));
  for l = 1:(Nx*Ny)
    [xl, yl] = ixy(l,Nx);
    Kuu(k,l) = int(int(bx{xl}*by{yl}*diff(lam,x)*bx{xk}*by{yk} + ...
      bx{xl}*by{yl}*lam*diff(bx{xk},x)*by{yk},x,-1,1) - ...
      subs(lam,x,1)*subs(bx{xk},x,1)*by{yk}*subs(bx{xl},x,1)*by{yl},y,-1,1);
    Kvu(k,l) = int(int(bx{xl}*by{yl}*W*bx{xk}*by{yk},x,-1,1) + ...
      Q*subs(lam,x,-1)*subs(bx{xk},x,-1)*by{yk}*subs(bx{xl},x,-1)*by{yl},y,-1,1);
    Kuv(k,l) = int(int(bx{xl}*by{yl}*TH*bx{xk}*by{yk},x,-1,1) + ...
      R*subs(mu,x,1)*subs(bx{xk},x,1)*by{yk}*subs(bx{xl},x,1)*by{yl},y,-1,1);
    Kvv(k,l) = -int(int(bx{xl}*by{yl}*diff(mu,x)*bx{xk}*by{yk} + ...
      bx{xl}*by{yl}*mu*diff(bx{xk},x)*by{yk},x,-1,1) + ...
      subs(mu,x,-1)*subs(bx{xk},x,-1)*by{yk}*subs(bx{xl},x,-1)*by{yl},y,-1,1);
  end
end

% full matrices for spectral modal approxiamtion
MM = [MMu, zeros(Nx*Ny); zeros(Nx*Ny), MMu];
KK = [Kuu, Kvu; Kuv, Kvv];
BB = [zeros(Nx*Ny,1); BBv];
CC = [zeros(1,Nx*Ny), CCv];
PP = [PP1; PP2];
end

% auxiliary indexing function
function [ix, iy] = ixy(k,Nx)
ix = mod(k-1,Nx)+1;
iy = floor((k-1)/Nx)+1;
end