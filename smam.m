function [MM, KK, BB, CC] = smam(m,Nx,lam,mu,W,TH,Q,R,F,g,m1,m2)
%SMAM Construct spectral modal approximation MM*\dot{a} = KK*a + BB*U with 
% output y = CC*a for m(2x2) system
%   m:  number of 2x2 systems
%   Nx: approximation order in x
%   lam, mum W, TH, Q, R, F,g, m1, m2: 2x2 system parameters

syms x
% Legendre polynomials up to order Nx-1
bx = cell(1,Nx);
for k = 1:Nx
  bx{k} = legendreP(k-1,x);
end

% initializations of auximilary matrices
MMu = zeros(m*Nx); Kuu = zeros(m*Nx); Kvv = zeros(m*Nx);
Kvu = zeros(m*Nx); Kuv = zeros(m*Nx);
BBv = zeros(m*Nx,1); CCv = zeros(1,m*Nx);

% fill auxiliary matrices
for j = 1:m
  for k = 1:Nx
    MMu((j-1)*Nx+k,(j-1)*Nx+k) = int(bx{k}*bx{k},x,-1,1);
    BBv((j-1)*Nx+k) = F(j)*subs(mu(j),x,1)*subs(bx{k},x,1);
    if j >= m1 && j <= m2
      CCv((j-1)*Nx+k) = g(j)*subs(bx{k},x,-1)/m;
    end
    for l = 1:Nx
      Kuu((j-1)*Nx+k,(j-1)*Nx+l) = int(bx{l}*bx{k}*diff(lam(j),x) + ...
        bx{l}*diff(bx{k},x)*lam(j),x,-1,1) - ...
        subs(lam(j),x,1)*subs(bx{k},x,1)*subs(bx{l},x,1);
      Kvu((j-1)*Nx+k,(j-1)*Nx+l) = int(bx{l}*bx{k}*W(j),x,-1,1) + ...
        Q(j)*subs(lam(j),x,-1)*subs(bx{k},x,-1)*subs(bx{l},x,-1);
      Kuv((j-1)*Nx+k,(j-1)*Nx+l) = int(bx{l}*bx{k}*TH(j),x,-1,1) + ...
        R(j)*subs(mu(j),x,1)*subs(bx{k},x,1)*subs(bx{l},x,1);
      Kvv((j-1)*Nx+k,(j-1)*Nx+l) = -int(bx{l}*bx{k}*diff(mu(j),x) + ...
        bx{l}*diff(bx{k},x)*mu(j),x,-1,1) - ...
        subs(mu(j),x,-1)*subs(bx{k},x,-1)*subs(bx{l},x,-1);
    end
  end
end

% full matrices for spectral modal approxiamtion
MM = [MMu, zeros(m*Nx); zeros(m*Nx), MMu];
KK = [Kuu, Kvu; Kuv, Kvv];
BB = [zeros(m*Nx,1); BBv];
CC = [zeros(1,m*Nx), CCv];
end