syms x y
% m(2x2) parameters as symbolic expressions of x
% instead of [0,1], the spatial domain is [-1,1] due to definition of
% Legendre polynomials, the respective change of variables is x'=(x+1)/2
m = 10;
cw = [-3,-3,-2,-2,1,0,1,2,3,6]'/20;
cth = [12,11,14,13,15,16,18,19,18,20]'/20;
q = [1,2,4,5,4,7,7,7,8,10]'/10;
r = [9,8,6,4,4,3,2,2,-1,0]'/10;
lamm = @(i) 2;
mum = @(i) 2;
Wm = @(i) cw(i)*(x+3)/2.*((x+1)/2); 
THm = @(i) cth(i)*(x+1)/2;
Qm = @(i) q(i);
Rm = @(i) r(i);
Fm = @(i) 1;
gm = @(i) 1;
g1m = 1; g2m = m;

% spectral approximation of plant
Nxm = 15; % spectral approximation order in x for m(2x2) system
% spectral modal approximation for the m(2x2) system
[MMm, KKm, BBm, CCm] = smam(m,Nxm,lamm,mum,Wm,THm,Qm,Rm,Fm,gm,g1m,g2m);

% ODE part with initial condition
[A, C, x0] = waveform('Q');
C = 1e3*C; % scaling of harmonics
n = numel(C);

for ex = 1:2 % loop over approximation orders 1 and 3 in y
% orders for polynomial approximnation of parameters
M = 3; 
switch ex
  case 1
    My = 1;
  case 2
    My = 3;
end
% construction of continuum parameters
% the y domain is also [-1,1] instead of [0,1], same as for x
lam = 2;
mu = 2;
[W, rdW] = L2polyfit2(Wm,m,M,My,false);
[TH, rdTH] = L2polyfit2(THm,m,M,My,false);
[Q, rdQ] = L2polyfit(Qm,m,My);
[R, rdR] = L2polyfit(Rm,m,My);
F = 1;
g = 1;
g1 = 2*(g1m-1)/m-1; g2 = 2*g2m/m-1;

% compute observer gains
N = 40; % order for power series approximation
% reduced order in y
switch ex
  case 1
    Ny = 4;
  case 2
    Ny = 8;
end
[G1, G2, rd] = GammaSolver(lam,mu,W,TH,Q,R,F,A,C,N,Ny);
% take real parts of G1s,G2s in case complex residuals due to round-off
GC = double(real(int(g*subs(G2,x,-1),g1,g2))); % aux vector
L = -lqr(A',GC',1e4*eye(n),1)'; % output injection gain by LQR

% spectral approximation of continuum observer
Nx = 15; Ny = My+1; % spectral approximation orders in x and y
% spectral modal approximation for the observer
[MMc, KKc, BBc, CCc, PPc] = smac(Nx,Ny,lam,mu,W,TH,Q,R,F,g,g1,g2,L,G1,G2);

% closed-loop system MM*\dot{ae} = KK*ae
MM = blkdiag(eye(n),MMm,eye(n),MMc);
KK = [A, zeros(n,2*m*Nxm+n+2*Nx*Ny); ...
  BBm*C, KKm, zeros(2*m*Nxm,n+2*Nx*Ny); ...
  zeros(n), -L*CCm, A, L*CCc; ...
  zeros(2*Nx*Ny,n), -PPc*CCm, BBc*C, KKc+PPc*CCc];

% initial conditions and conversion to basis coefficiens based on the
% method from https://github.com/Kristian-MJA/Satmodel (c) Kristian Asti 2020
% plant
um0 = @(x,i) cos(pi/2*x);
vm0 = @(x,i) cos(pi/2*x);
Ntm = 2*Nxm;
xkm = linspace(-1,1,Ntm).';
phimat = zeros(Ntm,Nxm);
for k = 1:Nxm
    phimat(:,k) = legendreP(k-1,xkm)';
end
% initial condition as basis function coefficients
uma0 = zeros(m*Nxm,1); vma0 = uma0;
for k = 1:m
  uma0((k-1)*Nxm+1:k*Nxm) = phimat\um0(xkm,k);
  vma0((k-1)*Nxm+1:k*Nxm) = phimat\vm0(xkm,k);
end
% same for observer, although the PDE part is initalized to zero
Ntx = 2*Nx;
Nty = 2*Ny;
xkc = linspace(-1,1,Ntx).';
ykc = linspace(-1,1,Nty).';
xo0 = ones(size(x0));
ucoa0 = zeros(Nx*Ny,1);
vcoa0 = zeros(Nx*Ny,1);
% closed-loop initial condition
a0 = [x0; uma0; vma0; xo0; ucoa0; vcoa0];

% simulate closed-loop system
tspan = [0 15];
options = odeset('mass', MM,'jacobian', KK);
sol = ode45(@(t,x) KK*x, tspan, a0, options);

% evaluate solution on desired time grid
tt = linspace(tspan(1), tspan(end), 513);
asol = deval(sol,tt);
xsol = asol(1:n,:); % ODE solution
xosol = asol(n+2*m*Nxm+1:n+2*m*Nxm+n,:); % ODE estimate
% % for plotting selected m(2x2) sollution component (not used)
% xg = linspace(-1,1, 129)';
% phimat = zeros(numel(xg),Nxm);
% for k = 1:Nxm
%     phimat(:,k) = legendreP(k-1,xg);
% end
% pj = 1; % subsystem index
% uj = phimat*asol(n+(pj-1)*Nxm+1:n+pj*Nxm,:);
% vj = phimat*asol(n+(pj+m-1)*Nxm+1:n+(pj+m)*Nxm,:);
% surf(tt,(xg+1)/2,uj)

switch ex
  case 1
    Yerr1 = C*xsol - C*xosol;
    figure(1)
  case 2
    Yerr3 = C*xsol - C*xosol;
    figure(2)
end
subplot(211)
plot(tt,C*xsol,tt,C*xosol,'linewidth',2)
set(gca,'tickdir','out','xtick',0:5:20,'ytick',-1:1,'fontsize',13, ...
  'position',[0.13 0.6 0.84 .36],'xticklabel',[])
ylim([-.1 1.3])
legend('$\mathbf{CX}(t)$', '$\mathbf{C}\hat{\mathbf{X}}(t)$',...
  'interpreter','latex','box','off','location','northeast','numcolumns',2)
subplot(212)
plot(tt,C*xsol - C*xosol,'linewidth',2)
set(gca,'tickdir','out','xtick',0:5:20,'ytick',-.5:.5:.5,'fontsize',13,...
  'position',[0.13 0.18 0.84 .36])
ylim([-.9 .1])
ylabel('$\mathbf{C}\mathbf{X}(t)-\mathbf{C}\hat{\mathbf{X}}(t)$',...
  'interpreter','latex','fontsize',13)
xlabel('$t$','interpreter','latex')
end
figure(3)
semilogy(tt,abs(Yerr1),tt,abs(Yerr3),'linewidth',2)
set(gca,'tickdir','out','xtick',0:5:20,'fontsize',13, ...
  'position',[0.13 0.25 0.84 .65])
xlabel('$t$','interpreter','latex')
legend('$M_y = 1$', '$M_y = 3$','interpreter','latex',...
   'box','off','location','northeast','numcolumns',2)
ylim([9e-4, 1e0])