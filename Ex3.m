syms x y z u v
% physical parameters
h = 0.5e-3; E = 4e5; b = 4/3;
rho = 1060; nu = 0.0035/rho; Kr = 8*pi*nu; beta = h*E*sqrt(pi)*b;
Ks = 0; RT = 1.33e8;
m = 10;
r0 = 5e-3 + (1:m)/2*1e-4; A0 = pi*r0.^2; 
kap = 2^9*Kr*(beta/rho./A0).^2;
% steady states (at rest)
Qin = 0; As = A0; Qs = 0; Vs = Qs./As;
% Riemann steady states
vs = Vs + 2*sqrt(2*beta./(rho*A0)).*As.^(1/4);
us = Vs - 2*sqrt(2*beta./(rho*A0)).*As.^(1/4);
% aux coefficients
d0 = (rho*A0/beta).^2/2^11;
d1 = RT/4^3*(rho*A0/beta).^2;
% aux functions
Psu = rho/2^4.*(us-vs);
Psv = rho/2^4.*(vs-us);
G = rho.*(u-v).^2 - 32*beta./sqrt(A0) - d1.*(u-v).^4.*(u+v);
Gu = diff(G,u); Gv = diff(G,v);
Gsu = zeros(1,m); Gsv = zeros(1,m);
for k = 1:m
  Gsu(k) = subs(Gu(k),[u v],[subs(us(k),x,0), subs(vs(k),x,0)]);
  Gsv(k) = subs(Gv(k),[u v],[subs(us(k),x,0), subs(vs(k),x,0)]);
end

% parameters of the linearized 2x2 system
lamd = -(5*us + 3*vs)/8;
mud = (3*us + 5*vs)/8;
S11 = 5/8*diff(us,x) - kap.*(3*us+5*vs)./(us-vs).^5;
S12 = 3/8*diff(us,x) + kap.*(5*us+3*vs)./(us-vs).^5;
S21 = 3/8*diff(vs,x) - kap.*(3*us+5*vs)./(us-vs).^5;
S22 = 5/8*diff(vs,x) + kap.*(5*us+3*vs)./(us-vs).^5;
rod = -Psu./Psv; fod = 1./Psv; qd = -Gsv./Gsu;
% parameters corresponding to change of variables
wd = S12.*exp(int(subs(S22./mud + S11./lamd,x,z),z,0,x));
thd = S21.*exp(-int(subs(S22./mud + S11./lamd,x,z),z,0,x));
rd = rod.*exp(subs(S22./mud + S11./lamd,x,1));
fd = fod.*exp(subs(S22./mud,x,1));
gd = ones(1,m); g1m = 1; g2m = m;
% conversion to functions of i and x \in [-1,1] from x \in [0,1]; see Ex1.m
lamm = @(i) 2*lamd(i);
mum = @(i) 2*mud(i);
Wm = @(i) subs(wd(i),x,(x+1)/2);
THm = @(i) subs(thd(i),x,(x+1)/2);
Qm = @(i) qd(i);
Rm = @(i) rd(i);
Fm = @(i) fd(i);
gm = @(i) gd(i);

% spectral approximation of plant
Nxm = 15; % spectral approximation order in x for m(2x2) system
% spectral modal approximation for the m(2x2) system
[MMm, KKm, BBm, CCm] = smam(m,Nxm,lamm,mum,Wm,THm,Qm,Rm,Fm,gm,g1m,g2m);

% ODE part with initial condition
[A, C, x0] = waveform('P');
n = numel(C);

% orders for polynomial approximnation of parameters
M = 3; My = 1;
[W, rdW] = L2polyfit2(Wm,m,M,My,false);
[TH, rdTH] = L2polyfit2(THm,m,M,My,false);
[Q, rdQ] = L2polyfit(Qm,m,My);
[R, rdR] = L2polyfit(Rm,m,My);
[F, rdF] = L2polyfit(Fm,m,My);
[lam, rdlam] = L2polyfit(lamm,m,My);
[mu, rdmu] = L2polyfit(mum,m,My);
[g, rdg] = L2polyfitg(gm,m,My,g1m,g2m);
g1 = 2*(g1m-1)/m-1; g2 = 2*g2m/m-1;

% compute observer gains
N = 40; % order for power series approximation
Ny = 4; % reduced order in y
[G1, G2, resd] = GammaSolver(lam,mu,W,TH,Q,R,F,A,C,N,Ny);
% take real parts of G1s,G2s in case complex residuals due to round-off
GC = double(real(int(g*subs(G2,x,-1),g1,g2))); % aux vector
L = -lqr(A',GC',100*eye(n),1)'; % output injection gain L

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
um0 = @(x,i) -1*exp(-double(S11(i))/lamm(i)*(x+1)/2);
vm0 = @(x,i) 4*exp(double(S22(i))/mum(i)*(x+1)/2);
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
xo0 = ones(size(x0));
ucoa0 = zeros(Nx*Ny,1);
vcoa0 = zeros(Nx*Ny,1);
% closed-loop initial condition
a0 = [x0; uma0; vma0; xo0; ucoa0; vcoa0];

% simulate closed-loop system
tspan = [0 7];
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

subplot(211)
plot(tt,C*xsol,tt,C*xosol,'linewidth',2)
set(gca,'tickdir','out','fontsize',13, ...
  'position',[0.13 0.59 0.84 .34],'xticklabel',[])
ylim([0 2.5]*1e4)
legend('$\mathbf{CX}(t)$', '$\mathbf{C}\hat{\mathbf{X}}(t)$',...
  'interpreter','latex','box','off','location','southeast','numcolumns',2)
subplot(212)
plot(tt,C*xsol - C*xosol,'linewidth',2)
set(gca,'tickdir','out','fontsize',13,...
  'position',[0.13 0.16 0.84 .34])
ylim([-1.2 1]*1e4)
ylabel('$\mathbf{C}\mathbf{X}(t)-\mathbf{C}\hat{\mathbf{X}}(t)$',...
  'interpreter','latex','fontsize',13)
xlabel('$t$','interpreter','latex')