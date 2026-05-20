function [G1s, G2s, rd] = GammaSolver(lam,mu,W,TH,Q,R,F,A,C,N,Ny)
%GAMMASOLVER Solves a 1.5D ODE for the output injection gains

syms x y
% Taylor series approximations for PDE parameters
Ls = taylor(lam,[x y],'order',N+1);
Ms = taylor(mu,[x y],'order',N+1);
Ws = taylor(W,[x y],'order',N+1);
THs = taylor(TH,[x y],'order',N+1);
Qs = taylor(Q,y,'order',N+1);
Rs = taylor(R,y,'order',N+1);
Fs = taylor(F,y,'order',N+1);

% initialize power series for 
NG = (2*N*Ny + 2*N + Ny - Ny^2 + 2)/2; % number of terms
KC = sym('k', [1 2*NG]); % symbols for coefficients
rind = 0;
% Kbar Taylor (2D)
Uk = 0;
Vk = 0;
for tm=0:Ny
  for tn=0:(N-tm)
    rind = rind + 1;
    Uk = Uk + KC(rind)*x^tn*y^tm;
    Vk = Vk + KC(NG+rind)*x^tn*y^tm;
  end
end

[V, E] = eig(A);
G1s = zeros(size(C)); G2s = G1s; % stores the solution
rd = zeros(size(C)); % residuals for least squares fits
for k = 1:numel(C)
  % solve over the eigenvalues and vectors of A
  sk = E(k,k); vk = V(:,k);
  % construct equations and boundary conditions with power series
  E1 = sk*Uk + Ls*diff(Uk,x) - Ws*Vk;
  E2 = sk*Vk - Ms*diff(Vk,x) - THs*Uk;
  BC1 = subs(Uk,x,-1) - Qs*subs(Vk,x,-1);
  BC2 = subs(Vk,x,1) - Rs*subs(Uk,x,1) - Fs*C*vk;
  % extract coefficients for different powers of spatial variables
  E1C = coeffs(E1,[x y]);
  E2C = coeffs(E2,[x y]);
  BC1C = coeffs(BC1,y);
  BC2C = coeffs(BC2,y);
  % construct set linear equations for the coefficients based on the above
  neqs = [numel(E1C), numel(E2C), numel(BC1C), numel(BC2C)];
  AA = zeros(sum(neqs), 2*NG); % initialize A...
  B = zeros(sum(neqs),1); % ...and b
  EQS = {E1C,E2C,BC1C,BC2C}; % cell array for the equation data
  % disp(['Parsing data: ',num2str(sum(neqs)),' equations for ',...
  %   num2str(2*NG),' unknonws']);
  rind = 0; % row index
  for ii = 1:4
   DAT = EQS{ii};
    for m=1:neqs(ii)
     rind = rind + 1;
      [cc, kc] = coeffs(DAT(m),KC);
      for mk=1:numel(kc) % check which K's are present in the coefficients
        str = char(kc(mk)); 
       sl = numel(str);
       if sl > 1 % look up the index and insert to AA
         AA(rind,str2double(str(2:sl))) = cc(mk);
       else % if no index, it is a contant; insert to B
         B(rind) = -cc(mk);
       end
      end
    end
  end
  % solve equations and compute residual
  Kc = AA\B;
  rd(k) = norm(AA*Kc-B);

  % insert the obtained values into the power series
  Uks = subs(Uk,KC(1:NG),Kc(1:NG)');
  Vks = subs(Vk,KC(NG+1:2*NG),Kc(NG+1:2*NG)');
  % sum over eigenvectors of AA for full solution
  G1s = G1s + Uks*vk.';
  G2s = G2s + Vks*vk.';
  % G1s and G2s may have small imaginary parts due to round-off errors,
  % can be elinated by taking real part at the end
end
end