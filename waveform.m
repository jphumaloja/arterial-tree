function [A, C, x0] = waveform(type)
%WAVEFORM Generates an ODE that outputs a desired flow/pressure waveform
%   type: 'P' or 'Q' for pressure or flow waveform, respectively.

bfreq = 2*pi; % base frequency of Fourier series
d1 = zeros(1,2*4);
d1(2:2:2*4) = (1:4)*bfreq;
A = diag(d1,1) - diag(d1,-1); % base signal generator
x0 = [1; repmat([0; 1],4,1)]; % initial condition

% generate C based on Fourier coefficients, add decay to pressure waveform
switch type
  case 'P'
    dr = log(17/17.5); % decay rate for pressure waveform
    A = A + diag(dr*ones(1,9)); % add decay to A
    % fit Fourier series to pressure data and insert coefficients to C
    xdat = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
    ydat = [12.5 11.25 10 13.75 17.5 15 15.5 14.75 14 13.25 12.5];
    F = fit(xdat',ydat',"Fourier4"); % Fourier fit of order 4
    C = 1e3*[F.a0 F.b1 F.a1 F.b2 F.a2 F.b3 F.a3 F.b4 F.a4];
  case 'Q'
    % insert flow flow waveform coefficients to C
    C = [0.86393, 1.3368, -0.88455, -1.228, -0.52515, 0.22459, 0.86471, ...
      0.22693, -0.26395]*1e-4;
end
end