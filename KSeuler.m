function [u,a,ok] = KSeuler(alpha,beta,gamma,tgrd,N,a0,s_f,x)
%
%  function [u,a,ok] = KSeuler(alpha,beta,gamma,tgrd,N,u0,s_f,x)
%
% alpha, beta, gamma : K-S equation parameters
% tgrd : time grid
% N : Fourier coefficients from -N to N
% a0 : initial condition of the FS coefficients
% s_f : space is [0, s_f), periodic
% x : space grid
%
% u : signal
% a : Fourier coefficients from -N to N
% ok : 1 if the simulation completes without blowing up, 0 otherwise
%

ok = 1;

% length of the time grid
T = length(tgrd);

% FS coefficients
a = zeros([2*N+1 T]);
a(1:2*N+1,1) = a0;

% Fundamental frequency & powers
w0 = 2*pi/s_f;
w02 = w0^2;
w03 = w0^3;
w04 = w0^4;
W0 = [w0 w02 w03 w04]';

% Freq. index & powers
idx = (-N:N)';
idx2 = ((-N:N).^2)';
idx3 = ((-N:N).^3)';
idx4 = ((-N:N).^4)';
Ix = [idx idx2 idx3 idx4];

% Imaginary unit
j = sqrt(-1);

% iFT matrix
FM = exp( j*w0*x*(-N:N) ); 
u = zeros([length(x) T]);
u(:,1) = FM*a0;

% time loop
for t = 2:T

    % time step
    ts = tgrd(t)-tgrd(t-1);

    % Euler
    idn = 1:N+1;
    dudt = KSdudt2(W0,Ix,a(:,t-1),alpha,beta,gamma,N);
    a(idn,t) = a(idn,t-1) + ts.*dudt;
    a(N+2:end,t) = conj(a(N:-1:1,t));

    % reconstruction of the signal
    u(:,t) = FM*a(:,t);

    % check for numerical stability
    if max(abs(a(:,t)))>1e+2
        ok = 0;
        return;
    end %if

end %t
