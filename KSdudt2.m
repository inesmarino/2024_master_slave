function d = KSdudt2(W0,Ix,a0,alpha,beta,gamma,N)
%
% function d = KSdudt2(W0,Ix,a0,alpha,beta,gamma,N)
%
% W0 : powers of w0, the fundamental freq.
% Ix : powers of the coefficient indices, from -N to N
% a0 : previous value of the FS coefficients (2N+1 x 1)
% alpha, beta, gamma : model parameters
% N : order of the FS approximation
%
% d : time derivative of the FS coefficients from -N to 0
%

idn = 1:N+1;

pa0 = [ zeros([N 1]); a0; zeros([N 1])];   % padded a(k)
py = conv(pa0,pa0);                        % padded y(m)
y = py(3*N+1:5*N+1);                       % y(m), -N≤m≤N

% time derivatives
M = [ W0(2).*Ix(idn,2).*a0(idn) 1i*W0(3).*Ix(idn,3).*a0(idn) -W0(4).*Ix(idn,4).*a0(idn) ];
d = -0.5*1i*W0(1).*Ix(idn,1).*y(idn) + M*[alpha beta gamma]';