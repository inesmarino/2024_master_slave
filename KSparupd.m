function par = KSparupd(W0,Ix,b,par0,N,ts,ae,D,mu)
%
% function par = KSparupd(W0,Ix,b,par0,N,ts,ae,D,mu)
%
% W0 : powers of w0, the fundamental freq.
% Ix : powers of the coefficient indices, from -N to N
% b : available value of the FS coefficients (2N+1 x 1) at time t-h
% par0 : model parameters at the previous step
% N : order of the FS approximation
% ts : time step
% no -> Phi : iFT matrix for the observations available
% no -> u : observed signal 
% ae : estimated FS coeffs, this is ae(:,t-1:t)
% D : coupling strength
% mu : adaptation rate
%
% par : updated parameter vector
%



pa0 = [ zeros([N 1]); b; zeros([N 1])];   % padded FS coeffs a(k)
py = conv(pa0,pa0);                       % padded y(m)
y = py(3*N+1:5*N+1);                      % y(m), -N≤m≤N
fb = -0.5*1i*W0(1).*Ix(:,1).*y;

% linear part on the parameters
M = [ W0(2).*Ix(:,2).*b 1i*W0(3).*Ix(:,3).*b -W0(4).*Ix(:,4).*b ];

% parameter update
par = ( eye(3) + mu*(ts^3)*(M'*M) ) \ ( par0 + mu*(ts^2)*M'*( ae(:,2) - (b + ts*fb + ts*D*(ae(:,1)-b) ) ) );