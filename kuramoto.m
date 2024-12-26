%% Simulation parameters
swa = [     0   %Euler slave, coupling only
            1   %parameter update
    ];
NR = 2;     % no. of repetitions (noise is independently generated at each run)

%% Model parameters
alpha = 1.15;
beta = -0.05;
gamma = 0.98;
t_f = 100;                          % final continuous time
T = 2e+4;                           % no. of nodes in the time grid (excluding 0)
s_f = 120;                          % final continuous space
tgrid = 0:(t_f/T):t_f-(t_f/T);      % time grid
N = 64;                             % Fourier coefficients (from -N to N) in the master model
L = N*ones([1 NR]);                 % Fourier coefficients in the slave model: we may have different a value for each repetition
mu = 2e+2;                          % adaptation rate in the slave model
SNR = 12;                           % signal to noise ratio in the observations (SNR≥100 => zero noise)
xstep = 1/4;                        % spatial resolution for plots
xgrid = (0:xstep:s_f)';             % observation grid
x_obs_gap = 2;                      % spatial gap between observations (in 'xstep' units) 
D = 0.5;                            % coupling parameter


%% Initialisation
u = cell([NR 1]);
v = cell([NR 1]);
w = cell([NR 1]);

a = cell([NR 1]);
a0 = cell([NR 1]);

b = cell([NR 1]);
b0 = cell([NR 1]);

c = cell([NR 1]);
c0 = cell([NR 1]);

x = cell([NR 1]);

rt_e = zeros([NR 1]);
rt_ep = rt_e;
rt_cp = rt_e;

param = cell([NR 1]);

alpha0 = zeros([NR 1]);
beta0 = zeros([NR 1]);
gamma0 = zeros([NR 1]);


%% Simulations
for nr = 1:NR

    % initialisation of the Fourier coefficients in the master model
    
    % initialisation master model with N=32
    % load init_a0_N32.mat a0_N32;
    % a0{nr} = a0_N32;

    % initialisation master model with N=64
    load init_a0_N64.mat a0_N64;
    a0{nr} = a0_N64;

    % master model: spectral (Fourier series) decomposition + Euler scheme
    t0 = tic;
    [u{nr},a{nr},ok] = KSeuler(alpha,beta,gamma,tgrid,N,a0{nr},s_f,xgrid);
    rt_e(nr) = toc(t0);
    fprintf(1,'K-S, Fourier decomposition, master, Euler, time=%6.3f s\n', rt_e(nr));

    % slave model

    if swa(1) % coupling only, no parameter estimation
        
        t0 = tic;

        % initialisation
        b0{nr} = zeros([2*L(nr)+1 1]);
        b0{nr}([L(nr) L(nr)+2]) = 0.5;

        % observation grid
        ygrid = xgrid(1:x_obs_gap:end);
        Uy = u{nr}(1:x_obs_gap:end,1:T);
        J = size(Uy,1);

        % observation noise
        Pu = mean(sum(real(Uy).^2))/J;
        if SNR<100
            s2y = Pu*10^(-SNR/10);
        else 
            s2y = 0;
        end
        noise = sqrt(s2y).*randn([J T]);

        % approx. coeffs from observations
        w0 = 2*pi/s_f;
        FM = exp( 1i*w0*ygrid*(-L(nr):L(nr)) );
        ae = (FM'*FM)\FM'*(Uy+noise);
        ae(N+1,:) = real(ae(N+1,:));

        % Euler with coupling, fixed parameters
        [v{nr},b{nr},ok] = KScoupled_euler(alpha,beta,gamma,tgrid,L(nr),b0{nr},s_f,xgrid,ae,D);

        % timing
        rt_cp(nr) = toc(t0);
        fprintf(1,'K-S, master-slave, identical parameters, L=%d, nr=%d/%d, time=%6.3f s\n', L(nr), nr, NR, rt_cp(nr));

    end %if


    if swa(2)   % parameter update
        t0 = tic;

        % parameters
        alpha0(nr) = 0;
        beta0(nr) = 0;
        gamma0(nr) = 0;

        % initialisation
        c0{nr} = zeros([2*L(nr)+1 1]);
        c0{nr}([L(nr) L(nr)+2]) = 0.5;

        % observation grid
        ygrid = xgrid(1:x_obs_gap:end);
        Uy = u{nr}(1:x_obs_gap:end,1:T);
        J = size(Uy,1);

        % observation noise
        Pu = mean(sum(real(Uy).^2))/J;
        if SNR<100
            s2y = Pu*10^(-SNR/10);
        else 
            s2y = 0;
        end
        noise = sqrt(s2y).*randn([J T]);

        % approx. coeffs from observations
        w0 = 2*pi/s_f;
        FM = exp( 1i*w0*ygrid*(-L(nr):L(nr)) );
        ae = (FM'*FM)\FM'*(Uy+noise);
        ae(L(nr)+1,:) = real(ae(L(nr)+1,:));

        % Euler with coupling and parameter update
        [w{nr},c{nr},param{nr},ok] = KScoupled(alpha0(nr),beta0(nr),gamma0(nr),tgrid,L(nr),c0{nr},s_f,xgrid,ae,D,mu);

        % timing
        rt_ep(nr) = toc(t0);
        fprintf(1,'K-S, parameter update, L=%d, nr=%d/%d, time=%6.3f s\n', L(nr), nr, NR, rt_ep(nr));
    end %if

end %nr


