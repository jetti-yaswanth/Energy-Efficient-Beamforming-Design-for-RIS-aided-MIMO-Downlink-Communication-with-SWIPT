clear variables
clc

for M = [50,100,150]

% Variable Intilization

sigma_sq_dB = -174;                           % Noise in dBm per Hz

bw = 1e6;                                     % Bandwidth in MHZ

sigma_sq = 10.^((sigma_sq_dB-30)/10)*bw;      % Noise varince

d_s = 2;                                      % Number of data streams at each IR

K_I = 2;                                      % Number of Information receivers

K_E = 4;                                      % Number of Energy Receivers

w_k=1;                                        % weight factor of each ERs, considered equal for all user

N_B = 4;                                      % Number of antenna at BS

N_I = 2;                                      % Number of antenna at infromation Rx

N_E =2;                                       % # antenna at Energy Rx

Q_bar = 2*10^-4;                              % Minimum harvested power threshold

eta = 0.5;                                    % energy harvesting efficeincy

alpha_l =1;                                   % enery weighting factor of ER

R_min = 5;                                    %  min rate constant

err_tol = 10^-4;                              % Error tolerance

x_ER = 5;                                     % x-coordinate of ER in m

x_IR = 400;                                   % x-coordinate of IR in m

x_RIS = 5;                                    % x-coordinate of IRS in m

%M = 50;                                       % Number of RIS elements

% Intilization for channel generation

beta=3;                                       % Rician channel

PL_0= -30;                                    % pathloss per meter in dB

PL_0_re =10.^(PL_0/10);                       % Pathloss in watt

C_b=[0,0];                                    % Distance at BS

C_l = [x_ER,0];                               % Distance from BS to ER

C_k = [x_IR,0];                               % Distance from BS to IR

C_I = [x_RIS,2];                              % Distance from  BS to IRS

% loop intilization

n_itr = 30;                                   % Iterations

sim=200;                                       % simulation

% Start Simulations

power_sim= zeros(sim,n_itr);

parfor j = 1:sim
    
    % Intial Random Precoder Generation
    
    phi =  exp(1i*pi*randn(1,M)');                 % Precoder at RIS
    
    PHI = diag(phi);
    
    F = randn(N_B,d_s,K_I)+1i*randn(N_B,d_s,K_I);  % Precoder at BS
    
    % Channel Generation with Intial Precoders
    
    [H_bk,H_rk,G_bl,G_rl,Z] = channel_random(beta,PL_0_re,C_b,C_l,C_k,C_I,N_I,N_B,K_I,M,N_E,K_E);
    
    % Equivalent channel generation
    
    [H_k_bar,G] = H_G_phi(G_bl,G_rl,Z,eta,alpha_l,H_bk,H_rk,PHI,K_I,K_E,N_B,N_I);
    
    % Rate, Decoding Matrix and Auxilary Matrix Generations
    
    [U,W,sum_rate] = rate(N_I,K_I,H_k_bar,F,sigma_sq,w_k,d_s);
    
    % Start Iteration
    
    power_itr = zeros(1,n_itr);
    
    for i = 1:n_itr
        
        % optimization algorithm for the transmit beamforming at BS
        
        F = algorithm2(W,U,H_k_bar,w_k,Q_bar,G,N_B,K_I,F,R_min,sigma_sq,d_s,err_tol);
        
        % Optimization algorithm for the phase shift at the RIS
        
        phi=algorithm4(G_rl,Z,F,G_bl...
            ,W,H_rk,H_bk,U,w_k,Q_bar,K_I,K_E,M,alpha_l,eta,phi);
        
        % Minimizing transmit power using optimal precoders
        
        pow_temp=0;
        
        for k=1:K_I
            
            pow_temp=pow_temp+norm(F(:,:,k),'fro');
            
        end
        
        % update the channels using the optimal functions
        
        [H_k_bar,G] = H_G_phi(G_bl,G_rl,Z,eta,alpha_l,H_bk,H_rk,PHI,K_I,K_E,N_B,N_I);
        
        % Updat the Decoding matrix, Auxilary matrix with obatined opt values of F, phi
        
        [U,W,sum_rate] = rate(N_I,K_I,H_k_bar,F,sigma_sq,w_k,d_s);
        
        % Get the opt Tranmit power in each iteration
        
        power_itr(i) = pow_temp;
        
    end
    
    % get the transmit power over simulations and iterations
    
    power_sim(j,:) = real(power_itr);
    
end

opt_power = mean(power_sim);

% figure for convergence

hold on
plot(1:n_itr,smooth(opt_power),'DisplayName',strcat('M=',num2str(M)));
xlabel('${Number of Iterations}$','interpreter','latex')
ylabel('${Transmit Power}$ (W)','interpreter','latex')
legend('show')

end