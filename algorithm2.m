%-------------------------------ALgorithm~2---------------------------------%

function [F_out_A2] = algorithm2(W_in,U_in,H_bar_in,w_k_in,Q_bar_in,G_in,N_B,K_I,F,R_min,sigma_sq,d_s,err_tol)

% Inputs from main

Q_bar = Q_bar_in;

W = W_in;

U = U_in;

H_bar = H_bar_in;

w_k=w_k_in;

%------------------- step 1------------------------------%

% Intialization

n_max= 20;       % maximum iterations

G = G_in;

F_0 = F;

A =zeros(N_B,N_B);

for m=1:K_I
    
    w_m =w_k;
    
    H_m_bar = H_bar;
    
    U_m = U;
    
    W_m = W;
    
    A =A+w_m*H_m_bar(:,:,m)'*U_m(:,:,m)*W_m(:,:,m)*U_m(:,:,m)'*H_m_bar(:,:,m);   % content related to eq. (15)
    
end

% Calculate Intial OF

obj_old = zeros(d_s,d_s);

for ki=1:K_I
    
    obj_old = obj_old+(F_0(:,:,ki)'*F_0(:,:,ki));
    
end

n=1;

s(n) = trace(obj_old);  

%------------------- step 2------------------------------%

F_temp = F_0;

while (n<n_max+1)
    
    sum_step2_alg2 = zeros(d_s,d_s) ;
    
    for kq = 1:K_I
        
        F_k = F_temp(:,:,kq);
        
        sum_temp = F_k'*G*F_k;
        
        sum_step2_alg2 = sum_step2_alg2+sum_temp;
        
    end
    
    Q_tilda = Q_bar+ trace(sum_step2_alg2);
    
%------------------- step 3------------------------------%
    
    [F_opt] = algorithm1(F_temp,W,U,H_bar,w_k,Q_tilda,G,N_B,K_I,A,R_min,sigma_sq,d_s);
    
%------------------- step 4------------------------------%
    
    sum_new = zeros(1,K_I);
    
    for ki1=1:K_I
        
        F_k = F_opt(:,:,ki1);
        
        sum_new(ki1) = trace(F_k*F_k');
        
    end
    
    s(n+1) = sum(sum_new);     
    
    error = abs(s(n+1)-s(n))/abs(s(n+1));       % calculate error
    
    if error<err_tol
        
        n=n_max+1;   % for termination
        
    else
        
        n=n+1;       % next iteration
        
    end
    
end

F_out_A2=F_opt;  

end