function [Uk,Wk,sum_rate] = rate(N_I,K_I,H_k_bar,F,sigma_sq,w_k,d_s)

Jkinner = zeros(N_I,N_I,K_I);

Jk = zeros(N_I,N_I,K_I);

SINR = zeros(N_I,N_I,K_I);

Uk = zeros(N_I,d_s,K_I);

Ek = zeros(d_s,d_s,K_I);

Wk = zeros(d_s,d_s,K_I);

RK = zeros(1,K_I);

for K=1:K_I
    
    Sum14_temp = zeros(N_I,N_I);
    
    for m=1:K_I
        
        if m~=K
            
            Jkinner(:,:,K)=Jkinner(:,:,K)+H_k_bar(:,:,K)*F(:,:,m)*F(:,:,m)'*H_k_bar(:,:,K)';
            
        end
        
        Sum14_temp(:,:) = Sum14_temp(:,:)+H_k_bar(:,:,K)*F(:,:,m)*F(:,:,m)'*H_k_bar(:,:,K)';
        
    end
    
    Sum14=Sum14_temp+sigma_sq*eye(N_I);
    
    Jk(:,:,K)=Jkinner(:,:,K)+(sigma_sq*eye(N_I));
    
    SINR(:,:,K)=H_k_bar(:,:,K)*F(:,:,K)*F(:,:,K)'*H_k_bar(:,:,K)'*inv(Jk(:,:,K));
    
    RK(K)=w_k*log2(det(eye(N_I)+SINR(:,:,K)));
    
    Uk(:,:,K)= inv(Jk(:,:,K)+H_k_bar(:,:,K)*F(:,:,K)*F(:,:,K)'*H_k_bar(:,:,K)')*H_k_bar(:,:,K)*F(:,:,K);
    
    Ek(:,:,K) = eye(d_s)-F(:,:,K)'*H_k_bar(:,:,K)'*inv(Sum14)*H_k_bar(:,:,K)*F(:,:,K);
    
    Wk(:,:,K) = inv(Ek(:,:,K));
       
end

sum_rate = real(sum(RK));

end

