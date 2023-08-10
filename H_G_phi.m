function [H_k_bar,G] = H_G_phi(G_bl,G_rl,Z,eta,alpha_l,H_bk,H_rk,PHI,K_I,K_E,N_B,N_I)

G = zeros(N_B,N_B);

for l=1:K_E           
    
    G_l_bar = G_bl(:,:,l)+G_rl(:,:,l)*PHI*Z;
    
    G_temp = alpha_l*eta*(G_l_bar'*G_l_bar);
    
    G=G+G_temp;
    
end

H_k_bar = zeros(N_I,N_B,K_I);

for k=1:K_I
    
    H_k_bar(:,:,k) = H_bk(:,:,k) +H_rk(:,:,k)*PHI*Z;       
    
end

end