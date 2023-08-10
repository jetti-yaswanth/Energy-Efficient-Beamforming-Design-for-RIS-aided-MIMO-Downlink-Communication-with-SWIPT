
%-----------------------Algorithm~4 ---------------------------%

function [phi_out_A4]=algorithm4(G_rl,Z,F,G_bl...
    ,W ,H_rk ,H_bk,U,w_k,Q_bar,K_I,K_E,M,alpha_l,eta,phi)

%---------------------- Step 1-------------------------------%

%Initialization

err_tol = 10^-6;       % error tolerance

n_max= 40;             % maximum iterations

lambda_max=10;


F_tilda=0;

for m = 1:K_I
    
    F_m = F(:,:,m);
    
    F_tilda = F_tilda + F_m*F_m';   
    
end

Gr=0;

Gb=0;

Gbr_temp =0;

for l=1:K_E
    
    Gr_temp = alpha_l*eta*G_rl(:,:,l)'*G_rl(:,:,l);
    
    Gr=Gr+Gr_temp;
    
    Gb_temp = alpha_l*eta*G_bl(:,:,l)'*G_bl(:,:,l);
    
    Gb=Gb+Gb_temp;
    
    Gbr_temp1 = alpha_l*eta*G_bl(:,:,l)'*G_rl(:,:,l);
    
    Gbr_temp = Gbr_temp+Gbr_temp1;
    
end

Gbr= Z*F_tilda*Gbr_temp;

C = Z*F_tilda*Z';   

Li=Gr.*C.';      

for k = 1:K_I
    
    B_k(:,:,k) = w_k*H_rk(:,:,k)'*U(:,:,k)*W(:,:,k)*U(:,:,k)'*H_rk(:,:,k);   
    
    %B = B + B_k;
    
    D_k(:,:,k) = w_k*Z*F_tilda'*H_bk(:,:,k)'*U(:,:,k)*W(:,:,k)*U(:,:,k)'*H_rk(:,:,k);    
    
    T_k(:,:,k) = w_k*Z*F(:,:,k)*W(:,:,k)*U(:,:,k)'*H_rk(:,:,k);      
    
end

%T_k = w_k*Z*F*W_k*U_k'*H_rk;      %Equation 32 and Related Content

B = sum(B_k,3);

Xi = B.*C.';

Sum1= sum(D_k,3);

Sum2= sum(T_k,3);

V=Sum1-Sum2;                         

v_star=conj(diag(V));                     

g_star=conj(diag(Gbr));                     

% step 1

  phi_temp=phi;

n=1;

f(n) = phi_temp'*Xi*phi_temp+2*real(phi_temp'*v_star);      

while n<=n_max
    
    %-------------- Step 2 --------------------%
    
    Q_hat= Q_bar-trace(Gb*F_tilda);         
    
    Q_cap = Q_hat+ phi_temp'*Li*phi_temp;   
    
    % Step 3
    
    q = (lambda_max*eye(M)-Xi)*phi_temp-v_star;   
    
    % Step 4
    
    phi_temp=algorithm3(g_star,Q_cap,q,Li,phi_temp);
    
    % Step 5
    
    f(n+1) = phi_temp'*Xi*phi_temp+2*real(phi_temp'*v_star);
    
    error = abs(f(n+1)-f(n))/f(n+1); %calculate error
    
    if error<= err_tol
        
        n=n_max+1; %for termination
        
    else
        
        n=n+1;
        
    end
    
end

phi_out_A4 = phi_temp;

end

