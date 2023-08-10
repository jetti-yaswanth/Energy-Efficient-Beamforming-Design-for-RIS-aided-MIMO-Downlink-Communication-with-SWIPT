%-----------------------------Algorithm 1---------------------------------%

function [F_out_A1] = algorithm1(F_temp1,W_in,U_in,H_bar_in,w_k_in,Q_tilda,G_in,N_B,K_I,A_in,R_min,sigma_sq,d_s)

% initialization

A=A_in;

w_k=w_k_in;

U=U_in;

W = W_in;

G=G_in;

H_bar = H_bar_in;

F_temp2=F_temp1;

%------------------- step 1------------------------------%

error_tol = 10^-6;

lambda_l =10^-8;

lambda_u=1e2;

while abs(lambda_l-lambda_u)>=error_tol
    
    %------------------- step 2------------------------------%
    
    lambda = (lambda_l+lambda_u)/2;
    
    %------------------- step 3------------------------------%
    
    sum25 = zeros(d_s,d_s);
    
    for ki=1:K_I
        
        F_k = F_temp2(:,:,ki);
        
        H_k_bar = H_bar(:,:,ki);
        
        U_k = U(:,:,ki);
        
        W_k = W(:,:,ki);
        
        F_k_star_0= pinv(eye(N_B)+A*lambda)*(lambda*w_k*H_k_bar'*U_k*W_k);
        
        sum25_temp = F_k'*G*F_k_star_0;
        
        sum25 = sum25+sum25_temp;
        
    end
    
    cond_25 = 2*real(trace(sum25));
    
    if cond_25>=Q_tilda
        
        mu=0;
        
    else
        
        sum26_Num = zeros(d_s,d_s);
        
        sum26_den = zeros(d_s,d_s);
        
        for ki2=1:K_I
            
            F_k = F_temp2(:,:,ki2);
            
            H_k_bar = H_bar(:,:,ki2);
            
            U_k = U(:,:,ki2);
            
            W_k = W(:,:,ki2);
            
            sum26_Num_temp = w_k*F_k'*G*inv(eye(N_B)+A*lambda)*lambda*H_k_bar'*U_k*W_k;
            
            sum26_Num =  sum26_Num+sum26_Num_temp;
            
            sum26_den_temp = F_k'*G*inv(eye(N_B)+A*lambda)*G*F_k;
            
            sum26_den =  sum26_den+ sum26_den_temp;
            
        end
        
        Num_26 = 2*real(trace(sum26_Num));
        
        den_26 = 2*trace(sum26_den);
        
        mu = (Q_tilda-Num_26)/(den_26);
        
    end
    
    %------------------- step 4------------------------------%
    
    F_k_star = zeros(N_B,d_s,K_I);
    
    R_temp1 = zeros(d_s,d_s);
    
    R_temp2 = zeros(d_s,d_s);
    
    R_temp3 = zeros(d_s,d_s);
    
    R_temp4 = zeros(d_s,d_s);
    
    R_temp5 = zeros(d_s,d_s);
    
    for ki1=1:K_I
        
        F_k = F_temp2(:,:,ki1);
        
        H_k_bar = H_bar(:,:,ki1);
        
        U_k = U(:,:,ki1);
        
        W_k = W(:,:,ki1);
        
        F_k_star(:,:,ki1)= pinv(eye(N_B)+A*lambda)*(w_k*lambda*H_k_bar'*U_k*W_k+mu*G*F_k);
        
        %------------------- step 5------------------------------%
        
        R_obj1= F_k_star(:,:,ki1)'*A*F_k_star(:,:,ki1);
        
        R_obj2=W_k*U_k'*H_k_bar*F_k_star(:,:,ki1);
        
        R_obj3= W_k*F_k_star(:,:,ki1)'*H_k_bar'*U_k;
        
        R_obj4= W_k;
        
        R_obj5= W_k*(U_k)'*U_k;
        
        R_temp1=  R_temp1+ R_obj1;
        
        R_temp2=  R_temp2+ R_obj2;
        
        R_temp3=  R_temp3+ R_obj3;
        
        R_temp4=  R_temp4+ R_obj4;
        
        R_temp5=  R_temp5+ R_obj5;
        
    end
    
    %------------------- step 6------------------------------%
    
    R_lambda = log(det(W(:,:,K_I)))-(trace(R_temp1)-trace(R_temp2)-trace(R_temp3)+trace(R_temp4)+sigma_sq*trace(R_temp5))+d_s;
    
    if R_lambda<=R_min
        
        lambda_l = lambda;
        
    else
        
        lambda_u = lambda;
        
    end
    
end

F_out_A1 = F_k_star;

end
