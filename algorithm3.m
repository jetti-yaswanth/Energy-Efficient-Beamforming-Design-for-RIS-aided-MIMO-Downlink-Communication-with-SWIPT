
%-------------------- Algorithm~3---------------------------%

function [phi_out_A3]= algorithm3(g_star,Q_cap,q,Li,phi_temp)

% inputs from Algorithm 4

phi=phi_temp;

%------------------------------- Step 1-----------------------------------%

phi_p_0= exp(1i*angle(q)); 

J_0= 2*real(phi_p_0'*(g_star+Li*phi));

if   J_0<=Q_cap 
    
    phi_out_A3 = phi; 
    
else
    
    %---------------------------- Step 2--------------------------------%
    
    % intializtion
    
    error_tol = 10^-6;  %error tolerance
    
    p_l = 10;
    
    p_u = 60;
    
    while abs(p_l-p_u)>= error_tol % condtiton for error tolerance
        
   %------------------------ Step 3 ------------------------------%
   
        p = (p_l+p_u)/2;
        
   %------------------------- step 4-------------------------------%
        
        phi_p = exp(1i*angle(q+p*(g_star+Li*phi)));      
        
        J_p = 2*real(phi_p'*(g_star+Li*phi));     
        
   %--------------------------- Step 5--------------------------------%
        
        if J_p>=Q_cap % condition for bisection search
            
            p_u = p;
            
        else
            
            p_l = p;
            
        end
        
    end
    
    phi_out_A3 = phi_p;
    
end

end
