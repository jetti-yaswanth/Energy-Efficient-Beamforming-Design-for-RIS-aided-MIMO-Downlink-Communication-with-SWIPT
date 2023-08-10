function [G_l]=rician_ch(N_E,M,K_E,beta,PL)

belta1 = sqrt((beta)/(1+beta));

belta2 = sqrt(1/(1+beta));

G_l = zeros(N_E,M,K_E);

for l = 1:K_E
    
    a_N_E = exp(1i*pi*(0:1:N_E-1)'*sin(pi*rand-pi/2));
    
    a_M = exp(1i*pi*(0:1:M-1)'*sin(pi*rand-pi/2));
    
    G_LOS = a_N_E*a_M';
    
    G_l(:,:,l) = (PL/sqrt(2)).*(belta1*G_LOS+belta2*(randn(N_E,M)+1i*randn(N_E,M)));
    
end

end
