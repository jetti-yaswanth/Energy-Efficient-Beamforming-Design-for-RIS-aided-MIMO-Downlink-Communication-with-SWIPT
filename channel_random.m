function [H_bk,H_rk,G_bl,G_rl,Z] = channel_random(beta,PL_0_re,C_b,C_l,C_k,C_I,N_I,N_B,K_I,M,N_E,K_E)

% pathloss exponent for each path

alpa_gbl =3.6;                                % BS-ER

alpa_grl =2.2;                                % RIS-ER

alpa_z = 2.2;                                 % BS-RIS

alpa_hbk =3.6;                                % BS-IR

alpa_hrk =2.4;                                % RIS-IR

% Distances  between devices

d_gbl =norm(C_b-C_l);

d_grl =norm(C_I-C_l);

d_z = norm(C_b-C_I);

d_hbk =norm(C_b-C_k);

d_hrk =norm(C_I-C_k);

% calculation of pathloss

PL_gbl =sqrt(PL_0_re*(d_gbl^(-alpa_gbl)));

PL_grl =sqrt(PL_0_re*(d_grl^(-alpa_grl)));

PL_z = sqrt(PL_0_re*(d_z^(-alpa_z)));

PL_hbk =sqrt(PL_0_re*(d_hbk^(-alpa_hbk)));

PL_hrk =sqrt(PL_0_re*(d_hrk^(-alpa_hrk)));

% Rayleigh fading channels

H_bk = (PL_hbk/sqrt(2)).*((randn(N_I,N_B,K_I)+1i*randn(N_I,N_B,K_I)));

H_rk = (PL_hrk/sqrt(2)).*((randn(N_I,M,K_I)+1i*randn(N_I,M,K_I)));

% Rician fading channels

G_bl = rician_ch(N_E,N_B,K_E,beta,PL_gbl);

G_rl = rician_ch(N_E,M,K_E,beta,PL_grl);

Z = rician_ch(M,N_B,1,beta,PL_z);

end