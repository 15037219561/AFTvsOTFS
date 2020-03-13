function [H_eq_AFT H1 H2] = H_eq_AFT_calc(N_AFT, c0, c1, c2,taps,delay_taps,Doppler_taps,chan_coef)
n = 0:N_AFT-1;
Gamma_c1_N = diag(exp(-1i*2*pi*(c1.*(0:N_AFT-1).^2)));
Gamma_c2_N = diag(exp(-1i*2*pi*(c2.*(0:N_AFT-1).^2)));
Gamma_c0_N = diag(exp(1i*2*pi*(c0.*(0:N_AFT-1))));
F_N = (1/sqrt(N_AFT))*dftmtx(N_AFT);
EyeMat = eye(N_AFT);
Per = [EyeMat(:, 2:end) EyeMat(:,1)];
H = zeros(N_AFT);
for l = 0:taps-1%Number of Cluster
    A_i_n                           = chan_coef(l+1);% It can be RV
    f_i_n                           = Doppler_taps(l+1);
    Delta_ki = diag(exp(-1i*2*pi.*n*f_i_n));
    onesCP = ones(1, N_AFT);
    onesCP(1:delay_taps(l+1)) = fliplr(exp(-1i*2*pi*c1*(N_AFT^2 - 2*N_AFT*(1:delay_taps(l+1)))));
    Gamma_CP = diag(onesCP);
    H = H + A_i_n*Delta_ki*Gamma_CP*Per^(delay_taps(l+1));
    if l == 0
        H1 = Delta_ki*Gamma_CP*Per^(delay_taps(l+1));
    else
        H2 = Delta_ki*Gamma_CP*Per^(delay_taps(l+1));
    end
end

H_eq_AFT = Gamma_c2_N*F_N*Gamma_c0_N*Gamma_c1_N*H*Gamma_c1_N'*F_N'*Gamma_c2_N';
H1 = Gamma_c2_N*F_N*Gamma_c0_N*Gamma_c1_N*H1*Gamma_c1_N'*F_N'*Gamma_c2_N';
H2 = Gamma_c2_N*F_N*Gamma_c0_N*Gamma_c1_N*H2*Gamma_c1_N'*F_N'*Gamma_c2_N';

