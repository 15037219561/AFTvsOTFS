function H_EQ = H_EQ_2d_AFT_Calc(N,M, c1, c2, c3, c4, H)
Gamma_c4_N= diag(exp(-1i*2*pi*(c4.*(0:N-1).^2)));
Gamma_c3_N= diag(exp(-1i*2*pi*(c3.*(0:N-1).^2)));
    
Gamma_c1_M= diag(exp(-1i*2*pi*(c1.*(0:M-1).^2)));
Gamma_c2_M= diag(exp(-1i*2*pi*(c2.*(0:M-1).^2)));

F_N = (1/sqrt(N))*dftmtx(N);
F_M = (1/sqrt(M))*dftmtx(M);

A_N = Gamma_c4_N*F_N*Gamma_c3_N;
A_M = Gamma_c2_M*F_M*Gamma_c1_M;

H_EQ = kron(A_N, A_M)*H*kron(A_N', A_M');
end