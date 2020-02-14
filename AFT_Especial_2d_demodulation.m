function y = AFT_Especial_2d_demodulation(N,M,r_2d_especial_AFT, c0, c1)
Gamma_c1_MN = diag(exp(-1i*2*pi*(c1.*(0:M*N-1).^2)));
Gamma_c0_MN = diag(exp(1i*2*pi*(c0.*(0:M*N-1))));
F_N = (1/sqrt(N))*dftmtx(N);
F_M = (1/sqrt(M))*dftmtx(M);

r_after_chirp = Gamma_c1_MN*Gamma_c0_MN*r_2d_especial_AFT;
R  = reshape(r_after_chirp, M, N);
y = F_M*R*F_N;
end