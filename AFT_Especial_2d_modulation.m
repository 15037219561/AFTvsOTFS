function s_with_cp = AFT_Especial_2d_modulation(N,M,X, N_CP, c1)

Gamma_c1_MN = diag(exp(-1i*2*pi*(c1.*(0:M*N-1).^2)));

F_N = (1/sqrt(N))*dftmtx(N);
F_M = (1/sqrt(M))*dftmtx(M);
s_mat = F_M'*X*F_N';
s = s_mat(:);
s_with_Chirp = Gamma_c1_MN'*s;
%Add CP
v = N_CP:-1:1;
s_with_cp = [s_with_Chirp(M*N-N_CP+1:M*N).*exp(-1i*2*pi*c1.*((M*N)^2 - 2*M*N*(v))); s_with_Chirp];
end