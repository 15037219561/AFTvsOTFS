clc
clear all

c1 = 1.1;
c2 = 2.2;
c3 = 0;
c4 = 4.4;
M = 4;
N =2;
L = [1 2];
P = length(L);%number of Pathes
power_channel_dB    = [0 0];
amp_k               = [1.5 2.58];
doppler_spread      =[3.12 2.13];
f0                  = doppler_spread(1);
f1                  = doppler_spread(2);
l0                  = L(1);
l1                  = L(2);
c0                  = (f0*l1 - f1*l0)/(l1 - l0)
c1                  = (f0-f1)/(2*(l1 - l0));


Num_OFDM_sym = 1;
N_CP = 0;

Gamma_c4_N = diag(exp(-1i*2*pi*(c4.*(0:N-1).^2)));
Gamma_c1_N = diag(exp(-1i*2*pi*(c1.*(0:N-1).^2)));
Gamma_c2_N = diag(exp(-1i*2*pi*(c2.*(0:N-1).^2)));
Gamma_c3_N = diag(exp(-1i*2*pi*(c3.*(0:N-1).^2)));

Gamma_c1_M = diag(exp(-1i*2*pi*(c1.*(0:M-1).^2)));
Gamma_c2_M = diag(exp(-1i*2*pi*(c2.*(0:M-1).^2)));
Gamma_c3_M = diag(exp(-1i*2*pi*(c3.*(0:M-1).^2)));
Gamma_c4_M = diag(exp(-1i*2*pi*(c4.*(0:M-1).^2)));


F_N = (1/sqrt(N))*dftmtx(N);
F_M = (1/sqrt(M))*dftmtx(M);

A1_N = Gamma_c2_N*F_N*Gamma_c1_N;
A2_N = Gamma_c4_N*F_N*Gamma_c3_N;

A1_M = Gamma_c2_M*F_M*Gamma_c1_M;
A2_M = Gamma_c4_M*F_M*Gamma_c3_M;

A_N = Gamma_c4_N*F_N*Gamma_c3_N;
A_M = Gamma_c2_M*F_M*Gamma_c1_M;


x = randn(N, M);

% X = fft(ifft(x).').'/sqrt(M/N); %%%ISFFT
% s_mat = ifft(X.')*sqrt(M); % Heisenberg transform

X = (A1_M*(A1_N'*x).').'; %%%ISFFT
s_mat = ifft(X.')*sqrt(M); % Heisenberg transform

r = s_mat;
r_mat = reshape(r,M,N);
Y = fft(r_mat)/sqrt(M); % Wigner transform
Y = Y.';
y = (A1_M'*(A1_N*Y).').'; %%%ISFFT



eye_MN = eye(M*N);
Gamma_mat = [eye_MN(:, 2:end) eye_MN(:, 1)];%Permutation Matrix
H = zeros(M*N);
for i= 0:P-1
    H =  H + (i+1)*Gamma_mat^L(i+1);
end
H_eq = kron(Gamma_c3_N, Gamma_c1_M)*H*kron(Gamma_c3_N, Gamma_c1_M)';



G_n_l       = zeros(M*N);
for n = 0:M*N-1
    for l = 0:length(L)-1%Number of Cluster
        A_i_n                           = amp_k(l+1);% It can be RV
        f_i_n                           = doppler_spread(l+1);
        G_n_l(n+1, L(l+1)+1)   = A_i_n*exp(-1i*2*pi*f_i_n*n);
    end
end



H_eq_test = zeros(M*N);
for l = 0:P-1
    for n = 0:N*M-1
        n;
        l;
        if n==3 && mod(n-L(l+1), M*N)==2
            ali = 1;
        end
        k = floor(n/M);
        if (n - L(l+1) < k*M)
            if k == 0
                Coeff  = exp(-1i*2*pi*(2*c1*(L(l+1) - M)*n - c3*(N-1)^2 -c1*(M-L(l+1))^2));
                Coeff_test  = exp(-1i*2*pi*c1*(M^2 + 2*M*(n - L(l+1))) + 1i*2*pi*c3*(N-1)^2);
            else
                Coeff  = exp(-1i*2*pi*(2*c1*(L(l+1) - M)*n + c3*(2*k-1) -c1*(L(l+1) - M)*(L(l+1) + (2*k-1)*M)));
                Coeff_test  = exp(-1i*2*pi*c1*(-(2*k-1)*M^2 + 2*M*(n-L(l+1))) + 1i*2*pi*c3*(2*k-1));
            end
        else
            Coeff  = exp(-1i*2*pi*(2*c1*n*L(l+1)-c1*L(l+1)*(2*k*M+L(l+1))));
            Coeff_test = 1;
        end
        A_i_n                           = amp_k(l+1);% It can be RV
        f_i_n                           = doppler_spread(l+1);
        f = exp(1i*2*pi*(-c1*L(l+1)*(2*k*M+L(l+1))));
        h = A_i_n*exp(-1i*2*pi*f_i_n*n)*Coeff*Coeff_test;
        H_eq_test(n + 1, mod(n-L(l+1), M*N) + 1) = H(n + 1, mod(n-L(l+1), M*N) + 1)*Coeff;
    end
end
H_eq_test;
abs(H_eq_test - H_eq);
[a,b] = find(abs(H_eq_test - H_eq)> 0.1);


%taking to account CP
H_eq_pure = zeros(M*N);
H_n_l = zeros(M*N);

Gamma_c1_MN = diag(exp(-1i*2*pi*(c1.*(0:M*N-1).^2)));
Gamma_c0_MN = diag(exp(1i*2*pi*(c0.*(0:M*N-1))));
for l = 0:P-1
    for n = 0:N*M-1
        n;
        l;
        k = floor(n/M);
        Coeff_AL2 = 1;
        if (n - L(l+1) < k*M)
            if k == 0
                Coeff  = exp(-1i*2*pi*c1*(M^2 + 2*M*(n - L(l+1))) + 1i*2*pi*c3*(N-1)^2);
                Coeff_AL2 = exp(-1i*2*pi*c1*((M*N)^2 + 2*M*N*(n - L(l+1))));
            else
                Coeff  = exp(-1i*2*pi*c1*(-(2*k-1)*M^2 + 2*M*(n-L(l+1))) + 1i*2*pi*c3*(2*k-1));
            end
        else
            Coeff = 1;
        end
        temp = exp(1i*2*pi*(-c1*L(l+1)*(2*k*M+L(l+1))))
        H_eq_pure(n + 1, mod(n-L(l+1), M*N) + 1) = G_n_l(n+1, L(l+1)+1)*Coeff;
        H_n_l(n + 1, mod(n-L(l+1), M*N) + 1) = G_n_l(n+1, L(l+1)+1)*Coeff_AL2;
    end
end
H_EQ = kron(Gamma_c3_N, Gamma_c1_M)*H_eq_pure*kron(Gamma_c3_N, Gamma_c1_M)'*Gamma_c0_MN

H_EQ_AL2 = Gamma_c0_MN*Gamma_c1_MN*H_n_l*Gamma_c1_MN'
kron(F_N, F_M)*H_EQ_AL2*kron(F_N', F_M')
kron(F_N, eye(M))*H_EQ_AL2*kron(F_N', eye(M))
kron(eye(N), F_M)*H_EQ_AL2*kron(eye(N), F_M')

%%
clc
N = 4;
M = 2;
X = [1 2 3 4; 5 6 7 8];
A = [-1 -2 -3 -4; -3 -5 -7 -9];
F_N = (1/sqrt(N))*dftmtx(N);
F_M = (1/sqrt(M))*dftmtx(M);


S = A.*(F_M'*X*F_N');
vecS = diag(A(:))*kron(F_N',F_M')*X(:);
vecS - S(:)
