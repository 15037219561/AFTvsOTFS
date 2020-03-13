clear all
clc




N = 10;
n = 0:N-1;
delay = [1 ];
amp_k = [1 ];
doppler_spread = [1/N];
F_l = doppler_spread(1);
L = 1;
c1 =  1/(2*N);
c2 =  1/(2*N);
Gamma_c1_N = diag(exp(-1i*2*pi*(c1.*(0:N-1).^2)));
Gamma_c2_N = diag(exp(-1i*2*pi*(c2.*(0:N-1).^2)));
F_N = (1/sqrt(N))*dftmtx(N);



EyeMat = eye(N);
Per = [EyeMat(:, 2:end) EyeMat(:,1)];
H = zeros(N);
for l = 0:length(delay)-1%Number of Cluster
    A_i_n                           = amp_k(l+1);% It can be RV
    f_i_n                           = doppler_spread(l+1);
    
    
    Delta_ki = diag(exp(-1i*2*pi.*n*f_i_n));
    onesCP = ones(1, N);
    onesCP(1:delay(l+1)) = fliplr(exp(-1i*2*pi*c1*(N^2 - 2*N*(1:delay(l+1)))));
    Gamma_CP = diag(onesCP);
    H = H + A_i_n*Delta_ki*Gamma_CP*Per^(delay(l+1));
    H_each(:, :, l+1) = A_i_n*Delta_ki*Gamma_CP*Per^(delay(l+1));
end
A = Gamma_c1_N*H*Gamma_c1_N';

for i = 0:N-1
    for j = 0:N-1
        test = 0;
        test2 = 0;
        test3 = 0;
        for n = 0:N-1
            test = test + exp(-1i*2*(pi/N)*(i*n - mod(n-delay(L), N)*j))*A(n+1, mod(n-delay(L), N)+1);
            test2 = test2 + exp(-1i*2*(pi/N)*(i*n - (n-delay(L))*j))*exp(1i*2*(pi/N)*(c1*N*(delay(L)^2 - 2*n*delay(L))-F_l*n*N));
            test3 = test3 + exp(-1j*(2*pi/N)*(i - j + N*F_l + 2*N*c1*delay(L))*n)
        end
        BB_L1(i+1, j+1)=1/N*test;
        BBB_L1(i+1, j+1)=1/N*test2;
        BBB_2(i+1, j+1) = (1/N)*exp(1j*2*(pi/N)*(N*c1*delay(L)^2 - j*delay(L)))*test3;
    end
end

B = F_N*A*F_N';
BB_L1;
max(max(B - BB_L1));
max(max(BBB_L1 - BB_L1));
max(max(B - BBB_2))

C = Gamma_c2_N*B*Gamma_c2_N';
%Gamma_c2_N*F_N*Gamma_c1_N*H_each(:, :, 2)*Gamma_c1_N'*F_N'*Gamma_c2_N';

%testing

for L = 1:16%100
    A = zeros(20,5);
    for n =32:32 %100:10000
        %N = 2^n;
        N = n;
        c2 = 1/(2*N);
        X = 1:N-1;
        x = 0;
        if ~isempty(find(mod(X.*(2*N*c2.*X + L/2), N) == 0))
            N;
            x = X(find(mod(X.*(2*N*c2.*X + L/2), N) == 0))
            N ;
            L
            a = min(length(find(factor(L) == 2)), n-1);
            A(n, :) = [n, L, length(x), length(x)/N, 2^a-1];
        end
    end
    A(:, :);
    [a b] = max(A(:, 3))
    if (a ~= 0)
        ali = 1;
    end
end

%%
p = 0.2;
n = 10;
a = 5;
nv = 5;
v = zeros (1, nv);
p1 = 0;
ready = false;
while ~ ready
    p1 = p1 + prod (p * (1-p).^v);
    % Update the index vector:
    ready = true;       % Assume that the WHILE loop is ready
    for k = nv: -1: 1
        v (k) = v (k) + 1;
        if v (k) <= n - a - sum (v (1: k-1))
            % This is your "0: (na -ijkl -...)" criterion
            ready = false;  % No, WHILE loop is not ready now
            break ;          % v (k) increased successfully, leave "for k" loop
        end
        v (k) = 0;  % v (k) reached the limit, reset it and iterate v (k-1)
    end
end
ali = 1


%%
x1 = [1 2 3 4]
x2 = [1 2 3 4]
for i = 1:length(x1)
    for j =1:length(x2)
        diff(i, j) = x1(i) - x2(j);
    end
end

load('H1.mat')
load('H2.mat')
N = 8
M_mod = 4
M_bits = log2(M_mod);
N_bits_perfram = M_bits*N;
cntr1 = 0;
cntr2 = 0;
for i = 1:1%M_mod^N
    data_info_bit = randi([0,1],N_bits_perfram,1);
        data_temp = bi2de(reshape(data_info_bit,N,M_bits));
        x1 = qammod(data_temp,M_mod,'gray');
    for j = 1:M_mod^N
        data_info_bit = randi([0,1],N_bits_perfram,1);
        data_temp = bi2de(reshape(data_info_bit,N,M_bits));
        x2 = qammod(data_temp,M_mod,'gray');
        delta = x1 - x2;
        Phi = [H1*delta H2*delta];
        if rank(Phi) == 1
            cntr1 = cntr1 + 1;
        else
            cntr2 = cntr2 + 1;
        end
        %x = reshape(x,M,N);
    end
end