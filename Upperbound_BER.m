clc
clear all
tic
% number of subcarriers
N = 16;
% size of constellation
M_mod = 2
M_bits = log2(M_mod);
N_bits_perfram = N*M_bits;
QAM_table = qammod(0:M_mod-1,M_mod,'gray');
n = N;            % Number of loops
v_1 = ones (1, n);
v_2 = ones (1, n);   % Index vector
Q = length(QAM_table);
Output = zeros (size(repmat (Q, n)));
ready_1 = false;
ready_2 = false;
metric = 100000;

SNR_dB =0:5:40;[0 5 10 15 20 25 30 35 40];0:5:30;
SNR = 10.^(SNR_dB/10);

PEP_bits = 0;
PEP_bits_ave_T = zeros(1, length(SNR));
c0 = 0;
c1 = 1/2/N;
c2 = 1/2/N;
taps_set = 2;
P = taps_set;
delay_taps_set =  [0 3];
% Doppler_taps_set = [1/N  0/N;
%     1/N  1/N;
%     1/N  2/N;
%     1/N  3/N;
%     1/N  4/N;
%     1/N  5/N;
%     1/N  6/N;
%     1/N  7/N];
Doppler_taps_set = [1/N 14/N;
    1/N 4/N;
    1/N 5/N];

h_iter = 1e3;

[taps,delay_taps,Doppler_taps,chan_coef] = OTFS_channel_gen(taps_set(1), delay_taps_set(1, :), Doppler_taps_set(1, :));
[H_eq_AFT H1 H2]= H_eq_AFT_calc(N, c0, c1, c2,taps,delay_taps,Doppler_taps,chan_coef);
ali = 1;
for iesn0 = 1:length(SNR_dB)
    ready_1 = false;
    while ~ ready_1
        index_1 = sub2ind(size (Output), v_1)
        x_i = QAM_table(index_1);
        x_i_demod = qamdemod(x_i,M_mod,'gray');
        x_i_bits = reshape(de2bi(x_i_demod,M_bits),N_bits_perfram,1);
        ali
        ali = 1;
        ready_2 = false;
        PEP_bits = 0;
        PEP_bits_h = 0;
        while ~ ready_2
            index_2 = sub2ind(size (Output), v_2);
            x_j = QAM_table(index_2);
            x_j_demod = qamdemod(x_j,M_mod,'gray');
            x_j_bits = reshape(de2bi(x_j_demod,M_bits),N_bits_perfram,1);
            d_xi_xj = sum(xor(x_i_bits,x_j_bits));
            delta = (x_i - x_j).';
            Phi_delta = [H1*delta H2*delta];
            lambda = svd(Phi_delta);
            if ~isempty(find(abs(lambda) < 1e-5))
                delta
                Phi_delta
                ali = ali + 1;
            end
            %%
            pow_prof = (1/taps) * (ones(1,taps));
            
            [U,S,V] = svd(Phi_delta'*Phi_delta);
            PEP_h_iter = 0;
            for i = 1:h_iter
                temp = 0;
                chan_coef = (sqrt(pow_prof).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)))).';
                h_tilde = U'*chan_coef;
                for j = 1:length(lambda)
                    temp = temp + lambda(j)^2*abs(h_tilde(j))^2;
                end
                PEP_h_iter = PEP_h_iter + qfunc(sqrt(SNR(iesn0)*temp/2));
            end
            PEP_h = PEP_h_iter/h_iter;
            %%
            PEP = 1;
            for j = 1:length(lambda)
                PEP = PEP*(1)/(1 + SNR(iesn0)*(lambda(j)^2)/(4*P));
            end
            %%
            PEP_bits = PEP_bits + (d_xi_xj/(N*M_bits))*PEP;
            PEP_bits_h = PEP_bits_h + (d_xi_xj/(N*M_bits))*PEP_h; 
            % Update the index vector:
            ready_2 = true;       % Assume that the WHILE loop is ready
            for k = 1: n
                v_2 (k) = v_2 (k) + 1;
                if v_2 (k) <= Q
                    ready_2 = false;
                    break ;          % v (k) increased successfully, leave "for k" loop
                end
                v_2 (k) = 1;         % v (k) reached the limit, reset it
            end
        end
        PEP_bits_ave(iesn0) = PEP_bits
        PEP_bits_ave_h(iesn0) = PEP_bits_h
        PEP_bitst_ave_T(iesn0) = PEP_bits_ave_T(iesn0) + PEP_bits_ave(iesn0);
        % Update the index vector:
        ready_1 = true;       % Assume that the WHILE loop is ready
        %     for k = 1: n
        %         v_1 (k) = v_1 (k) + 1;
        %         if v_1 (k) <= Q
        %             ready_1 = false;
        %             break ;          % v (k) increased successfully, leave "for k" loop
        %         end
        %         v_1 (k) = 1;         % v (k) reached the limit, reset it
        %     end
    end
end
PEP_bits_ave_ave = (1/(Q^N))*PEP_bits_ave
toc

