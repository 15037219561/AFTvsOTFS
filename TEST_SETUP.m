clc
clear all
close all

warning('off','all')
tic
%% Enabling Options
enable_AFT = true;
enable_OTFS = false;
enable_2d_AFT = false;
enable_Especial_2d_AFT = false;
enable_OTFS_LMMSE = false;
SC_AFT = false;
%% OTFS parameters

% number of symbol
N = 8
% number of subcarriers
M = 8
% size of constellation
M_mod = 64
M_bits = log2(M_mod);
% number of symbols per frame
N_syms_perfram = N*M;
% number of bits per frame
N_bits_perfram = N*M*M_bits;


%% AFT parameters

% number of subcarriers
N_AFT = M;
N_AFT_SC = N_AFT;
SC_fac = 2;
% number of AFT symbol
Num_AFT_sym = N;
% noise poser
SNR_dB = [0 5 15 17 20 22 25];0:5:30;
SNR = 10.^(SNR_dB/10);
noise_var_sqrt = sqrt(1./SNR);

% Singal Power--> calculated in iesn0 = 0
sig_energy_OTFS = 0;
sig_energy_2d_especial_AFT = 0;
sig_energy_AFT = 0;

rng(1)
N_fram = 1000;%10^4;
err_ber_OTFS = zeros(length(SNR_dB),1);
err_ber_2d_AFT = zeros(length(SNR_dB),1);
err_ber_2d_especial_AFT = zeros(length(SNR_dB),1);
err_ber_AFT = zeros(length(SNR_dB),1);
for iesn0 = 0:length(SNR_dB)
    for ifram = 1:N_fram
        %% random input bits generation
        
        data_info_bit = randi([0,1],N_bits_perfram,1);
        data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
        %x = qammod(data_temp,M_mod,0, 'gray');
        x = qammod(data_temp,M_mod,'gray');
        %x = reshape(x,N,M);
        x = reshape(x,M,N);
        %% channel generation
        
        [taps,delay_taps,Doppler_taps,chan_coef] = OTFS_channel_gen(N,M);
        N_CP = max(delay_taps);
        % for the moment, we assume two-tap delay channel
        if taps == 2
            [c0, c1, c2] = ComputeC0_C1_for2path(Doppler_taps, delay_taps);
            c3 = 1;
            c4 = 1;
            %             c1 = 1/(2*M);
            %             c2 = 1/(2*M);
            %             c0 = 0;
            %             c0 = 0;
            %             c1 = 0;
            %             c2 = 0;
        end
        
        %% Modulation
        
        if enable_OTFS
            % OTFS modulation
            s_OTFS = OTFS_modulation(N,M,x);
        end
        if enable_2d_AFT
            % 2d_AFT modulation
            s_2d_AFT = AFT_2d_modulation(N,M,x, c1, c2, c3, c4);
            %s_2d_AFT = AFT_2d_modulation(N,M,x, 1/(2*N), 1/(2*N), 0, 0);
        end
        if enable_Especial_2d_AFT
            % 2d_AFT modulation
            s_2d_Especial_AFT = AFT_Especial_2d_modulation(N,M,x, N_CP, c1);
            %s_2d_AFT = AFT_2d_modulation(N,M,x, 1/(2*N), 1/(2*N), 0, 0);
        end
        if enable_AFT
            % AFT modulation
            if SC_AFT
                N_AFT = N_AFT_SC*SC_fac;
                x = fft(x, N_AFT);
            end
            s_AFT = AFT_modulation(N_AFT,Num_AFT_sym, N_CP, c1, c2, x);
        end
        %% Calculate the Signal Energy
        
        sig_energy = 0;
        if iesn0 == 0
            if enable_OTFS
                sig_energy = OTFS_Sig_energy(N,M,taps,delay_taps,Doppler_taps,chan_coef,s_OTFS);
                sig_energy_OTFS = sig_energy_OTFS + sig_energy;
            end
            if enable_2d_AFT
                sig_energy = OTFS_Sig_energy(N,M,taps,delay_taps,Doppler_taps,chan_coef,s_OTFS);
                sig_energy_OTFS = sig_energy_OTFS + sig_energy;
            end
            if enable_Especial_2d_AFT
                sig_energy = AFT_especial_2d_Sig_energy(N,M,taps,delay_taps,Doppler_taps,chan_coef,s_2d_Especial_AFT);
                sig_energy_2d_especial_AFT = sig_energy_2d_especial_AFT + sig_energy;
            end
            
            
            if enable_AFT
                % AFT
                sig_energy = AFT_Sig_energy(N_AFT, Num_AFT_sym, taps, delay_taps, Doppler_taps, chan_coef,s_AFT);
                sig_energy_AFT = sig_energy_AFT + sig_energy;
            end
            continue;
        end
        %% channel output
        
        if enable_OTFS
            % OTFS
            % H_OTFS_eq is the equivalent channel matrix which is used for
            % the MMSE equalizer
            [r_OTFS, H_OTFS_eq] = OTFS_channel_output(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2_OTFS(iesn0),s_OTFS);
        end
        if enable_2d_AFT
            [r_2d_AFT, H_2d_AFT] = OTFS_channel_output(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2_OTFS(iesn0),s_2d_AFT);
        end
        
        if enable_Especial_2d_AFT
            r_2d_especial_AFT = AFT_channel_output(N*M, 1, taps, delay_taps, Doppler_taps, chan_coef,sigma_2_2d_especial_AFT(iesn0),s_2d_Especial_AFT.'); % OTFS
            r_2d_especial_AFT = r_2d_especial_AFT.';
        end
        if enable_AFT
            % AFT
            r_AFT = AFT_channel_output(N_AFT, Num_AFT_sym, taps, delay_taps, Doppler_taps, chan_coef,sigma_2_AFT(iesn0),s_AFT); % OTFS
        end
        
        
        %% OTFS demodulation
        
        if enable_OTFS
            % MMSE
            if enable_OTFS_LMMSE
                %r_OTFS = H_OTFS_eq'*(H_OTFS_eq*H_OTFS_eq' +sigma_2_OTFS(iesn0)/sig_energy_OTFS_sqrt^2*eye(M*N))^(-1)*r_OTFS;
                %
                F_N = (1/sqrt(N))*dftmtx(N);
                H_EQ_OTFS = kron(F_N, eye(M))*H_OTFS_eq*kron(F_N', eye(M));
            end
            y_OTFS = OTFS_demodulation(N,M,r_OTFS);
        end
        if enable_2d_AFT
            %if enable_OTFS_LMMSE
                %r_2d_AFT = H_2d_AFT'*(H_2d_AFT*H_2d_AFT' +sigma_2_OTFS(iesn0)/sig_energy_OTFS_sqrt^2*eye(M*N))^(-1)*r_2d_AFT;
                %
                H_EQ_2d_AFT = H_EQ_2d_AFT_Calc(N,M, c1, c2, c3, c4, H_2d_AFT);
            %end
            y_2d_AFT = AFT_2d_demodulation(N,M,r_2d_AFT, c1, c2, c3, c4);
        end
        if enable_Especial_2d_AFT
            y_2d_especial_AFT =  AFT_Especial_2d_demodulation(N,M,r_2d_especial_AFT, c0, c1);
        end
        if enable_AFT
            y_AFT = AFT_demodulation(N_AFT,Num_AFT_sym, c0, c1, c2,r_AFT);
        end
        
        
        %% detector
        
        if enable_OTFS
            if enable_OTFS_LMMSE
                %x_est_OTFS = y_OTFS;
                %
                x_est_OTFS  = H_EQ_OTFS'*(H_EQ_OTFS*H_EQ_OTFS' +sigma_2_OTFS(iesn0)/sig_energy_OTFS_sqrt^2*eye(M*N))^(-1)*y_OTFS(:);
            else
                x_est_OTFS = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2_OTFS(iesn0),y_OTFS.');
                x_est_OTFS = x_est_OTFS.';
            end
        end
        if enable_2d_AFT
            %if enable_OTFS_LMMSE
                %x_est_2d_AFT = y_2d_AFT;
                x_est_2d_AFT = H_EQ_2d_AFT'*(H_EQ_2d_AFT*H_EQ_2d_AFT' +sigma_2_OTFS(iesn0)/sig_energy_OTFS_sqrt^2*eye(M*N))^(-1)*y_2d_AFT(:);
%             else
%                 x_est_2d_AFT = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2_OTFS(iesn0),y_2d_AFT.');
%                 x_est_2d_AFT = x_est_2d_AFT.';
%             end
        end
        if enable_Especial_2d_AFT
            x_est_2d_especial_AFT = especial_2d_AFT_mp_detector(M, N, c0, c1, taps,delay_taps,Doppler_taps,chan_coef,...
                y_2d_especial_AFT, sigma_2_2d_especial_AFT(iesn0), sig_energy_2d_especial_AFT_sqrt);
        end
        if enable_AFT
            if enable_ML
                x_est_AFT = AFT_ML_detector(N_AFT, Num_AFT_sym, c0, c1, c2,taps,delay_taps,Doppler_taps,chan_coef, y_AFT);
            else
                x_est_AFT = AFT_mp_detector(N_AFT, Num_AFT_sym, c0, c1, c2,taps,delay_taps,Doppler_taps,chan_coef, y_AFT);
            end
            if SC_AFT
                x_est_AFT = transpose(ifft(transpose(x_est_AFT)));
                N_AFT = N_AFT/SC_fac;
                x_est_AFT = x_est_AFT(:, 1:N_AFT);
            end
        end
        %% output bits and errors count
        
        if enable_OTFS
            % OTFS
            %data_demapping = qamdemod(x_est_OTFS,M_mod,0, 'gray');
            data_demapping = qamdemod(x_est_OTFS,M_mod,'gray');
            data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
            errors = sum(xor(data_info_est,data_info_bit));
            err_ber_OTFS(iesn0) = errors + err_ber_OTFS(iesn0);
        end
        if enable_2d_AFT
            % OTFS
            %data_demapping = qamdemod(x_est_OTFS,M_mod,0, 'gray');
            data_demapping = qamdemod(x_est_2d_AFT,M_mod,'gray');
            data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
            errors = sum(xor(data_info_est,data_info_bit));
            err_ber_2d_AFT(iesn0) = errors + err_ber_2d_AFT(iesn0);
        end
        if enable_Especial_2d_AFT
            data_demapping = qamdemod(x_est_2d_especial_AFT,M_mod,'gray');
            data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
            errors = sum(xor(data_info_est,data_info_bit));
            err_ber_2d_especial_AFT(iesn0) = errors + err_ber_2d_especial_AFT(iesn0);
        end
        if enable_AFT
            % AFT
            x_est_AFT_serial           = reshape(transpose(x_est_AFT) ,[1,size(x_est_AFT,1)*size(x_est_AFT,2)]);
            %data_demapping = qamdemod(x_est_AFT_serial, M_mod,0, 'gray');
            data_demapping = qamdemod(x_est_AFT_serial, M_mod, 'gray');
            data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
            errors = sum(xor(data_info_est,data_info_bit));
            err_ber_AFT(iesn0) = errors + err_ber_AFT(iesn0);
        end
        if mod(ifram, 100) == 0
            ifram
        end
    end
    if iesn0 ==0
        sig_energy_OTFS_sqrt = sqrt(sig_energy_OTFS/N_fram);
        sig_energy_2d_especial_AFT_sqrt = sqrt(sig_energy_2d_especial_AFT/N_fram);
        sig_energy_AFT_sqrt = sqrt(sig_energy_AFT/N_fram);
        sigma_2_OTFS = abs(sig_energy_OTFS_sqrt*noise_var_sqrt).^2;
        sigma_2_2d_especial_AFT = abs(sig_energy_2d_especial_AFT_sqrt*noise_var_sqrt).^2;
        sigma_2_AFT = abs(sig_energy_AFT_sqrt*noise_var_sqrt).^2;
    end
end
if enable_OTFS
    err_ber_fram_OTFS = err_ber_OTFS/N_bits_perfram./N_fram
    semilogy(SNR_dB, err_ber_fram_OTFS,'-*','LineWidth',2);
    title(sprintf(['N = ' num2str(N) ', M = ' num2str(M) ', ' num2str(M_mod) 'QAM']))
    legend('OTFS')
    ylabel('BER'); xlabel('SNR in dB');grid on
    hold on
end
if enable_2d_AFT
    err_ber_fram_2d_AFT = err_ber_2d_AFT/N_bits_perfram./N_fram
    semilogy(SNR_dB, err_ber_fram_2d_AFT,'-*','LineWidth',2);
    title(sprintf(['N = ' num2str(N) ', M = ' num2str(M) ', ' num2str(M_mod) 'QAM']))
    legend('2d AFT')
    ylabel('BER'); xlabel('SNR in dB');grid on
    hold on
end
if enable_Especial_2d_AFT
    err_ber_fram_2d_especial_AFT = err_ber_2d_especial_AFT/N_bits_perfram./N_fram
    semilogy(SNR_dB, err_ber_fram_2d_especial_AFT,'-*','LineWidth',2);
    title(sprintf(['N = ' num2str(N) ', M = ' num2str(M) ', ' num2str(M_mod) 'QAM']))
    legend('2d AFT')
    ylabel('BER'); xlabel('SNR in dB');grid on
    hold on
end
if enable_AFT
    err_ber_fram_AFT = err_ber_AFT/N_bits_perfram./N_fram
    semilogy(SNR_dB, err_ber_fram_AFT,'-*','LineWidth',2);
    legend('1d AFT')
end
legend('OTFS', '2d OFDM','2d especial AFT' ,'1d AFT')
% if enable_OTFS_LMMSE
%     legend('OTFS MMSE', 'AFT');
% else
%     legend('OTFS Message Passing', 'AFT');
% end
toc
