clc
close all

A = [];
N = 8;
for i = 0:N-1
    L2 = 4;
    open(['Diversity_N8_L_0_' num2str(L2) '_F_1_' num2str(i) '.fig']);
    %open(['Diversity_N8_L_0_' num2str(L2) '_F_1_' num2str(i) '_C2_sqrt_2.fig']);
    a = get(gca,'Children');
    xdata = get(a, 'XData');
    ydata = get(a, 'YData');
    close
    figure(1)
    semilogy(xdata, ydata,'-*','LineWidth',2);
    hold on
    A = [A; ['L 0-' num2str(L2) ', F 1-' num2str(i) ', diff ' num2str(abs(mod(L2+i, N) - 1))]];
end
legend(A)
xlabel('SNR')
ylabel('BER')
title('N = 8, P = 2, BPSK, ML detection')
%print  Diversity_N8_L_0_4_F_1_All_c2_1_2N -depsc2