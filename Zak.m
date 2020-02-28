clc
clear all
close all

a = 9;
K = 1/3;
L = 21;
b = 1/2;

a = 4;K=1, L=20;b= .5;


N = K*L^2;
n = 0:N-1;

x = exp(1j*2*pi*((a/2/L^2)*n.^2 + (b/L).*n));
x_mat_KL_L = reshape(x, K*L, L);
x_mat_L_KL = x_mat_KL_L.';
F_L = (1/sqrt(L))*dftmtx(L);
F_N = (1/sqrt(N))*dftmtx(N);
X_T = x_mat_KL_L*F_L;
X = F_L*x_mat_L_KL;



% Check
[kron(F_L, eye(K*L))*x(:) X_T(:)]
imtool(abs(X_T))

