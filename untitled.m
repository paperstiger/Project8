clear all;close all;clc;
%9/9/2016 Fei Tan
%signal
N = 2^11;
t_l = 0:1/N:1;   %make range theta_tl,theta similar??
theta_tl = 40 * pi * t_l +.2 * cos(6 * pi * t_l);
a_tl = ones(size(t_l));
f_tl = cos(theta_tl/6) + 2 * sin(2 * theta_tl/8);

%constants
L_theta =(theta_tl(end) - theta_tl(1))/2/pi; 
K = 30;
theta_bar = [0:N-1] / N * 2 * pi;
theta_bar_hat = (theta_tl-theta_tl(1))/L_theta;
f_theta = spline(theta_bar_hat,f_tl,theta_bar);

omega = -N/2:N/2-1;
f_theta_omega = fftshift(fft(f_theta));

f_theta_omega_rep = repmat(f_theta_omega',1,K+1);

for k = 1:K+1
w0 = (k-1-1/2)*L_theta;
wf = (k-1+1/2)*L_theta - 1;
ind = (w0 + N/2 + 1) : (wf + N/2 + 1);
mul = zeros(N,1);
mul(ind) = ones(L_theta,1);
f_theta_omega_rep(:,k) = f_theta_omega_rep(:,k) .* mul;
end
f_theta_k = ifft(ifftshift(f_theta_omega_rep));

%Assemble F_theta
F_theta_real = real(f_theta_k);
F_theta_imag = imag(f_theta_k);
F_theta = [F_theta_real,F_theta_imag(:,2:end)];

%SVD
[U,S,V] = svd(F_theta);
a_theta = U(:,1);
c = V(:,1)';

k = 1:K;
s_tl = c(1) + c(2:K+1) * cos(k' * theta_bar) + c(K+2:end) * sin(k' * theta_bar);

a_tl_cal = spline(theta_bar,a_theta,theta_bar_hat);