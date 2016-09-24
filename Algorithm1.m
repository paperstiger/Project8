clear all;close all;clc;
%9/9/2016 Fei Tan
%signal
N = 2^10;
t_l = 0:1/N:(1-1/N);   %make range theta_tl,theta similar??
theta_tl = 40 * pi * t_l +2 * cos(6 * pi * t_l);
a_tl = 1 ./ (2 + sin(2 * pi * t_l));
f_tl = a_tl .* (sin(theta_tl/2) + cos(theta_tl/4));

%constants
L_theta = (40 * pi * 1 + 2 * cos(6 * pi * 1) - theta_tl(1)) / (2 * pi);
K = 20;%round(L_theta/2);% K<N/L_theta/2

theta_bar = (0:N-1) / N *2*pi;
theta = L_theta * theta_bar + theta_tl(1);
f_theta = spline(theta_tl,f_tl,theta);

%fourier coefficients
omega = -N/2:N/2-1;
%f_theta_omega = zeros(1,N);
%for l = 1:N
%f_theta_omega(l) = sum(f_theta .* exp((-1i) * omega(l) * theta_bar));
%end
f_theta_omega = f_theta * exp((-1i) * theta_bar' * omega);


%compute f_theta_k
%以下存疑 待优化
f_theta_k = zeros(K+1,N);
%for t = 1:N
%    for k = [0:K]+1
%        index_omega = find(omega >= ((k - 0.5) * L_theta) & omega <= ((k + 0.5) * L_theta - 1));
%        f_theta_k(k,t) = sum(f_theta_omega(index_omega) .* exp((1i) * omega(index_omega) * theta_bar(t)));
%    end
%end
index = find(omega >= (-1/2) * L_theta & omega < (K + 1/2) * L_theta);
temp_mat = f_theta_omega' * ones(1,length(index)) .* exp((1i) * theta_bar' * omega(index)) ;
for k = 0:K
    f_theta_k(k+1,:) = sum(temp_mat(:,(k * L_theta + 1): ((k + 1) * L_theta)),2);
end
%以上存疑 待优化

%Assemble F_theta
F_theta_real = real(f_theta_k);
F_theta_imag = imag(f_theta_k);
F_theta = [F_theta_real;F_theta_imag(2:end,:)]';

%SVD
[U,S,V] = svd(F_theta);
a_theta = U(:,1);
c = V(:,1)';
c_k = c(1:K+1)+[0,c(K+2:end)] * 1i;
c_k_minus = c(2:K+1)+ c(K+2:end) * -1i;
c_total = [c_k_minus,c_k];
%shape function
%tl = theta_bar;
k = -K:K;
%for l = 1:N
%s(l) = sum(c_total .* exp(1i * k * tl(l) ) );
%end
s_tl = c_total * exp(1i * k' * theta_bar);
plot(theta_bar, s_tl)
a_tl_cal = interp1(theta,a_theta,theta_tl);