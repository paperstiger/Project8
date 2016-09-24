%% The main function to process a given signal with given 
%  Input signal, phase function, and the bandwidth K
%  Output: shape function s and envelope a(t)
function main
clear;
clc;
%% Generate the signal
N = 2^12;
t = 0:1/N:1-1/N;
y = SignalFun(t');
% plot(t, y)
%% Interpolate
% Get phi at different points
Vtheta = PhiFun(t);
theta0 = PhiFun(0);
thetaf = PhiFun(1);
Ltha = (thetaf - theta0)/2/pi;
theta_bar = (0:N-1) / N *2*pi;
theta = Ltha * theta_bar + theta0;
Vftheta = spline(Vtheta, y, theta);
%Then scale Vthetabar to [0, 2 * pi]
Vthetabar = 1/N*(0:N-1);
Vthetabar = Vthetabar';
%% Compute the FFT of f in theta-space
Ff = fftshift(fft(Vftheta));
%% Compute f_{theta, k}
K = 100;
MatF = zeros(N, 2*K + 1);
for k = 0:K
    w0 = (k - 1/2) * Ltha;
    wf = (k + 1/2) * Ltha - 1;
    ind = w0 + 1 + N/2: wf + 1 + N/2;
    tmpary = zeros(N, 1);
    tmpary(N/2 + 1 - Ltha/2 : N/2 + Ltha/2) = Ff(ind);
    f_theta_k = ifft(ifftshift(tmpary));
    if k == 0
        MatF(:, 1) = real(f_theta_k);
    else
        MatF(:, 1 + k) = real(f_theta_k);
        MatF(:, K + 1 + k) = imag(f_theta_k);
    end
end
%% Use SVD to calculate the optimal solution
[U,S,V] = svd(MatF);
a = U(:, 1);
c = V(:, 1);
fprintf('The second SV is %f\n', S(2, 2));
% c = V(1, :)';
%% Compute the shape function and envelope function
% compute shape function using ifft
% Calculate s(tl)
stl = zeros(N, 1);
for i = 1 : N
    stl(i) = c(1) + 2*sum(c(2:K+1)'.*cos((1:K)*2*pi*Vthetabar(i)) - c(K+2:2*K+1)'.*sin((1:K)*2*pi*Vthetabar(i)));
end
InterpTheta = (Vtheta - theta0)/(thetaf - theta0) * (N - 1) / N;
atl = spline(Vthetabar, a, InterpTheta');
%Normalize
[val, ind] = max(abs(stl));
if stl(ind) < 0
    val = -val;
end
stl = stl / val;
%Compare the shape function
figure(1)
subplot(2, 1, 1)
plot(2*pi*Vthetabar, stl)
title('Shape Function')
subplot(2, 1, 2)
tmpy = 1./(1.1 + cos(theta_bar + cos(2*theta_bar)));
tmpy = tmpy/(max(tmpy));
plot(theta_bar, tmpy);
%ezplot('1/(1.1 + cos(phi + cos(2*phi)))', [0, 2*pi])
%% Rebuild the signal
%Repeat the stl for Ltha times
ShpFun = repmat(stl, Ltha, 1);
ShpFun = ShpFun(1 : Ltha : end);
Sig = val * atl .* ShpFun * S(1, 1);
figure(2)
subplot(2, 1, 1)
plot(t, Sig)
title('The Signal')
subplot(2, 1, 2)
plot(t, y)
end

%% The phase function, return \phi(t), vector form
function phi = PhiFun(t)
phi = 40 * pi * t + 2 * cos(2*pi*t);
end
%% The function which return signals with noises
function y = SignalFun(t)
a = 1./(2 + sin(2*pi*t));
% a = 2;% Use a constant in order to simplify
phi = PhiFun(t);
y = a./(1.1 + cos(phi + cos(2*phi))) + 0.0*randn(length(t), 1);
%y = a .* (cos(phi)+sin(2*phi));
% y = a .* sin(3*phi);
end