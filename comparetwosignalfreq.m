%% See the two shape function if have similar distribution of frequency
%If so, the weird results can be explained, otherwise
N = 2^10;
phi = (0:N-1)/N * 2*pi;
% y1 = 1./(1.1 + cos(phi + cos(2*phi)));
% y2 = 1./(2 + sin(3*phi));
y1 = cos(phi);
y2 = sin(3*phi);
f1 = fftshift(fft(y1));
f2 = fftshift(fft(y2));
figure(250)
plot(abs(f1))
hold on
plot(abs(f2))