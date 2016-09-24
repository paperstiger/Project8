clear;
close all;

% This code solves the general mode decomposition problem under the weak
% well-separation and well-different condition in the paper "Synchrosqueezed
% wave packet transforms and diffeomorphism based spectral analysis for 1D
% general mode decompositions", H. Yang, Applied and Computational Harmonic
% Analysis, 2014.
%
% Test spike wave shape functions
%
% By Haizhao Yang


test = 2; % test = 1, real example by Hou and Shi's paper, PROBLEM: need to debug
% test = 2, complex example by our intuition, not bad



if(1)
    %set up data
    N = 2048/4;
    x = (0:N-1)/N;
    w = (-N/2:1:(N/2-1));
    F2 = 20;
    amp = 0.005;
    yy = x + amp*cos(2*pi*x);
    switch test
        case 1
            f = ( cos(2*pi*F2*yy) + cos(2*pi*2*F2*yy) )/2;%gen_shape(F2*yy,2);
            shapeTrue = cos(2*pi*x)+cos(2*pi*2*x);
            shapeTrue = shapeTrue(:);
            [val,pos] = max(shapeTrue);
            shapeTrue = circshift(shapeTrue,-pos);
            shapeTrue = shapeTrue/(norm(shapeTrue)/sqrt(length(shapeTrue)));
        case 2
            % Note: test different cases
            %s = @(t) exp(2*pi*i*t)+0.3*exp(2*pi*i*2*t);%
            %s = @(t) exp(2*pi*i*t)+0.3*exp(2*pi*i*2*t) + 0.4*exp(2*pi*i*3*t);
            s = @(t) exp(2*pi*1i*t)+0.3*exp(2*pi*1i*2*t) + 0.4*exp(2*pi*1i*3*t) + 0.2*exp(2*pi*1i*4*t);
            %s = @(t) exp(2*pi*i*t)+0.3*exp(2*pi*i*2*t) + 0.4*exp(2*pi*i*3*t) + 0.2*exp(2*pi*i*4*t) + 0.05*exp(2*pi*i*5*t);
            f = s(F2*yy);%gen_shape(F2*yy,2);
            shapeTrueReal = real(s(x));
            shapeTruReale = shapeTrueReal(:);
            [val,pos] = max(shapeTrueReal);
            shapeTrueReal = circshift(shapeTrueReal,-pos);
            shapeTrueReal = shapeTrueReal/(norm(shapeTrueReal)/sqrt(length(shapeTrueReal)));
            shapeTrueImag = imag(s(x));
            shapeTrueImag = shapeTrueImag(:);
            [val,pos] = max(shapeTrueImag);
            shapeTrueImag = circshift(shapeTrueImag,-pos);
            shapeTrueImag = shapeTrueImag/(norm(shapeTrueImag)/sqrt(length(shapeTrueImag)));
    end
    
    
    freq = F2*(1-2*pi*amp*sin(2*pi*x) );
    phase = F2*yy;
end


if (1)
    % step 1: interpolation
    K = 100;
    uniform_sample = x;
    nonuniform_sample = phase/F2;
    switch test
        case 1
            f_theta = spline(nonuniform_sample,f,uniform_sample);
        case 2
            f_theta = spline(nonuniform_sample,real(f),uniform_sample);
            f_theta = f_theta + sqrt(-1)*spline(nonuniform_sample,imag(f),uniform_sample);
    end
    % step 2: FFT
    f_theta_hat = fftshift(fft((f_theta)));
    figure;plot(w,abs(f_theta_hat));
    % step 3: pad zero and ifft
    fh = f_theta_hat';
    Fh = zeros(N,K+1);% repmat(fh,1,K+1);
    Ltha = F2;
    for k = 1:K+1
        st = (N/2 + 1-Ltha/2)+(k-1)*Ltha; ed = (N/2 + 1 + Ltha/2 -1)+(k-1)*Ltha;
        if ed<=N
            Fh(st:ed, k) = fh(st:ed);
        end
    end
    f_theta_k = zeros(N,K+1);
    for cnt = 1:K+1
        f_theta_k(:,cnt) = ifft(ifftshift(Fh(:,cnt)));
    end
    switch test
        case 1
            % Get the huge matrix
            Ftheta = [real(f_theta_k(:, 1)), real(f_theta_k(:, 2 : K + 1)), imag(f_theta_k(:, 2 : K + 1))];
            [U, S, V] = svd(Ftheta);
            s1 = U(:, 1); % problem: we should get a constant vector
            c = V(:,1)'; % problem: we should get two nonzero entries
            ctrue = zeros(1,2*K+1);
            ctrue(K+1) = c(1);
            for cnt = 1:K
                ctrue(cnt+K+1) = c(1+cnt)+sqrt(-1)*c(1+cnt+K);
                ctrue(2*K+1-(cnt+K)) = ctrue(cnt+K+1)';
            end
            ctruetrue = zeros(N,1);
            ctruetrue(N/2-K+1:N/2+K+1) = ctrue;
            %shape = real(fft(fftshift(ctruetrue)));
            shape = ifft(ifftshift(ctruetrue));
            [val,pos] = max(shape);
            shape = circshift(shape,-pos);
            shape = shape/(norm(shape)/sqrt(length(shape)));
            figure;subplot(1,2,1);plot(shapeTrue);title('true shape');subplot(1,2,2);plot(real(shape));title('recovered shape');
        case 2
            eps = 1e-1;
            [U,S,V] = svd(f_theta_k);
            pos = find(diag(abs(S))/abs(S(1,1))>eps);
            numComp = length(pos);
            S = diag(U(1,:))*S; % Note: need to use U(1,:) to adjust the initial phase
            V = S*V'; % Note: need to combine S and V, Hou and Shi forgot this scalar in their paper
            c = zeros(1,K+1);
            for cnt = 1:numComp % Note: Hou and Shi might make a mistake in the derivation, we need multiple components after SVD
                c = V(cnt,:)+c;
            end
            ctrue = zeros(1,2*K+1);
            ctrue(K+1) = c(1);
            for cnt = 1:K
                ctrue(cnt+K+1) = c(1+cnt);
            end
            ctruetrue = zeros(N,1);
            ctruetrue(N/2-K+1:N/2+K+1) = ctrue;
            shapeReal = real(ifft(fftshift(ctruetrue)));
            [val,pos] = max(shapeReal);
            shapeReal = circshift(shapeReal,-pos);
            shapeReal = shapeReal/(norm(shapeReal)/sqrt(length(shapeReal)));
            
            shapeImag = imag(ifft(fftshift(ctruetrue)));
            [val,pos] = max(shapeImag);
            shapeImag = circshift(shapeImag,-pos);
            shapeImag = shapeImag/(norm(shapeImag)/sqrt(length(shapeImag)));
            figure;subplot(2,2,1);plot(shapeTrueReal);title('true shapeReal');subplot(2,2,2);plot(real(shapeReal));title('recovered shapeReal');
            subplot(2,2,3);plot(shapeTrueImag);title('true shapeImag');subplot(2,2,4);plot(real(shapeImag));title('recovered shapeImage');
        case 3
            % Get the huge matrix
            Ftheta = [real(f_theta_k(:, 1)), real(f_theta_k(:, 2 : K + 1)), imag(f_theta_k(:, 2 : K + 1))];
    end
end

