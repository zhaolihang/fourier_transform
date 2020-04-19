clear
f = 1000;            % Sampling frequency
T = 1/f;             % Sampling period
L = 1000;            % Length of signal
t = (0:L-1)*T;       % Time vector  0-1

S = 0.9*sin(2*pi*50*t) + 1.8*sin(2*pi*120*t); %standard signal
X = S + 1*randn(size(t)); %there are some noise mixing in the signal.
%X = S;

%draw the curve
figure(1);
plot(t,X);            % show the original signal
xlabel('Ê±¼ä')
ylabel('Õñ·ù')


% fft

Cn = fft(X);           % Fourier Transform of the signal
Cn = Cn*2 / L;           % calculate the amplitude of Fourier Transform of the signal

figure(2);
plot(1:L/2,abs(Cn(1:L/2)));% show the Fourier Transform of the signal
xlabel('Hz');
ylabel('Õñ·ù');


Y = fft(X);           % Fourier Transform of the signal
threadhold = 150;                        %setting the filtering threadhold 
Y(threadhold:(L-threadhold)) = 0;      %filtering

Y = ifft(Y);                          %Inverse Fourier transform
figure(3);
plot(t,Y);
title('signal after filtering');
xlabel('time');
ylabel('amplitude');