xn=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; %����Ĳ�������
N=16;  %���г���
f=fft(xn,N);
L=0:1:N-1;  %���᳤��
subplot(1,1,1);  %��ͼ
stem(L,abs(f));  %����Դ
xlabel('k');  %���ú�������
ylabel('X(k)');   %������������
title('DFT N=16');