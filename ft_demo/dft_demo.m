xn=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; %输入的采样序列
N=16;  %序列长度
f=fft(xn,N);
L=0:1:N-1;  %横轴长度
subplot(1,1,1);  %画图
stem(L,abs(f));  %数据源
xlabel('k');  %设置横轴名称
ylabel('X(k)');   %设置纵轴名称
title('DFT N=16');