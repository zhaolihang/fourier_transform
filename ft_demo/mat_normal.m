%矩阵操作
N=100;
Nx=N; Ny=N;
xn = zeros(Nx,Ny);
sz=size(xn);

index=0;
for r=1:Nx
    for c=1:Ny
        xn(r,c)=index;
        index=index+1;
    end
end

max_num = max(max(xn));
xn = xn./max_num; % 0-1化

figure(1);
imshow(xn, 'InitialMagnification', 'fit');%显示图像

f2=fft2(xn, Nx, Ny);
f2=abs(f2);
f2=log(f2+1);

f2 = f2./max(max(f2)); % 0-1化

figure(2);
imshow(f2, 'InitialMagnification', 'fit');%显示图像

figure(3);
f3 = fftshift(f2);
imshow(f3, 'InitialMagnification', 'fit');%显示图像


figure(4);
f5=ifft2(fft2(xn, Nx, Ny));
imshow(f5, 'InitialMagnification', 'fit');%显示图像


figure(5);
title('log(x) 和 log(x+1)');
% syms x;
% y=log(x);
% y1=log(x+1);
% ezplot(y);
% ezplot(y1);

% 绘图 https://jingyan.baidu.com/article/200957619f7b3ccb0721b483.html
% 蓝色：“b”、绿色：“g”、红色：“r”、蓝绿：“c”、紫红：“m”、黄色：“y”、黑色：“k”。
% 实线：“-”、点线：“:”、点划线：“-.”、虚线：“--”。
x=0:0.01:10;
y=log(x);
x1=0:0.01:10;
y1=log(x+1);
plot(x,y,'-r', x1,y1,'-g');
%hold on;
%加网格
grid on;
legend('log(x)','log(x+1)');
xlabel('x');
ylabel('y');

%figure; %figure 开启新的图片显示窗口

