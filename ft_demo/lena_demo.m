image_path = 'lena-gray-512.bmp';
image_info = imfinfo(image_path);
img = imread(image_path);
img = double(img);
img = img/255;

figure(1);
title('原图');
imshow(img, 'InitialMagnification', 'fit');%显示图像


figure(2);
title('频谱');

f=fft2(img);
f=fftshift(f);
fimage=log(abs(f)+1);   %取模并进行缩放
imshow(fimage,[], 'InitialMagnification', 'fit');%显示图像

% 此时如要显示图像则需要先用abs取复数矩阵的模，再进行显示。
% 取模后图像矩阵的数值一般会很大，直接用imshow函数是无法显示的，
% 此时可以用log函数取其对数，如log(abs(F)+1)，这样就可以对频谱进行缩放。
% 至于为什么用+1，log(x)对于（0，1）之间的x值经过取对数后会变成负值，
% 而log（x+1）则将所有的x值映射到正数范围内。
% 原文链接：https://blog.csdn.net/weixin_43637490/article/details/89196212
