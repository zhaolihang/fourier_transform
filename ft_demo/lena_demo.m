image_path = 'lena-gray-512.bmp';
image_info = imfinfo(image_path);
img = imread(image_path);
img = double(img);
img = img/255;

figure(1);
title('ԭͼ');
imshow(img, 'InitialMagnification', 'fit');%��ʾͼ��


figure(2);
title('Ƶ��');

f=fft2(img);
f=fftshift(f);
fimage=log(abs(f)+1);   %ȡģ����������
imshow(fimage,[], 'InitialMagnification', 'fit');%��ʾͼ��

% ��ʱ��Ҫ��ʾͼ������Ҫ����absȡ���������ģ���ٽ�����ʾ��
% ȡģ��ͼ��������ֵһ���ܴ�ֱ����imshow�������޷���ʾ�ģ�
% ��ʱ������log����ȡ���������log(abs(F)+1)�������Ϳ��Զ�Ƶ�׽������š�
% ����Ϊʲô��+1��log(x)���ڣ�0��1��֮���xֵ����ȡ��������ɸ�ֵ��
% ��log��x+1�������е�xֵӳ�䵽������Χ�ڡ�
% ԭ�����ӣ�https://blog.csdn.net/weixin_43637490/article/details/89196212
