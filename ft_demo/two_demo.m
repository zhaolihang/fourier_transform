xn = [0, 1, 2, 3; 4, 5, 6, 7; 8, 9, 10, 11; 12, 13, 14, 15];
Nx=4; Ny=4;
f = fft2(xn, Nx,Ny);
figure(1);
imshow(xn, [], 'InitialMagnification', 'fit');
figure(2);
F2=log(abs(f)+1);
imshow(F2, [], 'InitialMagnification', 'fit');