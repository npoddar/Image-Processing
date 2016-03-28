%%
%Roberts Operator

clear all
load lenna.mat

rOperator = lenna;
for i=1:255
    for j=1:255
        imageMat = lenna(i:i+1, j:j+1); 
        maskx = [+1 0; 0 -1];
        masky = [0, 1; -1 0];
        fx = multiply(imageMat, maskx, 2);
        fy = multiply(imageMat, masky, 2);
        rOperator(i, j) = sqrt(fx*fx+fy*fy);
     end
end
colormap(gray(256));
image(rOperator);

%%
%Sobel Operator
clear all
load lenna.mat

sOperator = lenna;
for i=1:254
    for j=1:254
        imageMat = lenna(i:i+2, j:j+2); 
        maskx = [-1 0 1; -2 0 2; -1 0 1];
        masky = [-1 -2 -1; 0 0 0; 1 2 1];
        fx = multiply(imageMat, maskx, 3);
        fy = multiply(imageMat, masky, 3);
        sOperator(i+1, j+1) = sqrt(fx*fx+fy*fy);
     end
end
colormap(gray(256));
image(sOperator);

%%
%Prewitt Operator
clear all
load lenna.mat

pOperator = lenna;
for i=1:254
    for j=1:254
        imageMat = lenna(i:i+2, j:j+2); 
        maskx = [-1 0 1; -1 0 1; -1 0 1];
        masky = [-1 -1 -1; 0 0 0; 1 1 1];
        fx = multiply(imageMat, maskx, 3);
        fy = multiply(imageMat, masky, 3);
        pOperator(i+1, j+1) = sqrt(fx*fx+fy*fy);
     end
end
colormap(gray(256));
image(pOperator);


%%
%Laplace of Gaussian Operator
clear all
load lenna.mat

lOperator = lenna;
for i=1:251
    for j=1:251
        imageMat = lenna(i:i+4, j:j+4); 
        maskx = [0 0 -1 0 0;0 -1 -2 -1 0; -1 -2 16 -2 -1; 0 -1 -2 -1 0; 0 0 -1 0 0];
        %masky = [-1 -1 -1; 0 0 0; 1 1 1];
        fx = multiply(imageMat, maskx, 5);
        %fy = multiply(imageMat, masky, 5);
        lOperator(i+2, j+2) = fx ;
     end
end
colormap(gray(256));
image(lOperator);

%%
%Laplacian Operator
clear all
load lenna.mat

lOperator = lenna;
for i=1:254
    for j=1:254
        imageMat = lenna(i:i+2, j:j+2); 
        maskx = [0 -1 0; -1 4 -1;0 -1 0];
        fx = multiply(imageMat, maskx, 3);
        %fy = multiply(imageMat, masky, 3);
        lOperator(i+1, j+1) = fx;
     end
end
colormap(gray(256));
image(lOperator);
