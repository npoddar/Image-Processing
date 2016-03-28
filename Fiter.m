%%
clear all
load lenna.mat
colormap(gray(256))
image(lenna)

%%
x=1:256;
y=zeros(1,256);
for i=1:256;
    for j=1:256;
        y(lenna(i,j)) = y(lenna(i,j)) + 1;
    end
end
bar(x,y)
title('Intensity Histogram')
xlabel('Intensity')
ylabel('Numbers')

%%
Pb = zeros(1,256);
sum = 0;
output_intensity = 0;
sum_pixel = 0;

for i=1:256
   Pb(i) = ( (y(i) * 100.0) / (256.0*256.0) );
   sum = Pb(i) + sum;
   sum_pixel = y(i) + sum_pixel;
   if (sum >= 20)
       output_intensity = i;
       break;
   end
 end
disp(output_intensity);
disp(sum_pixel);
%output_intensity = 47;
%sum_pixel = 13293;

%%
load lenna.mat
x = randn(256);
y = 2*x + 2;
newLenna = zeros(256, 256);
for i=1:256;
    for j=1:256;
        newLenna(i,j) = lenna(i,j) + y(i,j);
    end
end
colormap(gray(256));
image(newLenna);

%%
smoothFilter = newLenna;

for i = 1 : 254;
    for j = 1 : 254;
        smoothFilter(i+1, j+1) = ((newLenna(i,j) + newLenna(i+1, j) + newLenna(i+2, j) + newLenna(i, j+1) + newLenna(i, j+2) + newLenna(i+1, j+2) + newLenna(i+2, j+2) + newLenna(i+2, j+1))/8.0);
    end
end

colormap(gray(256));
image(smoothFilter);

%%
%load mysort.m;

medianFilter = zeros(256, 256);
medianFiler = newLenna;
for i = 1:254;
    for j=1:254;
        array = zeros(1, 9);
        
        array(1) = newLenna(i, j);
        array(2) = newLenna(i, j+1);
        array(3) = newLenna(i, j+2);
        
        array(4) = newLenna(i+1, j);
        array(5) = newLenna(i+1, j+1);
        array(6) = newLenna(i+1, j+2);
        
        array(7) = newLenna(i+2, j);
        array(8) = newLenna(i+2, j+1);
        array(9) = newLenna(i+2, j+2);
        
        median= mysort(array, 9);
        medianFilter(i+1, j+1) = median(5);
    end
end
colormap(gray(256));
image(medianFilter);


    
