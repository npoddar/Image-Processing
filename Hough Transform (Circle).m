
%%
%Prewitt Operator
%Circle%

clear all

lenna1 = imread('circle.png');
lenna = rgb2gray(lenna1);
[row1, col1] = size(lenna);

pOperator = lenna;
for i=1:row1-2
    for j=1:col1-2
        G_x = lenna(i+2, j) + 1*lenna(i+2, j+1) + lenna(i+2, j+2) - lenna(i, j) - 1*lenna(i, j+1) - lenna(i, j+2) ;
        G_y = lenna(i, j+2) + 1*lenna(i+1, j+2) + lenna(i+2, j+2) - lenna(i, j) - 1*lenna(i+1, j) - lenna(i+2, j) ;
        
        G_x = abs(G_x);
        G_y = abs(G_y);

        %pOperator(i+1, j+1) = sqrt(double(G_x*G_x+G_y*G_y));
        pOperator(i+1, j+1) = (G_x) + (G_y);
     
       end
end
i = find ( pOperator <= 100 );
j = find ( pOperator > 100 );

colormap(gray(256));
image(pOperator);

%%
%Circle

edgeImage = pOperator;
[yIndexVal xIndexVal] = find(pOperator > 0);

max_x = max(xIndexVal);
max_y = max(yIndexVal);

min_x = min(xIndexVal);
min_y = min(yIndexVal);

N = 15; %spacing
M = 1;
accumulator_circle = zeros(ceil((max_x-min_x)/N)+5, ceil((max_y-min_y)/N)+5, ceil(max_y)+5);

len_yIndex = size(yIndexVal, 1);

for c1=1:N:max_x
    for c2=1:N:max_y
        for coor=1:len_yIndex
            y = yIndexVal(coor);
            x = xIndexVal(coor);
            radius = sqrt((y - c2)^2 + (x - c1)^2);
            radius_i = ceil((radius));
            if(radius_i == 0)
                radius_i = 1;
            elseif (radius_i > max_x)
                      continue;
            else
                %disp('Good Radius');
            end
     
             if((c1>max_x) || (c2>max_y))
                 disp('entered here');
             end
      
             c1_cell = ceil((c1-min_x)/N)+1;
             c2_cell = ceil((c2-min_y)/N)+1;
                          
             accumulator_circle(c1_cell, c2_cell, radius_i) = accumulator_circle(c1_cell, c2_cell, radius_i) + 1;
        end
    end
end

%%
%Parameters of Circle
max = 0;

for c1=1:(size((accumulator_circle), 1))
    for c2=1:(size((accumulator_circle), 2))
        for rad=1:(size((accumulator_circle), 3))
            if ((accumulator_circle(c1,c2, rad) > max))
                    max=accumulator_circle(c1,c2, rad);
                    max_c1 = c1;
                    max_c2 = c2;
                    max_rad = rad;
            end
        end
    end
end

output_c1 = (max_c2-1)*N + min_x ;
output_c2 = (max_c1-1)*N + min_x ;

disp(output_c1);
disp(output_c2);
disp(max_rad);

%% Draw Circle
%clear all;

parameter = [1006, 1006, 849];

[yIndexVal xIndexVal] = find(pOperator > 0);
max_x = max(xIndexVal);
max_y = max(yIndexVal);

X = 0:1:max_x;
%Y = (sqrt((849*849) - ((X-1006)*(X-1006)))) + 1006;

Y = zeros(size(X, 2));
count = 1;

for x=0:max_x
    Y(count) = (sqrt((849*849) - ((x-1006)*(x-1006)))) + 1006;
    count = count + 1;
end

plot(X, Y); hold on;

count =1;
for x=0:max_x
    Y(count) = -(sqrt((849*849) - ((x-1006)*(x-1006)))) + 1006;
    count = count + 1;
end

plot(X, Y); hold off;
