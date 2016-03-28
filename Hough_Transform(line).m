%%
%Prewitt Operator

clear all

lenna1 = imread('Triangle-lines.png');
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
%pOperator(i) = 0;
%pOperator(j) = 255;

colormap(gray(256));
image(pOperator);

%%
%Hough Transform Line
%x-y plane

edgeImage = pOperator;
[size_y size_x] = size(edgeImage);
[yIndexVal xIndexVal] = find(pOperator > 0);

diff = 0.1;

a = -1:diff:1;
b = -size_y:1:size_y;

accumulator = zeros(size((a),2), size((b),2));

for loop=1:size(yIndexVal)
    x_val = xIndexVal(loop);
    y_val = yIndexVal(loop);
    for a_i=1:size(a,2)
        a_val = a(a_i);
        b_val = (-1 * x_val * a_val) + y_val;
        b_val = floor(b_val);
        b_i = find(b==b_val);
        accumulator(a_i, b_i) = accumulator(a_i, b_i) + 1;
    end
end

%%
%Find Local Maxima

a_mat = -1:diff:1;
b_mat = -size_y:1:size_y;

N = 10;
max_prev = 0;
for i=1:N
    max= 0;
    for i=1:size(accumulator, 1)
        for j=1:size(accumulator, 2)
            if(max_prev ~= 0)
                if ( (accumulator(i,j) > max) && (accumulator(i,j) < max_prev) )
                    max=accumulator(i,j);
                    index = [i j];
                    max_a= a_mat(i);
                    max_b= b_mat(j);
                end
            else
               if ( (accumulator(i,j) > max) )
                    max=accumulator(i,j);
                    index = [i j];
                    max_a= a_mat(i);
                    max_b= b_mat(j);
                end
            end
        end
    end
    disp('max');
    disp(max);
    disp('max-a');
    disp(max_a);
    disp('max-b');
    disp(max_b);
    disp('index');
    disp(index);
    max_prev = max;
end

%%
%Hough Transform - Parameter Space

dbstop if caught error;

edgeImage = pOperator;

max_theta = 90;
min_theta = 0;
diff_theta = 1;
cells_theta = ceil((max_theta - min_theta)/diff_theta);

[row1, col1] = size(edgeImage);

max_rho = ceil(sqrt(row1*row1 + col1*col1)) ;
min_rho = 0;
diff_rho = max_rho / cells_theta;

cells_rho = cells_theta;

accumulator = zeros(cells_theta, cells_rho);

for i = 1:row1
    for j=1:col1
        if (edgeImage(i, j) > 0)
            for theta_index=1:max_theta
                rho = i*cosd(theta_index) + j*sind(theta_index);
                rho_index = ceil(rho / diff_rho);
                accumulator(rho_index, theta_index) =  accumulator(rho_index, theta_index) + 1;
            end
        end
    end
end

theta_range = 1:1:90;
rho_range = min_rho:diff_rho:max_rho-1;
pcolor(theta_range,rho_range,accumulator);
shading flat;
xlabel('Theta');
ylabel('Rho');
colormap('gray');
%---------------------------------------%

colormap(gray(256));
%image(pOperator);


max_prev = 10;
for i=1:max_prev
    max= 0;
    for i=1:cells_theta
        for j=1:cells_rho
            if(max_prev ~= 0)
                if ( (accumulator(i,j) > max) && (accumulator(i,j) < max_prev) )
                    max=accumulator(i,j);
                    max_theta=i;
                    max_rho=j*diff_rho;
                end
            else
               if ( (accumulator(i,j) > max) )
                    max=accumulator(i,j);
                    max_theta=i;
                    max_rho=j*diff_rho;
                end
            end
        end
    end
    disp('max');
    disp(max);
    disp('max-theta');
    disp(max_theta);
    disp('max-rho');
    disp(max_rho);
    max_prev = max;
end

