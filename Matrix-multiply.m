function m = multiply(image, mask, size)
    sum = 0;
    for i= 1:size;
        for j=1:size;
            sum = image(i,j) * mask(i,j) + sum;
        end
    end
    m = abs(sum);
end
    
