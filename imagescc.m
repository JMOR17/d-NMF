function [] = imagescc(c)
    d = size(c,1);
    imagesc(reshape(c,[sqrt(d) sqrt(d)]));
end