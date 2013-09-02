

% image file
image = 'Images/Water/merg.png';
scale = 0.4;

% used when generating the prior: width of the patch that is examined
% to determine the dark channel value for a point
patch_size = 15;

% denoted by lowercase omega in He, Sun, Tang paper
% fraction of haze to remove
% leaving some (removal_amount < 1.0) preserves 'aerial perspective'
removal_amount = 0.95;%0.95;

im = imresize(imread(image), scale);
addpath(genpath('RDCP'));
imshow(im);
im_rest = rdcp(im,removal_amount,patch_size);
%rmpath RDCP;
fprintf("Restarou\n");
%fprintf(im_rest);
im_rest
figure;
imshow(im_rest);
title('restaurada');
