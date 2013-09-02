%% Configuration Parameters
close all;
clear all;

% image file
image = 'Images/truccore.png';
scale = 1;
% used when generating the prior: width of the patch that is examined
% to determine the dark channel value for a point
patch_size = 15;

% denoted by lowercase omega in He, Sun, Tang paper
% fraction of haze to remove
% leaving some (removal_amount < 1.0) preserves 'aerial perspective'
removal_amount = 0.95;%0.95;

%% Some initial setup

im = imresize(imread(image), scale);
%im = imcrop(im)
im = double(im)./255;
dims = size(im);
pixels = dims(1)*dims(2);

grayIm = rgb2gray(im);
blurO = wBlurEstimation(grayIm);
contrastO = imgContrast(im);
noiseO = wNoise(grayIm);

break;


im = im.*255;
[wR,wG,wB,imWB]=general_cc(im,0,-1,0);
imWB = imWB./255;
 figure;
 imshow(imWB);
  title('White balance Prior');
im = im./255;

%imwrite(imWB,strcat(path,'wbprior.png'));

blurWB = wBlurEstimation(rgb2gray(imWB));
contrastWB = imgContrast(imWB);
noiseWB = wNoise(rgb2gray(imWB));
  

% 
% figure;
% imshow(im(:,:,1));
% title('Imagem Original R');
% 
% figure;
% imshow(im(:,:,2));
% title('Imagem Original G');
% 
% figure;
% imshow(im(:,:,3));
% title('Imagem Original B');

%% Compute the dark channel (non-normalized)

[dark,R,G,B,invR,invG,invB] = waterchannel(im, patch_size);

figure;
imshow(dark);
title('Dark Channel');

figure;
imshow(R);
title('R Channel');

figure;
imshow(invR);
title('invR Channel');

figure;
imshow(G);
title('G Channel');

figure;
imshow(invG);
title('invG Channel');

figure;
imshow(B);
title('B Channel');

figure;
imshow(invB);
title('invB Channel');

'Aqui'

%% Calculate atmospheric light using dark channel

% fraction of pixels to consider for the atmospheric light calculation
count = ceil(.01*pixels);

[~, ix] = sort(reshape(dark, 1, []), 'descend'); % create a vector ix with the dark channel ordered
ix = ix(1:count);
flatimage = reshape(im, [], 3);

newimage = im;
for index = ix
    [ind1,ind2] = ind2sub(size(dark),index);
    newimage(ind1,ind2,:) = [1 1 0];
end
figure;
imshow(newimage);
title('Imagem Com Binf');

%intensity_sum = zeros(1,3);
indi=1;
for index = ix
    %intensity_sum = intensity_sum + flatimage(index, :);
    vetorzinho(indi)=sum(flatimage(index, :));
    indi = indi + 1;
end

%A = intensity_sum/count

[C,I] = max(vetorzinho);
A = flatimage(ix(I),:);


'Aqui2'

%% Calculate estimated transmission map

A_matrix = zeros(dims(1), dims(2), 3);
%A_matrix(:,:,1) = max(max(im(:,:,1)));
A_matrix(:,:,1) = 1;
A_matrix(:,:,2) = A(2);
A_matrix(:,:,3) = 1;


[dark,R,G,B,invR,invG,invB] = waterchannelB(im./A_matrix, patch_size);
A_matrix(:,:,1) = A(1);
A_matrix(:,:,3) = A(3);
h=figure;

%im(:,:,1) = im(:,:,1).*(A(3)/A(1));
%im(:,:,2) = im(:,:,2).*(A(3)/A(2));
%im = wbalance(im);
imshow(im);
title('Imagem arrumada');
%imwrite(im./A_matrix,'ImagemSemB.png','png');


figure;
imshow(dark);
title('Dark Channel/B');

figure;
imshow(R);
title('R Channel/B');

figure;
imshow(invR);
title('invR Channel/B');

figure;
imshow(G);
title('G Channel/B');

figure;
imshow(invG);
title('invG Channel/B');

figure;
imshow(B);
title('B Channel/B');

figure;
imshow(invB);
title('invB Channel/B');


% eqn 12 from HST paper
trans_est = 1 - removal_amount * dark;

figure;
imshow(trans_est);
title('Transmissao');

'Aqui3' 
tic
%% Matting the transmission map
% Code in this section was partly inspired by code originally from
% http://www.soe.ucsc.edu/classes/ee264/Winter10/Proj6-Code.zip

win_size = 3; % window size
win_els = win_size.^2; % number of window elements
l_els = win_els.^2; % number of elements calculated in each iteration

win_bord = floor(win_size./2); 

e = 0.000001; 

[m,n,c] = size(im);
numpix = m*n;

k = reshape(1:numpix, m, n);
U = eye(win_size);
D = eye(win_els);

num_els = l_els*(m-2*win_bord)*(n-2*win_bord); 

ind_i  = ones(1,num_els);
ind_j = ind_i;

els = zeros(1,num_els);

count = 0;
time1 = toc
'Aqui4'
tic 
for x = (1 + win_bord):(n - win_bord)
    for y = (1 + win_bord):(m - win_bord)
		
        wk = reshape(im(y-win_bord:y+win_bord,x-win_bord:x+win_bord,:), win_els, c);
            
        w_ind = reshape(k(y-win_bord:y+win_bord,x-win_bord:x+win_bord), 1, win_els);
            
        [i j] = meshgrid(w_ind, w_ind);
        
        i = reshape(i,1,l_els);
        j = reshape(j,1,l_els);
        
        ind_i((count*l_els + 1):(count*l_els+l_els)) = i;
        ind_j((count*l_els + 1):(count*l_els+l_els)) = j;

        win_mu = mean(wk)';

        win_cov = wk'*wk/win_els-win_mu*win_mu';

        dif = wk' - repmat(win_mu,1,win_els);
        
        elements = D - (1 + dif(:,1:win_els)'*inv(...
            win_cov + e./win_els.*U)*dif(:,1:win_els))...
            ./win_els;

        els((count*l_els + 1):(count*l_els+l_els)) = ...
            reshape(elements,1,l_els);
        
        count = count + 1;
    end
end
time2 = toc
'Aqui5'
tic 
L = sparse(ind_i, ind_j, els, numpix, numpix);

%% generate refined transmission map
time3 =toc
'Aqui6'
% recommended value from HST paper
lambda = .0001;
% equation 15 from HST
tic
a=trans_est(:) .* lambda;
'Aqui6_1'
b=lambda .* speye(size(L));
'Aqui6_2'
soma=L + b;
'Aqui6_3'
t = (soma) \ a;
time4 = toc
tic
'Aqui61'
t = t - min(t(:));

'Aqui62'
t = t ./ (max(t(:)));
'Aqui63'
t = t .* (max(trans_est(:)) - min(trans_est(:))) + min(trans_est(:));
'Aqui64'
t = reshape(t, size(trans_est));

'Aqui7'
time5 = toc
%% generate and show the image

t_sub = 0.30; % transmission map lower bound

dehazed = zeros(size(im));

figure;
imshow(t);
title('Transmissao apos Matting');
%imwrite(t,strcat(path,'transmatte.png'));
'Aqui8'
% equation 16 from HST
for c = 1:3
    dehazed(:,:,c) = min(1,max(0,(im(:,:,c) - A_matrix(:,:,c))./(max(t, t_sub)) + A_matrix(:,:,c)));
end
figure;
imshow(dehazed);
title('Dehazed');

%imwrite(dehazed,strcat(path,'dehazed.png'));
%for i = 1:4
%    
%end
%b =0.45;
%a =0.4;

%rTrans = (((t./exp(-b*(log(t)./(-a-b))))).^(1/8)).*exp(-b*(log(t)./-(a+b)));

%figure;
%imshow(rTrans);
%title('rTrans');
% equation 16 from HST
%dehazed(:,:,1) = (im(:,:,1) - A_matrix(:,:,1))./(max(rTrans, t_sub)) +
%A_matrix(:,:,1);
%dehazed = dehazed./A_matrix;

figure;

imshow(dehazed);
title('Dehazed Nova');


%figure
%grayEdge(:,:,1) = Dx_plus(dehazed(:,:,1));
%grayEdge(:,:,2) = Dx_plus(dehazed(:,:,2));
%grayEdge(:,:,3) = Dx_plus(dehazed(:,:,3));
%imshow(grayEdge);


mink_norm=5;    % any number between 1 and infinity
sigma=2;        % sigma 
diff_order=1;   % differentiation order (1 or 2)
dehazed = dehazed.*255;
[wR,wG,wB,out4]=general_cc(dehazed,0,-1,0);
out4 = out4./255;
figure;imshow(out4);
title('Max-RGB');
dehazed = dehazed./255;
blurD = wBlurEstimation(rgb2gray(dehazed));
contrastD = imgContrast(dehazed);
noiseD = wNoise(rgb2gray(dehazed));
%imwrite(out4,strcat(path,'dehazedwb.png'));
 figure;
 dehazedwb = RealGWbal(dehazed,dehazed);
 imshow(dehazedwb);
  title('White balance GWA');

  
blurF = wBlurEstimation(rgb2gray(out4));
contrastF = imgContrast(out4);
noiseF = wNoise(rgb2gray(out4));
edgeF = wEdgePreserv(grayIm,rgb2gray(out4),0);  
  
'Aqui9'
% perform histogram adjustment before showing the image
adjusted = applycform(out4, makecform('srgb2lab'));
adjusted(:,:,1) = adapthisteq(adjusted(:,:,1)/100)*100;
adjusted = applycform(adjusted, makecform('lab2srgb'));
%imwrite(out4,strcat(path,'adjusted.png'));
blurA = wBlurEstimation(rgb2gray(adjusted));
contrastA = imgContrast(adjusted);
noiseA = wNoise(rgb2gray(adjusted));
edgeA = wEdgePreserv(grayIm,rgb2gray(adjusted),0);  
  

'Aqui10'
% show the image after applying enhancement

figure;
imshow(adjusted);
title('Ajustado');

%image = 'Images/DSCN6533_R.png';
%im2 = imresize(imread(image), scale);
%im2 = double(im2)./255;



%for c = 1:3
%    dehazed(:,:,c) = min(1,max(0,im2(:,:,c)./(max(t, t_sub))));
%end
%figure;
%imshow(dehazed);
%title('Dehazed Recovered');


