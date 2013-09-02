function im_rest = rdcp(im,removal_amount,patch_size) 

    if (nargin <1)
        error("Plz set input data");
    end
    if (nargin < 2)
        % seting default removal amount to 95 percent
        removal_amount = 0.95;
    end
    if (nargin < 3)
        % set default patch size to 15
        patch_size = 15
    end

    %% Some initial setup
    fprintf("Removing Water Haze \n");

    %im = imcrop(im)
    im = double(im)./255;
    dims = size(im);
    pixels = dims(1)*dims(2);

    %% Compute the dark channel (non-normalized)

    [dark,R,G,B,invR,invG,invB] = waterchannel(im, patch_size);



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
    %figure;
    %imshow(newimage);
    %title('Imagem Com Binf');

    %intensity_sum = zeros(1,3);
    indi=1;
    for index = ix
        %intensity_sum = intensity_sum + flatimage(index, :);
        vetorzinho(indi)=sum(flatimage(index, :));
        indi = indi + 1;
    end


    [C,I] = max(vetorzinho);
    A = flatimage(ix(I),:);
    %% A is the calculated airlight A.K.A  B inf



    %% Calculate estimated transmission map

    A_matrix = zeros(dims(1), dims(2), 3);
    %A_matrix(:,:,1) = max(max(im(:,:,1)));
    A_matrix(:,:,1) = A(1);
    A_matrix(:,:,2) = A(2);
    A_matrix(:,:,3) = A(3);


    [dark,R,G,B,invR,invG,invB] = darkchannel(im./A_matrix, patch_size);


    % Now the new part proposed by stephen.


    % eqn 12 from HST paper
    trans_est = removal_amount * dark;


    t = softmatting(trans_est,im);

    fprintf("Aplicou Soft Matting");
    %% generate and show the image

    t_sub = 0.25; % transmission map lower bound

    im_rest = zeros(size(im));

    
    % equation 16 from HST
    for c = 1:3
        im_rest(:,:,c) = min(1,max(0,(im(:,:,c) - A_matrix(:,:,c))./(max(t, t_sub)) + A_matrix(:,:,c)));
    end

end



 


