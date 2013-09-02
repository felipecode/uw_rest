function [J,R,G,B,invR,invG,invB] = waterchannelB(I, patch_size)
    [y x ~] = size(I);
    
    J = zeros(y, x);
    R = zeros(y, x);
    G = zeros(y, x);
    B = zeros(y, x);
    invR = zeros(y, x);
    invG = zeros(y, x);
    invB = zeros(y, x);
    
    % add some space around the edges
    I = padarray(I, [floor(patch_size/2) floor(patch_size/2)], 'symmetric');
    for m = 1:y
        for n = 1:x
             patch = I(m:(m+patch_size-1), n:(n+patch_size-1),:);
             %R(m,n)= max(max(patch(:,:,1)));
             R(m,n) = min(min(patch(:,:,1)));
             invR(m,n)=1 -min(min(patch(:,:,1)));
             G(m,n)=min(min(patch(:,:,2)));
             invG(m,n)=min(min(1-patch(:,:,2)));
             B(m,n)=min(min(patch(:,:,3)));
             invB(m,n)=min(min(1-patch(:,:,3)));
             J(m,n) = min([invR(m,n) G(m,n) ]);
             %J(m,n) = max([J(m,n) R(m,n)]);
        end
    end
    
end



