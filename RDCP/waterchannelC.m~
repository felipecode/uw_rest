function [J] = waterchannelC(I, patch_size,A)
    [y x ~] = size(I);
    
    J = zeros(y, x);

    
    % add some space around the edges
    I = padarray(I, [floor(patch_size/2) floor(patch_size/2)], 'symmetric');
    for m = 1:y
        for n = 1:x
             patch = I(m:(m+patch_size-1), n:(n+patch_size-1),:);
             J(m,n) = min(min(eudistance(patch,A))); 

             %J(m,n) = max([J(m,n) R(m,n)]);
        end
    end
    
end






