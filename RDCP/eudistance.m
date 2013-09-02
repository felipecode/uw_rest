function [J] = eudistance(I,A)
    [y x ~] = size(I);
    
    J = zeros(y, x);


    for m = 1:y
        for n = 1:x

             %eudist = sqrt((A(1) - I(m,n,1))^2 + (A(2) - I(m,n,2))^2 + (A(3) - I(m,n,3))^2);
             eudist = (A(2) - I(m,n,2))^2;
             J(m,n) = eudist;

             %J(m,n) = max([J(m,n) R(m,n)]);
        end
    end
    
end



