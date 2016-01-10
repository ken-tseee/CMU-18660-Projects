function alpha = OMP(A, B, lambda)
F = B;
n = size(A,2);
omega = zeros(n,1);
alpha = zeros(n,1);
sizeOmega = size(omega,1);

for p=1:lambda
    theta = abs(A' * F);
    [~, indexTheta] = max(theta);
    omega(indexTheta,:) = indexTheta;
    temp_A = [];
    for i=1:sizeOmega
        if(omega(i,1)~=0)
            temp_A = [temp_A A(:,i)];
        end
    end
    temp_alpha = temp_A \ B;
    
    F = B;
    idx_ta = 0;
    for i=1:sizeOmega
       if(omega(i,1)~=0)
           idx_ta = idx_ta + 1;
           alpha(i,1) = temp_alpha(idx_ta,1);
           F = F - alpha(i,1) * A(:,i);
       else
           alpha(i,1) = 0;
       end
    end   
end
end
