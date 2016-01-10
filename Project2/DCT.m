function T = DCT(imgIn)
[hB,wB] = size(imgIn);
n = hB * wB;
T = zeros(n,n);
for p=1:n
    x = ceil(p / hB);
    y = p - (x - 1) * wB;
    au = sqrt(2 / hB);
    bv = sqrt(2 / wB);
    for q=1:n
        u = ceil(q / hB);
        v = q - (u - 1) * wB;
        angX = pi * (2 * x - 1) * (u - 1) / (2 * hB);
        angY = pi * (2 * y - 1) * (v - 1) / (2 * wB);
        T(p,q) = cos(angX) * cos(angY) * au * bv;
    end
end
T(:,1:wB) = T(:,1:wB) * sqrt(1 / 2);
for h=1:hB
    T(:,(h - 1) * hB + 1) = T(:,(h - 1) * hB + 1) * sqrt(1 / 2);
end
end
