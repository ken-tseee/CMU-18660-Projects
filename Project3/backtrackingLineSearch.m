function [Z] = backtrackingLineSearch(Z,delta_Z,X,Y)
s=1;
numFeature=size(X,1);

while (true)
    tempZ=Z+s*delta_Z;
    tempW=tempZ(1:numFeature,1);
    tempC=tempZ(numFeature+1,1);
    tempXi=tempZ(numFeature+2:end,1);
    temp=tempW'*X.*Y+tempC*Y+tempXi'-1;
    if ((isempty(find(temp<=0,1))) && (isempty(find(tempXi<=0,1))))
        break;
    else
        s=0.5*s;
    end
end

Z=Z+s*delta_Z;

end
