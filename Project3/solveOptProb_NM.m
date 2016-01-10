function [optZ] = solveOptProb_NM(init_Z,Lambda,t,X,Y,Epsilon)
% Compute the optimal solution using Newton method
%
% INPUTS:
%   costFcn: Function handle of F(Z)
%   init_Z: Initial value of Z
%   tol: Tolerance
%
% OUTPUTS:
%   optSolution: Optimal soultion
%   err: Errorr
%
% @ 2011 Kiho Kwak -- kkwak@andrew.cmu.edu

Z = init_Z;
Sigma = 1;

first=1;
while (Sigma/2) > Epsilon
    if (first==1)
        first=0;
    else
        Z=backtrackingLineSearch(Z,delta_Z,X,Y);
    end
    [G, H] = costFcn(Lambda,t,Z,X,Y);
    delta_Z=-inv(H)*G;
    Sigma=-(G'*delta_Z);
end

optZ=Z;

end
