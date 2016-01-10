function [ Temperature ] = thermalsimCholesky( p, mediumX, mediumY, leftBound, rightBound, topBound, bottomBound )
%THERMALSIMCHOLESKY solves the 2D steady state thermal problem using
%Cholesky factorization
%   INPUT:
%   p:  discretized power density
%   mediumX:    x-dimension of the medium
%   mediumY:    y-dimension of the medium
%	leftBound:	Temperature at the left boundary (x=0), leftBound(j) means
%	the temperature at T(0,j)
%	rightBound:	Temperature at the right boundary (x=N+1)
%	topBound:	Temperature at the top boundary (y=M+1)
%	bottomBound:	Temperature at the bottom boundary (y=0)
%
%   OUTPUT:
%   Temperature: solved thermal map

tic
[n,m]=size(p);
hx=(mediumX./n).^2;
hy=(mediumY./m).^2;

a=2*(hx+hy)./(hx*hy)*ones(m,1);
b=-1./hy*ones(m-1,1);
c=-1./hx*ones(m,1);
d=zeros(m,m);

v1=diag(a)+diag(b,1)+diag(b,-1);
v2=diag(c);
v3=d;

total=n*m;
A=zeros(total,total);

for i=1:m:total
   for j=1:m:total
       if(i==j)
           A(i:i+m-1,j:j+m-1)=v1;
       elseif((abs(i-j)==m))
           A(i:i+m-1,j:j+m-1)=v2;
       else
           A(i:i+m-1,j:j+m-1)=v3;
       end
   end
end

f0=zeros(n*m,1);
for i=1:n
    f0((i-1)*m+1:i*m,1)=p(i,:);
end

k=157;
f=f0./k;
T=zeros(n*m,1);
T(1:m,1)=leftBound./hx;
T((n-1)*m+1:n*m,1)=rightBound./hx;

for i=1:n
    for j=1:m
      index=mod(j,m+1)+m*(mod(i,n+1)-1);
      if(mod(index,m)==1)
          T(index,1)=T(index,1)+bottomBound(i,1)./hy;
      elseif(mod(index,m)==0)
          T(index,1)=T(index,1)+topBound(i,1)./hy;
      end
    end
end
B=f+T;

L=zeros(total,total);
for i=1:total
   for j=1:total
       if(i==j)
           L(i,j)=sqrt(A(i,j)-sum(L(j,1:j-1).^2));
       else
           L(i,j)=(A(i,j)-sum(L(i,1:j-1).*L(j,1:j-1)))./L(j,j);
       end
   end
end

V=zeros(total,1);
V(1,1)=B(1,1)./L(1,1);
for i=2:total
   summation=0;
   for j=1:i-1
      summation=summation+L(i,j).*V(j,1); 
   end
   V(i,1)=(B(i,1)-summation)./L(i,i);
end

t=zeros(total,1);
t(total,1)=V(total,1)./L(total,total);
L=L';
for i=total-1:-1:1
   summation=0;
   for j=total:-1:i+1
      summation=summation+L(i,j).*t(j,1);
   end
   t(i,1)=(V(i,1)-summation)./L(i,i);
end

Temperature=reshape(t,m,n);
Temperature=Temperature';
toc

X=A\B;
x=reshape(X,m,n);
x=x';
sum1=0;
sum2=0;
for i=1:n
    for j=1:m
        sum1=sum1+(x(i,j)-Temperature(i,j)).^2;
        sum2=sum2+(x(i,j)).^2;
    end
end
error=sqrt(sum1./sum2)
end
