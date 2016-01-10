clc;
g1=[12.683683 12.831034 12.992619];
g2=[126.144976 125.544777 125.710479];
g3=[1147.385936 1141.724503 1151.625462];

c1=[12.535561 12.572393 12.712288];
c2=[98.339132 93.569628 94.963733];
c3=[782.649608 775.091424 766.032004];

ga1=mean(g1);
ga2=mean(g2);
ga3=mean(g3);

ca1=mean(c1);
ca2=mean(c2);
ca3=mean(c3);

x=[1 2 3];
y1=[ga1 ga2 ga3]
y2=[ca1 ca2 ca3]

plot(x,y1);
hold on
plot(x,y2);
grid on
title('Plot for runtimes');
xlabel('Case');
ylabel('Runtime');
legend('Gaussian elimination','Cholesky algorithm');

