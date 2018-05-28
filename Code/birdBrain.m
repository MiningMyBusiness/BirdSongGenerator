% 
%The CPG (Central Pattern Generator) Neural Network
%Modeling the RA brain region of a songbird
%From Laje and Mindlin(2002), Diversity within a Birdsong, Physical Review
%Letters
%the parameter rho2 is varied where all other parameters are fixed in
%generating birdsong syllables
function yprime = birdBrain(t,y,flag,rho2)
rho1=0;
rho3=6;
A=10;
B=10;
C=10;
D=-2;
E=4;
alpha=2;
beta=20;
%xp
yprime(1) = 30 * ( - y(1)  + ( 1 / ( 1 + exp(-1 * (rho1  + A * y(1) - B *y(2)))) )) ;
%y
yprime(2) = 30 * ( - y(2)  + ( 1 / ( 1 + exp (-1 * ( rho2 + C * y(1) - D*y(2) + alpha * y(3))))));
%xk
yprime(3) = 120 * ( -y(3)  + ( 1 / ( 1 + exp ( -1 * (rho3 + E* y(3) - beta*y(2))))));
yprime = yprime';