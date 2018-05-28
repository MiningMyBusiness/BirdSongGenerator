%The bird syrinx model used in
%From Laje and Mindlin(2002), Diversity within a Birdsong, Physical Review
%Letters
%pt = interlabial air pressure
%k:labia elasticity
function yprime = syrinx(t,y,flag,tra,p,k)
 
b=1000;  %(dissipation-friction coefficient) parameter taken from the paper
d=power(10,8);  %(nonlinear dissipation coefficent for bounded motion) parameter taken from the paper
 
yprime(1) = y(2) ;
 
pt=interp1(tra,p,t,'spline');
kt=interp1(tra,k,t,'spline');

yprime(2) = (pt-b)*y(2) - kt*y(1) - d*(power(y(1),2)*y(2));
yprime = yprime';