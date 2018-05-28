function [t,song] = singSyllable(rho2)
Fs = 22050;  %sampling
t = 0:1/Fs:0.5; 
 
 
[t1 y] = ode15s('birdBrain',t,[0.01, 0.01, 0.01],[],rho2);
 
kt=1.4*power(10,9)*y(:,3) + 4.8 * power(10,8);
pt=7000*y(:,1) - 2200;
 
[m z1] = ode15s('syrinx',t,[0.01 0.01],[],t,pt,kt);
song=z1(:,1);