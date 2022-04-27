clear all; close all; clc;
K0   = 1.165;
tau0 = 9.102*10^(-2);
T0   = 3.774*10^(-3);
num  = K0;
den  = [tau0 1];
g    = tf(num,den,'InputDelay',T0);
bode(g,'blue');
hold on
K1    = 1;
K2    = 33.11;
sys1  = feedback(K1*g,1,-1);
sys2  = feedback(K2*g,1,-1);
display(sys1);
bode(sys1,'magenta');
hold on
bode(sys2,'green');
box on;
legend('model teòric','feedback amb K=1','feedback amb Kp')
grid on
ax = gca; 
ax.GridColor = [0,0,0];
ax.GridAlpha = 0.3;
ax.MinorGridColor = [0,0,0];
ax.MinorGridAlpha = 0.3;
saveas(gcf,'Bode.fig')
print('Bode1.pdf','-dpdf','-fillpage')