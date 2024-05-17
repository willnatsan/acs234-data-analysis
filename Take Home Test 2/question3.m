k = 2; m = 1; L = 2;
par = [k, m, L];

y0 = [20*pi/180, 0];
h = 0.1;
T = 10;

[tt, yE] = Euler4odesys(@(t,y) pendulum_model(t,y,par),[0, T], y0, h);
[tt, yH] = Heun4odesys(@(t,y) pendulum_model(t,y,par),[0, T], y0, h);

upper = 1*ones(size(tt));
lower = -1*ones(size(tt));

figure;
plot(tt,180/pi*yE(:,1), tt, 180/pi*yH(:,1), tt, upper, tt, lower); grid;
xlabel('Time(s)');
ylabel('\theta(t) (degree)');

legend('Euler', 'Heun', 'upper = 1', 'lower = -1')