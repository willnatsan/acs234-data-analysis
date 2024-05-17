function dydt = pendulum_model(t, y, par)
   k = par(1);
   m = par(2);
   L = par(3);
   g = 9.81;
   theta = y(1);
   theta_dot = y(2);
   theta_ddot = -(g/L) * sin(theta) -(k/m) * theta_dot;
   dydt = [theta_dot; theta_ddot];
end