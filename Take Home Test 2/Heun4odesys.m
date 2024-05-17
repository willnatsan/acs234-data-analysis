function [tout, yout] = Heun4odesys(FunFcn, tspan, y0, ssize)

% Heun	
%       Integrates a system of ordinary differential equations using
%	the Heun's method.
%	[t,y] = Heun('yprime', tspan, y0) integrates the system
%	of ordinary differential equations described by the M-file
%	yprime.m over the interval tspan = [t0,tfinal] and using initial
%	conditions y0.
%	[t, y] = Heun(F, tspan, y0, ssize) uses step size ssize
%
% INPUT:
% F     - String containing name of user-supplied problem description.
%         Call: yprime = fun(t,y) where F = 'fun'.
%         t      - Time (scalar).
%         y      - Solution vector.
%         yprime - Returned derivative vector; yprime(i) = dy(i)/dt.
% tspan = [t0, tfinal], where t0 is the initial value of t, and tfinal is
%         the final value of t.
% y0    - Initial value vector.
% ssize - The step size to be used. (Default: ssize = (tfinal - t0)/100).
%
% OUTPUT:
% t  - Returned integration time points (column-vector).
% y  - Returned solution, one solution row-vector per tout-value.
%
% The result can be displayed by: plot(t,y).
%
%
%
% Initialization
% Written by Prof. Deborah Sulsky
%
t0=tspan(1);
tfinal=tspan(2);
pm = sign(tfinal - t0);  % Which way are we computing?
if (nargin < 4), ssize = abs(tfinal - t0)/100; end
if ssize < 0, ssize = -ssize; end
h = pm*ssize;  
t = t0;
y = y0(:);
tout = t;
yout = y.';

% To compute the number of steps.

dt = abs(tfinal - t0);
N = floor(dt/ssize) + 1;
if (N-1)*ssize < dt
  N = N + 1;
end

% Initialize the output.

tout = zeros(N,1);
tout(1) = t;
yout = zeros(N,size(y,1));
yout(1,:) = y.';
k = 1;

% The main loop

while k < N
  if pm*(t + h - tfinal) > 0
    h = tfinal - t;
    tout(k+1) = tfinal;
  else
    tout(k+1) = t0 +k*h;
  end
  k = k+1;
  % Compute the slope
  s1 = feval(FunFcn, t, y); s1 = s1(:); 
  ytemp=y+h*s1;  
  s2= feval(FunFcn,t+h,ytemp); s2 = s2(:);
  y = y + h*(s1+s2)/2;		
  t = tout(k);
  yout(k,:) = y.';
end;