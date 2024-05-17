function [ x, y ] = Euler4odesys( f, x_range, y_initial, step_size )

%
% OUTPUTS:
% x, the value of the independent variable at nstep+1 equally spaced points; 
% y, the value of the dependent variable at each x.
%
% INPUTS:
%
% f, the name of the derivative function, which will have to be in quotes; 
% x_range, a row vector of two values, the starting and stopping values of the independent variable; 
% y_initial, the value of the dependent variable at the initial time; 
% nstep, the number of Euler steps to take
%
x(1) = x_range(1);
%dx = ( x_range(2) - x_range(1) ) / nstep;
dx = step_size;
nstep = ( x_range(2) - x_range(1) )/dx;
y(:,1) = y_initial;

for i = 1 : nstep
x(i+1) = x(i) + dx;
y(:,i+1) = y(:,i) + dx * feval ( f, x(i), y(:,i) );
end 
x=x';
y = y';