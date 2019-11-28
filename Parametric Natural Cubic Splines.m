% HW06 Q3 parametric natural cubic splines ------------------------------------

clc;
clear;

n = 12; % total of n points
x = [25.0 19.0 13 9.0 5.0 2.2 1.0 3.0 8.0 13.0 18.0 25.0];
y = [5.0 7.5 9.1 9.4 9.0 7.5 5.0 2.1 2.0 3.5 4.5 5.0];

t(1) = 0;
for i = 2:n
    t(i) = t(i-1) + ( (x(i)-x(i-1))^2 + (y(i)-y(i-1))^2 ) ^ 0.5;
end
z = linspace(-t(1), t(12), 200); %vector containing the points at which the spline has to evaluated

% x(t) = a1(i) + b1(i)*(t-t(i)) + c1(i)*(t-t(i))^2 + d1(i)*(t-t(i))^3
[a1, b1, c1, d1, iflag] = cfspline( t, x );
s1 = csplineeval( t, a1, b1, c1, d1, z ); % return the value of x

% y(t) = a2(i) + b2(i)*(t-t(i)) + c2(i)*(t-t(i))^2 + d2(i)*(t-t(i))^3
[a2, b2, c2, d2, iflag] = cfspline( t, y );
s2 = csplineeval( t, a2, b2, c2, d2, z ); % return the value of y

plot(s1, s2, 'k'); hold on;
plot(x, y, 'o'); hold on;
title('Q3 parametric natural cubic splines'); xlabel('x'); ylabel('y'); hold on; 
legend('parametric natural cubic splines', 'interpolation points')
