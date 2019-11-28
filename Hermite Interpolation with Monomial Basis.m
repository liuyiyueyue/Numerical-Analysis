% HW06 Q1a Hermite Interpolation with monomial Basis------------------------------------------------------

clc;
clear;
set(0, 'defaultaxesfontsize',15,'defaultaxeslinewidth',1.5,...
       'defaultlinelinewidth',1.5,'defaultpatchlinewidth',1.5,...
       'defaulttextfontsize',15)

% calculate the coeff a_1...a_2n
V = [1 0 0 0 0 0;
     1 2 4 8 16 32;
     1 3 9 27 81 243;
     0 1 0 0 0 0;
     0 1 4 12 32 80;
     0 1 6 27 108 405];
 
y = [1; 0; 2; 0; 1; -1];
a = V\y
 
 % evalute the interpolation polynomial p_2n-1(x) on [0,3]
 x0 = [0 2 3]';
 y0 = [1 0 2]';
 x = linspace(0,3,100)';
 pVal = HornerV(a,x);
 plot(x0,y0,'o');
 hold on;
 plot(x,pVal);
 axis([-1 4 -1 3])
 title('Q1a Hermite Interpolation with monomial Basis'); xlabel('x'); ylabel('y');
 legend("interploation points", "Hermite interpolation polynomial with monomial Basis");
 
 
 function pVal = HornerV(a,z)
% pVal = HornerV(a,z)
% evaluates the Vandermonde interpolant on z where
%   a is an n-vector
%   z is an m-vector.
%   pVal is a vector the same size as z, where p(x) = a(1) + .. +a(n)x^(n-1)

     n = length(a);
     m = length(z);
     pVal = a(n) * ones(size(z));
     for k = n-1:-1:1
         pVal = z .* pVal + a(k);
     end
 end
 
 
