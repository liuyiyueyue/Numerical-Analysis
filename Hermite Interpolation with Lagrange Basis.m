% HW06 Q1b Hermite Interpolation with Lagrange Basis---------------------------------

clc;
clear;

x=[0 2 3];
y=[1 0 2];
dy=[0 1 -1];
x_testing=linspace(0,3,100);
[p a b]=Hermite_Interpolation_Lagrange(x_testing,x,y,dy);



function [p a b]=Hermite_Interpolation_Lagrange(x_testing,x,y,dy)
% x_testing: discrete data points; 
% x: [x_1,...,x_n]
% y: [y_1,...,y_n]
% dy: derivatives of y
   n=length(x);          % number of interpolating points
   k=length(x_testing);  % number of testing points
   li=ones(n,k);         % Lagrange basis polynomials
   a=zeros(n,k);         % basis polynomials H_2i-1(x)
   b=zeros(n,k);         % basis polynomials H_2i(x)
   p=zeros(1,k);         % Hermie interpolation polynomial H(x)
   for i=1:n       
       dl=0;             % derivative of Lagrange basis
       for j=1:n    
           if j~=i
               dl=dl+1/(x(i)-x(j));
               li(i,:)=li(i,:).*(x_testing-x(j))/(x(i)-x(j)); 
           end
       end
       l2=li(i,:).^2;
       b(i,:)=(x_testing-x(i)).*l2;
       a(i,:)=(1-2*(x_testing-x(i))*dl).*l2;
       p=p+a(i,:)*y(i)+b(i,:)*dy(i);          % Hermite polynomial p(x)
   end 
   
    plot(x,y,'o');
    hold on;
    plot(x_testing,p)
    title('Q1b Hermite Interpolation with Lagrange Basis'); xlabel('x'); ylabel('y');
    legend("interploation points", "Hermite interpolation polynomial with Lagrange Basis");
end
