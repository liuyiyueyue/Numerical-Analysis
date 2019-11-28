% HW06 Q1c Hermite Interpolation with Newtone Basis------------------------------------------------------

clc;
clear;

n=2; % total of n+1 interpolation points
% (x,y,y')= (0,1,0), (2,0,1), (3,2,-1)
 x=[0 2 3];
 y=[1 0 2]; 
 q=[1 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 1 0 0 0 0;
    2 0 0 0 0 0;
    0 -1 0 0 0 0];
 %{
    q = [f(x1) 0 0 0 0 0;
         0 f'(x1) 0 0 0 0;
         f(x2) 0 0 0 0 0;
         0 f'(x2) 0 0 0 0;
         f(x3) 0 0 0 0 0;
         0 f'(x3) 0 0 0 0]
    q is the starting matrix of Divided Differences table      
 %}
    
 z = zeros(1,2*n+2);
 for i = 0:n
    z(2*i+1) = x(i+1);
    z(2*i+2) = x(i+1);
    q(2*i+2,1) = q(2*i+1,1);
    if i ~= 0 
       q(2*i+1,2) = (q(2*i+1,1)-q(2*i,1))/(z(2*i+1)-z(2*i));
    end
 end

 k = 2*n+1;
 for i = 2:k
    for j = 2:i
       q(i+1,j+1) = (q(i+1,j)-q(i,j))/(z(i+1)-z(i-j+1));
    end
 end

 fprintf('Coefficients of the Hermite polynomial are:\n ');
 for i = 0:k
    fprintf(' %11.8f \n', q(i+1,i+1));
 end
 
 x_testing = linspace(0,3,100)';
 for i=1:100
     s(i) = interpolationEval(x_testing(i), k, q, z);
 end
 
 plot(x,y,'o');
 hold on;
 plot(x_testing, s);
 title('Q1c Hermite Interpolation with Newtone Basis'); xlabel('x'); ylabel('y');
 legend("interploation points", "Hermite interpolation polynomial with Newtone Basis");
  
 
 function s = interpolationEval(xx, k, q, z)
     s = q(k+1,k+1)*(xx-z(k));
     for i = 2:k
        j = k-i+1;
        s = (s+q(j+1,j+1))*(xx-z(j));
     end
     s = s + q(1,1);
 end


