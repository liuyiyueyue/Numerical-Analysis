% Q2b Piecewise Cubic Hermite Interpolating Polynomial (PCHIP) on Runge's ------------------------

clc;
clear;

func = @(x) 1/(1+x^2);
funcDeriv = @(x) -2*x/((x^2+1)^2);

n = 10; % or n = 10
pointx = linspace(-5,5,n);
pointy = arrayfun(func,pointx);
dy = arrayfun(funcDeriv,pointx);    % derivative of the function at pointx
x_testing = -5:0.01:5;
y = Piecewise_Cubic_Hermite(x_testing,pointx,pointy,dy);


function p = Piecewise_Cubic_Hermite(x,pointx,pointy,dy)
% Piecewise Hermite cubic interpolation between 2 points
%
% INPUT: 
% x = an arbitrary vector that will be interpolated
% pointx = data points of the independent variable 
%          (The points do not have to be equally spaced)
% pointy = data points of the dependent variable. pointy is the value of
%          the function at pointx
% dy = data points of the dependent variable's derivative. yprime is the
%          derivative of the function at pointx
%
% OUTPUT:
% p = the piecewise interpolation polynomial of a vector "x". 

	n = length(pointx);                 % number of interpolating points
	p = zeros(1,length(x));
	dx = (x(end)-x(1))/(length(x)-1);   % number of pieces(segments)
	start = 1; 
    
	for s = 1:n-1                       % loops through all segments and interpolates in the segments
        interval = round((pointx(s+1)-pointx(s))/dx);    
        xi = x(start:interval+start);   % testing points between two interpolation points 
        H0 = -(((xi-pointx(s+1)).^2.*(2*xi-3*pointx(s)+pointx(s+1)))./(pointx(s)-pointx(s+1))^3); % H0
        H1 = ((xi-pointx(s)).^2.*(2*xi+pointx(s)-3*pointx(s+1)))./(pointx(s)-pointx(s+1))^3;      % H1
        h0 = ((xi-pointx(s)).*(xi-pointx(s+1)).^2)./(pointx(s)-pointx(s+1))^2;                    % in fact, it is h0*(x_i+1 - x_i)
        h1 = ((xi-pointx(s)).^2.*(xi-pointx(s+1)))./(pointx(s)-pointx(s+1))^2;                    % in fact, it is h1*(x_i+1 - x_i)
        p(start:interval+start) = (pointy(s)*H0)+(pointy(s+1)*H1)+(dy(s)*h0)+(dy(s+1)*h1);
        start = start+interval;
    end

	plot(pointx,pointy,'o'); 
	hold on;
	plot(x,p,'-k'); 
    hold on;
	title('Piecewise Cubic Hermite Interpolation on Runge function'); xlabel('x'); ylabel('y'); 
    legend("interploation points", "PCHIP");
    
end
