% HW06 Q4 least square problems ----------------------------------------------------
clc;
format long;


func = @(x) cos(4*x);
pointx = linspace(0,1,50);
b = arrayfun(func, pointx)';
A = zeros(50,12);
[m,n] = size(A);

for i=1:50
    for j=1:12
        A(i,j) = pointx(i)^(j-1);
    end
end


% 1. --------------------------------------------------------------------------
% By solving the normal equations 
% Rearranging, we have x = A?b := (A?A)^(?1) A? b,
Moore_Penrose_pseudoinverse = (transpose(A)*A)^(-1) * transpose(A);
x1 = Moore_Penrose_pseudoinverse * b


% 2. ---------------------------------------------------------------------------
% By using the QR factorization built by the code mgs_qr.m on Canvas
% (modified Gram-Schmidt). A = Q * R
[Q2,R2] = mgs_qr(A);
c2 = Q2'*b;
x2 = inv(R2(1:n,1:n)) * c2(1:n)

% 3. -------------------------------------------------------------------------------
% By using the QR factorization built by the code slow_householder_qr.m on Canvas (Householder reflectors).
[Q3,R3] = slow_householder_qr(A);
c3 = Q3'*b;
x3 = inv(R3(1:n,1:n)) * c3(1:n)


% 4. --------------------------------------------------------------------------------
% By using the QR factorization built by MATLAB's qr.
[Q4,R4] = qr(A);
c4 = Q4'*b;
x4 = inv(R4(1:n,1:n)) * c4(1:n)

% 5. --------------------------------------------------------------------------------
% By using backslash (just compute x = A\b).
x5 = A\b

% 6. --------------------------------------------------------------------------------
% By using the SVD built by MATLAB's svd. 
[U,S,V] = svd(A);
s_diag = diag(S);
% determine the effective rank r of A using singular values
r = 1;
while( r < size(A,2) & s_diag(r+1) >= max(size(A))*eps*s_diag(1) )
    r = r+1; 
end
temp = U'*b;
x6 = V * ( [temp(1:r)./s_diag(1:r); zeros(n-r,1) ] ) % max norm solution
