function [miu,x] = SymPowerMethod(A,x0,Tol,N)

k = 1; 
nX2 = norm(x0);
x = x0/nX2;

while k <= N
    y = A*x;
    miu = transpose(x)*y;
    nY2 = norm(y);
    if nY2 == 0
        disp ('Eigenvector',x);
        disp("A has the eigenvalue 0, selecta new vector x and restart");
        break
    end
    error = x - y/nY2;
    x = y/nY2;
    if error < Tol
        disp(miu);
        disp(x);
        break
    end
    k = k+1;
    if k == N
        disp("The maximum number of iterations exceeded");
        break
    end
end