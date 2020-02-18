function [miu,x] = PowerMethod(A,x0,Tol,N)

k = 1;
xp = norm(x0,Inf);
x = x0/xp;

while k < N
    y = A*x;
    yp = norm(y,Inf);
    miu = yp;
    if yp ==0
        disp("Eigenvector",x);
        disp("A has the eigenvalue 0, selecta new vector x and restart");
        break
    end
    err = x-(y/yp);
    error = norm(err,Inf);
    x = y/yp;
    if error < Tol
        printMiu = sprintf('Miu:',miu);
        printX = sprintf('Eigenvector',x);
        %printK = sprintf('Iteration',k);
        disp (printMiu);
        disp (printX);
        %disp (printK);
        break
    end
    k = k+1;
    if k == N
        disp("The maximum number of iterations exceeded");
        break
    end
end







