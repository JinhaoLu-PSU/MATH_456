function [x,nit] = NewtonMethod(A,x0,Tol,MaxIter)
syms a b 
nit = 1;
%t = zeros(2,2);
while nit < MaxIter
    Jx = jacobian(A);
    Fx = transpose(A);
    s = inv(Jx)*(-Fx);
    s = subs(s,[a,b],[x0(1,1),x0(2,1)]);
    x = x0+s;
    x0 = double(x);
    if norm(s)<Tol
        break;
    end
    nit = nit+1;
end
end
