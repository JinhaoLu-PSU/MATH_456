function [x,nit] = ModifyNewton(A,x0,Tol,MaxIter)
syms a b 
nit = 1;
%t = zeros(2,2);
Jx = jacobian(A);
JLU = subs(Jx,[a,b],[x0(1,1),x0(2,1)]);
while nit < MaxIter
    Fx = transpose(A);
    [L,U] = lu(JLU);
    y = inv(L)*(-Fx);
    s = inv(U)* y ; 
    s = subs(s,[a,b],[x0(1,1),x0(2,1)]);
    x = x0+s;
    x0 = x;
    if norm(s)<Tol
        break;
    end
    nit = nit+1;
end
end
