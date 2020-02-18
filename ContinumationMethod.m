function x = ContinumationMethod(A,N,x0)
syms a b
fx = subs(A,[a,b],[x0(1,1),x0(2,1)]);
i = 0;
h = 1/N;
B = -h*fx;
x = x0;
while i < N
    Jx = jacobian(A);
    k1= inv(Jx)*B;
    x = x+k1;
    x = subs(x,[a,b],[x0(1,1),x0(2,1)]);
    x0 = x;
    i = i+1;
end
end
