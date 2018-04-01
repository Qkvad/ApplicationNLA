function [x, iter, vecres] = sor(A, b, x0, tol, omega)
D=diag(diag(A));
U=triu(A-D);
L=tril(A-D);
iter=1;
normb=norm(b,2);
vecres(iter)=(norm(b - A*x0,2))/normb;
while vecres(iter)>tol
    x=(inv(D+omega*L))*(((1-omega)*D-omega*U)*x0+omega*b);
    iter=iter+1;
    vecres(iter)=(norm(b - A*x,2))/normb;
    x0=x;
end
end

