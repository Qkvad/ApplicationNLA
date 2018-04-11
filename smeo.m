function smeo
m=[2;5;3;6]; %masses
M=inv(diag(sqrt(m)));
K=[24 -9 -5 0; -9 22 -8 -5; -5 -8 25 -7; 0 -5 -7 18];
A=M*K*M;
I=eye(4);

u0=ones(4,1);
u=u0/norm(u0);
ro=u'*A*u;

tol=1e-8; 
maxit=4;
k=1;
res(k)=norm(A*u-ro*u);
while res(k)>= tol
    y=minres(A-4*I, u, tol, maxit);
    u=y/norm(y);
    ro=u'*A*u;
    
    figure(1)
    plot(u,'-*');
    axis([0.5 4.5 -1 1]);
    xlabel('Indeks komponente');
    ylabel('Komponenta');
    title(sprintf('Aproksimacija svojstvenog vektora: iteracija %d',k));
    grid on
    pause(0.5);
    k=k+1;
    res(k)=norm(A*u-ro*u);
end

figure(2)
semilogy(1:k,res(1:k))
xlabel('iteracija');
ylabel('norma reziduala');
ro=u'*A*u

end