function [omegamin]=sor_konvergencija(A)
D=diag(diag(A));
U=triu(A-D);
L=tril(A-D);
ind=1;
omega(ind)=0;
T=(inv(D+omega(ind)*L))*((1-omega(ind))*D-omega(ind)*U);
ro(ind)=max(abs(eig(T)));
while omega(ind)<=2
   ind=ind+1;
   omega(ind)=omega(ind-1)+0.01;
   T=(inv(D+omega(ind)*L))*((1-omega(ind))*D-omega(ind)*U);
   ro(ind)=max(abs(eig(T)));
end
plot(omega,ro)
xlabel('omega')
ylabel('ro')

[romin,ind]=min(ro);
omegamin=omega(ind);
end
