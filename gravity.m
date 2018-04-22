function gravity
  load('model_gravitacija_g.mat')
  m=15; hs=1/(m-1); s=0:hs:1;
  n=15; ht=1/(n-1); t=0:ht:1;
  d=0.25;
  %ksi=normrnd(0,0.1,[1,m]);
  ksi=0.1*ones(m,1);
  g=g-ksi;
  %exact solution
  for j=1:n
    f(j)=sin(pi*t(j)) + 0.5*sin(2*pi*t(j));
  end
    
  for i=1:m
    K(i,1)=(d/((d^2+(s(i)-t(j))^2)^(3/2)))/(2*(n-1));
    K(i,n)=(d/((d^2+(s(i)-t(j))^2)^(3/2)))/(2*(n-1));
  end
  
  for i=1:m
    for j=2:(n-1)
      K(i,j)=(d/(d^2+(s(i)-t(j))^2)^(3/2))/(n-1);
    end
  end
      
  [U,S,V]=svd(K);
  sigma=diag(S);
  
  figure(1)
  plot(sigma,'o')
  ylabel('Singular values of K')
    
  figure(2)
  % solution x=V(:,1:r)*inv(S(1:r,1:r))*U(:,1:r)'*b
  for p=1:14
    %x=V(:,1:p)*(S(1:p,1:p)\(U(:,1:p)'*g)); 
    x=V(:,1:p)*inv(S(1:p,1:p))*U(:,1:p)'*g;
    plot(t,x,'r',t,f)
    legend('approx','exact')
    title(sprintf('rank %d',p+1))
    pause(0.2)
  end
  
  
end
