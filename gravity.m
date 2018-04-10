function gravity
  load('model_gravitacija_g.mat')
  m=15; n=15; d=0.25; %discretization s=[0,1] with m and n dots, depth d
  stddeverr=0.1; %standard deviation error
  ksi=normrnd(0,0.1,[1,m]); 
  for i=1:m
    s(i)=(i-1)/(m-1);
  end
  for j=1:n
    t(j)=(j-1)/(n-1);    
  end
  for i=1:m
    for j=1:n
      K(i,j)=d/(( d^2 + (s(i)-t(j))^2 )^(3/2));  
    end
  end
  [U,S,V]=svd(K);
  sigma=diag(S);
  
  figure(1)
  int=1:15;
  plot(int,sigma,'-o')
  ylabel('Singular values of K')
  axis([1 15])
  
  %exact solution
  for j=1:n
    f(j)=sin(pi*t(j)) + 0.5*sin(2*pi*t(j));
  end
  
  b=g-ksi;
  
  figure(2)
  % solution x=V(:,1:r)*inv(S(1:r,1:r))*U(:,1:r)'*b
  for p=1:15
    x=V(:,1:p)*inv(S(1:p,1:p))*U(:,1:p)'*b;
    plot(int,x,int,f)
    legend('approx','exact')
    axis([1 15 -6 6])
    title(sprintf('rank %d',p))
    pause(0.2)
  end
  
  
end