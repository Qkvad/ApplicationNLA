function phillips_eq
  m=150; n=121;
  for i=1:m
    s(i)=-5.925 + (i-1) * 11.85/(m-1);     
  end
  for j=1:n
    t(j)=-3 + (j-1) * 6/(n-1);    
  end
  for i=1:m
    y(i)=1/6 * ( (6-abs(s(i)))*(1+ 1/2 * cos((pi*s(i))/3)) + 9/(2*pi) * sin((pi*abs(s(i)))/3) );
  end
  for j=1:n
    x(j)=1 + cos((pi*t(j))/3);    
  end
  
  for i=1:m
    if abs(t(1)-s(i))<=3
      K(i,1)= 3/(n-1) * ( 1/6 * ( 1 + cos( (pi*(t(j)-s(i)))/3 ) ) );
    else
      K(i,1)=0;
    end
    if abs(t(n)-s(i))<=3
      K(i,n)=3/(n-1) * ( 1/6 * ( 1 + cos( (pi*(t(j)-s(i)))/3 ) ) );
    else
      K(i,n)=0;
    end
  end
  
  for i=1:m
    for j=2:(n-1)
      K(i,j)=6/(n-1);
      if abs(t(j)-s(i))<=3
        K(i,j)=K(i,j)*( 1/6 * ( 1 + cos( (pi*(t(j)-s(i)))/3 ) ) );
      else
        K(i,j)=0;
      end
    end
  end
  
  for i=1:m
    zeta(i)=1e-4 * y(i);    
  end
  S=diag(zeta);
  
  eta=normrnd(0,1,[1,m]); eta=eta';
  y=y';
  A=inv(S)*K; 
  b=inv(S)*y+eta;
  
  rangA=rank(A);
  [mA,nA]=size(A);
  [mb,nb]=size(b);
  
  [U,S,V]=svd(A);
  xNK=V*inv(S(1:n,:))*U(:,1:n)' *b;
  
  x=x'; 
  normxNK=norm(A*xNK-b)
  normdiff=norm(xNK-x)
  
  figure(1)
  plot(-3:6/(n-1):3,x,-3:6/(n-1):3,xNK)
  axis([-3 3 -100 100])
  xlabel('t')
  ylabel('x_{NK}(t)')
  legend('x(t)','x_{NK}(t)')
  
  %-----------------------------------------------------------------------------
  %                         regularisation considered
  %-----------------------------------------------------------------------------
  lambdaopt=0.748;
  In=eye(n);
  for j=1:n
    yopt(j)=1/6 * ( (6-abs(t(j)))*(1+ 1/2 * cos((pi*t(j))/3)) + 9/(2*pi) * sin((pi*abs(t(j)))/3) );
  end
  Aopt=[A;lambdaopt*In]; 
  yopt=yopt';
  bopt=[b;lambdaopt*yopt];
    
  [Uopt,Sopt,Vopt]=svd(Aopt);
  [mUopt,nUopt]=size(Uopt)
  [mSopt,nSopt]=size(Sopt)
  [mVopt,nVopt]=size(Vopt)
  xopt=Vopt*inv(Sopt(1:n,:))*Uopt(:,1:n)'*bopt;
  xopt=xopt'; 
  normxopt=norm(A*xopt'-b)
  normdiffopt=norm(xopt' -x)
  figure(2)
  plot(-3:6/(n-1):3,x,-3:6/(n-1):3,xopt)
  axis([-3 3 -0.5 3])
  xlabel('t')
  ylabel('x_{opt}(t)')
  legend('x(t)','x_{opt}(t)')
  
  %-----------------------------------------------------------------------------
  %                         statistical lambda
  %-----------------------------------------------------------------------------
  lambdanaj=77.5;
  Anaj=[A;lambdanaj*In]; 
  bnaj=[b;lambdanaj*yopt];
    
  [Unaj,Snaj,Vnaj]=svd(Anaj);
  xnaj=Vnaj*inv(Snaj(1:n,:))*Unaj(:,1:n)'*bnaj;
  xnaj=xnaj'; 
  normxnaj=norm(A*xnaj'-b)
  normdiffnaj=norm(xnaj' -x)
  figure(3)
  plot(-3:6/(n-1):3,x,-3:6/(n-1):3,xnaj)
  axis([-3 3 -0.5 3])
  xlabel('t')
  ylabel('x_{naj}(t)')
  legend('x(t)','x_{naj}(t)')
end
