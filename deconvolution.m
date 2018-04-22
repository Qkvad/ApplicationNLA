function deconvolution
  load('model_dekonvolucija_uy.mat') %input
  load('model_dekonvolucija_hh.mat') %exact solution
    
  m=103; n=19;
  
  for j=1:m
    U(j,:)=up(:,n+j-1:-1:j);   % UR{103x19}
  end

  Y = yp;      
  C = [U Y];
  s1 = svd(U); [s1m,s1n]=size(s1);
  s2 = svd(C); [s2m,s2n]=size(s2);
  if s1(s1m) > s2(s2m)
    [P,S,R]=svd(C,0);
    P2 = P(:,s2m);
    S2 = S(s2m,s2m);
    R12 = R(1:s1m,s2m);
    R22 = R(s2m,s2m);
    Q = [R12' R22'];
    E0R0 = -P2*S2*Q;
    E0 = E0R0(:,1:s1m);
    R0 = E0R0(:,s2m);
    X = (U+E0)\(Y+R0);
  end

  plot(1:s1m,hh,'b',1:s1m,X,'r');
  legend('exact solution','computed solution');  
  axis([1 19 -1 1])
end
