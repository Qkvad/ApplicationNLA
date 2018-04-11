function deconvolution
  load('model_dekonvolucija_uy.mat') %input
  load('model_dekonvolucija_hh.mat') %exact solution
    
  m=102; n=18;
  
  for j = 1:m+1
    U(j,:) = up(:,n+j:-1:j);   % matrica U dimenzije 103x19
  end

  Y = yp;     % matrica Y dimenzije 103x1;     sustav je Y = UH
  C = [U Y];
  s1 = svd(U);
  s2 = svd(C);
  if s1(19) > s2(20)
    [P,S,R]=svd(C,0);
    P2 = P(:,20);
    S2 = S(20,20);
    R_12 = R(1:19,20);
    R_22 = R(20,20);
    Z = [R_12' R_22'];
    RJ = -P2*S2*Z;
    E0 = RJ(:,1:19);
    R0 = RJ(:,20);
    X = (U+E0)\(Y+R0);
    
  end

  plot(hh,'r');
  hold on;
  plot(X,'g');
  legend('egzaktnoRJ','izracunatoRJ');
  hold off;

end
