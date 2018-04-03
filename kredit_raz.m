function kredit_raz (n)
%n is the max number of iterations

load('model_kredit_A.mat')    %matrix
u0=zeros(8,1); u0(1)=1;       %u0=[1 0 0 0 0 0 0 0]^T

u=A*u0;
for k=0:n
  p=plot(1:8,u);
  title(sprintf("Iteracija broj %d",k));
  axis([0.5 8.5 0.0 1.0])
  grid on
  xlabel('Stanje')
  ylabel('Gustoća')
  set(gca,'XTick',1:8)
  set(gca,'XTickLabel', { 'AAA', 'AA', 'A','BBB', 'BB', 'B', 'CCC', 'D' } );
  
  pause(0.2);
  
  u=A*u;
end

[vec,vrij]=eig(A)

end
