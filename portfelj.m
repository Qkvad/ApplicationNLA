## Copyright (C) 2018 Petra
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} portfelj (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Petra <petra@petra>
## Created: 2018-03-13

function [rj1, rj2, rj1sor, rj2sor] = portfelj (C,m)
mip=0.05;
[msize,n]=size(C);
e=ones(n,1);
x0=zeros(n,1);
tol=1e-8;

%lijeva i desna strana prvog sustava (za omega min)
invC=inv(C);
b1=invC*e;
M1=e' * b1;

%pomocne vrijednosti za drugi sustav za omega mip
a1=invC*m;
a2=m.' * a1;
a3=e.' * a1;

%lijeva i desna strana drugog sustava (za omega mip)
M2=M1*a2-(a3).^2;
b2=(a2-mip*a3)*b1 + (mip*M1-a3)*a1;

%sor konvergencija i sor rjesavac 1.sustav
omega1=sor_konvergencija(M1);
[rj1sor,iter1sor,resvec1sor]=sor(M1,b1,x0,tol,omega1);

%sor konvergencija i sor rjesavac 1.sustav
x0=zeros(n,1);
omega2=sor_konvergencija(M2);
[rj2sor,iter2sor,resvec2sor]=sor(M2,b2,x0,tol,omega2);

%pcg sustav 1
[rj1,flag1,relres1,iter1,resvec1,eigest1]=pcg(M1,b1,tol,n,[],[],[]);

%pcg sustav 2
[rj2,flag2,relres2,iter2,resvec2,eigest2]=pcg(M2,b2,tol,n,[],[],[]);


interval=1:max(iter1sor,iter2sor);

figure
semilogy(interval, resvec1sor,'r')
hold on
semilogy(interval, resvec2sor,'b')


endfunction
