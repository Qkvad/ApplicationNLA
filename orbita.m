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
## @deftypefn {} {@var{retval} =} orbita (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Petra <petra@petra>
## Created: 2018-03-17

function [A,b,a] = orbita
load model_orbita_polozaji.mat
P=[p1;p2;p3;p4;p5];
k=1;
for i=1:5
  A(i,1)=P(k)*P(k+1);
  A(i,2)=P(k+1)*P(k+1);
  A(i,3)=P(k);
  A(i,4)=P(k+1);
  A(i,5)=1;
  b(i)=-(P(k)*P(k));
  k=k+2;
end
b=b';

tol=1e-4;
x0=zeros(5,1);
[a]=gmres(A,b,[],tol,5,[],[],x0);

xstep=1; ystep=2;
for i=1:5
  x(i)=P(xstep);
  xstep=xstep+2;
  y(i)=P(ystep);
  ystep=ystep+2;
end
f=fja_orbita(x,y,a);

x1=-5:0.1:5;
y1=x1;
[X,Y]=meshgrid(x1,y1);
Z=fja_orbita(X,Y,a);
contour(x1,y1,Z,[0 0]);
grid on

surf(x1,y1,Z)
endfunction
