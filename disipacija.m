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
## @deftypefn {} {@var{retval} =} disipacija (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Petra <petra@petra>
## Created: 2018-03-14

function [x,iter,resvec]=disipacija 
  load var_disipacija
  [m,n]=size(A);
  x0=ones(n,1); 
  tol=1e-8;
  [x,flag,relres,iter,resvec]=gmres(@(x0)mdAx(x0),b,[],tol,n,[],[],x0);
  %interval=1:iter
  plot(resvec)

endfunction
