## Copyright (C) 2023 Leonardo Araujo <leolca@gmail.com>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING. If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {@var{L} = } isminphase. (@var{b}, @var{a})
## @deftypefnx {Function File} {@var{L} = } isminphase (@var{sos})
##
## Determine whether a digital filter is minimum phase. The filter might be defined 
## by the numerator coefficients, @var{b}, and the denominator coefficients, 
## @var{a}, or, alternatively, by a matrix of second-order sections, @var{sos}. 
##
## Example:
## @example
## a = [1 0.5]; b = [3 1];
## isminphase (b, a)
## @end example
##
## Ref [1] Oppenheim, Alan, and Ronald Schafer. Discrete-Time Signal Processing. 
## 3rd edition, Pearson, 2009.
## @end deftypefn

function flag = isminphase (b, a)
    if (nargin < 1 || nargin > 2 )
	print_usage;
    endif

    if (nargin == 2 && ( ! isrow(a) || ! isrow(b) ))
        error ( "coefficient array should be a row vector" );
    endif

    if (nargin == 1 && isrow (b))
        a = 1;
    endif

    if (nargin == 2 && isrow(a) && isrow (b)) || (nargin == 1 && isrow (b))
        zm = abs (roots (b));
	pm = abs (roots (a));
	flag = (all (zm < 1) || isempty (zm)) && (all (pm < 1) || isempty (pm));
    elseif (nargin == 1 && all(size(b)>[1 1]))
        [b, a] = sos2tf(b);
        flag = isminphase(b, a);
    endif

endfunction

%!demo
%! f = isminphase(b, a)

%! ## test input validation
%!error n = isminphase ()
%!error n = isminphase (1, 1, 1)
%!error n = isminphase (1, 1, 1, 1)
%!error n = isminphase (1, 1, 1, 1, 1)
%!error n = isminphase ([1:10]', 1)
%!error n = isminphase (1, [1:10]')
%!error n = isminphase ([1:10]', [1:10]')
%!error n = isminphase (1:10, 1:10, 1:10)
%!error n = isminphase (ones(3), ones(3))

%!test
%! b = [3 1];
%! a = [1 .5];
%! f = isminphase (b, a);
%! assert (f, true)

%!test
%! [b, a] = butter (1, .5); 
%! f = isminphase (b, a);
%! assert (f, false)

%!test
%! [b, a] = butter (8, .5);
%! f = isminphase (b, a);
%! assert (f, false)

%!test
%! b = 1.25^2 * conv( conv( conv([1 -0.9*e^(-j*0.6*pi)], [1 -0.9*e^(j*0.6*pi)]), [1 -0.8*e^(-j*0.8*pi)] ), [1 -0.8*e^(j*0.8*pi)] );
%! a = 1;
%! f = isminphase(b, a);
%! assert (f, true)
