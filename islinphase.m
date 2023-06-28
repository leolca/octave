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
## @deftypefn  {Function File} {@var{flag} = } islinphase (@var{b}, @var{a})
## @deftypefnx {Function File} {@var{flag} = } islinphase (@var{sos})
## @deftypefnx {Function File} {@var{flag} = } islinphase (@dots{}, @var{tol})
##
## Determine whether a digital filter has linear phase. Returns a flag
## equal to true if the given filter has linear phase and false otherwise.
## The filter might be defined by the numerator coefficients, @var{b}, and the 
## denominator coefficients, @var{a}, or, alternatively, by a matrix of 
## second-order sections, @var{sos}. A tolerance @var{tol} might be given to
## define when two numbers are close enough to be considered equal.
##
## Example:
##
## @example
## b = [1 2 4 4 2 1];
## freqz (b);
## islinphase (b)
## @end example
##
## References:
##
## [1] Oppenheim, Alan, and Ronald Schafer. Discrete-Time Signal Processing. 
## 3rd edition, Pearson, 2009.
## 
## [2] Paquelet, Stéphane, and Vincent Savaux. “On the Symmetry of FIR Filter 
## with Linear Phase.” Digital Signal Processing, vol. 81, Oct. 2018, pp. 57–60. 
## https://doi.org/10.1016/j.dsp.2018.07.011.
## @end deftypefn

function flag = islinphase (b, a, tol)
    if (nargin < 1 || nargin > 3 )
	print_usage;
    endif

    if (nargin == 2 && ( ! isrow (a) || ! isrow (b) ))
        error ( "coefficient array should be a row vector" );
    endif

    if (nargin == 1 && isrow (b)), a = 1; endif

    if nargin < 3, 
        tol = eps^(3/4); 
    elseif length (tol) > 1,
	error ( "a scalar is expected as the tolerance value" );
    endif

    if (nargin > 1 && isrow (a) && isrow (b)) || (nargin == 1 && isrow (b))
        % remove leading and trailing zeros
	b = b(find (b, 1):end);
	a = a(find (a, 1):end);
	b = b(1:find (b,1,"last"));
	a = a(1:find (a,1,"last"));
        % True linear phase characteristics are not achievable with traditional 
	% IIR (Infinite Impulse Response) filters due to their recursive nature.
	if length (a) == 1 % FIR filter
	    % The corresponding impulse response of a linear phase system may
	    % have even symmetry about α, if 2α is an integer; i.e., 
	    % h[2α − n] = h[n].
	    % If it is a FIR filter, it suffices to check if the numerator is symmetric.
	    if mod ( length(b), 2) == 0,
	       if all (b(1:length(b)/2) == fliplr (b(length(b)/2+1:end)))
	           flag = true;
	       else
	           flag = false;
	       endif
	    else
	       if all (b(1:floor (length(b)/2)) == fliplr (b(ceil (length(b)/2)+1:end)))
	           flag = true;
	       else
		   flag = false;
	       endif
	    endif
	
	else
	    % Generalized linear phase systems have constant group delay. 
	    % Use grpdelay and check if it is constant up to a given tolerance.
	    gd = grpdelay(b, a, 128);
	    mgd = mean (gd);
	    if all (abs (gd - mgd) < tol)
	        flag = true;
	    else
		flag = false;
	    endif
	endif
    elseif (nargin == 1 && all(size (b) > [1 1]))
        [b, a] = sos2tf (b);
        flag = islinphase (b, a, tol);
    endif

endfunction

%!demo
%! f = islinphase (b, a)

%! ## test input validation
%!error n = islinphase ()
%!error n = islinphase (1, 1, 1, 1)
%!error n = islinphase (1, 1, 1, 1, 1)
%!error n = islinphase ([1:10]', 1)
%!error n = islinphase (1, [1:10]')
%!error n = islinphase ([1:10]', [1:10]')
%!error n = islinphase (1:10, 1:10, 1:10)
%!error n = islinphase (ones (3), ones (3))

%!test
%! a = 1; b = [1, 2, 3, 4, 3, 2, 1];
%! f = islinphase (b, a);
%! assert (f, true)

%!test
%! a = 1; b = [1, 2, 3, 4, 4, 3, 2, 1];
%! f = islinphase (b, a);
%! assert (f, true)

%!test
%! b = [1 -1]; a = [1 -2];
%! f = islinphase (b, a, 1E-1);
%! assert (f, false)

%!test
%! [b, a] = butter (1, .5); 
%! f = islinphase (b, a);
%! assert (f, true)

%!test
%! [b, a] = butter (1, .25); 
%! f = islinphase (b, a);
%! assert (f, false)

%!test
%! [b, a] = butter (8, .5);
%! f = islinphase (b, a);
%! assert (f, false)

