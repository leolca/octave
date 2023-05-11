## Copyright (C) 2012 Leonardo Araujo <leolca@gmail.com>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {@var{qidx} =} dpcmenco (@var{sig}, @var{codebook}, @var{partition}, @var{predictor})
## @deftypefnx {Function File} {[@var{qidx}, @var{q}] =} dpcmenco (@var{sig}, @var{codebook}, @var{partition}, @var{predictor})
## @deftypefnx {Function File} {[@var{qidx}, @var{q}, @var{d}] =} dpcmenco (@dots{})
## Encode using differential pulse code modulation (DPCM).
##
## @table @code
## @item qidx = dpcmenco (sig, codebook, partition, predictor)
## Determine position of the prediction error in a strictly monotonic table (partition).
## The predictor vector describes a m-th order prediction for the
## output according to the following equation
## y(k) = p(1)sig(k-1) + p(2)sig(k-2) + ... + p(m-1)sig(k-m+1) + p(m)sig(k-m) ,
## where the predictor vector is given by
## predictor = [0, p(1), p(2), p(3),..., p(m-1), p(m)].
##
## @item [qidx, q] = dpcmenco (sig, codebook, partition, predictor)
## Also return the quantized values.
##
## @item [qidx, q, d] = dpcmenco (...)
## Also compute distortion: mean squared distance of original sig from the
## corresponding quantized values.
##
## @end table
## @seealso{dpcmdeco, dpcmopt, quantiz}
## @end deftypefn

function [indx, quants, distor] = dpcmenco (sig, codebook, partition, predictor)

  if (nargin != 4)
    print_usage ();
  endif

  y = zeros (size (sig));
  indx = []; quants = [];
  for i = 1:length (y)
    ## use last predicted value to find the error
    y(i) = y(max (i - length (predictor) + 1, 1):i) * predictor(end:-1:max (end-i+1, 1))'; # convolution
    e(i) = sig(i) - y(i); # error
    [indx(i), quants(i)] = quantiz (e(i), partition, codebook); # quantize the error
    y(i) +=  quants(i); # update prediction value
  endfor

  ## compute distortion
  if (nargout > 2)
    sigq = dpcmdeco (indx, codebook, predictor);
    distor = sumsq (sig(:) - sigq(:)) / length (sig);
  endif

endfunction

%!function y = my_sawtooth (t, width)
%!  ## sawtooth function, so not to require the signal package for the demo and test
%!    if (nargin == 1)
%!      width = 1;
%!  endif
%!  t = mod (t / (2 * pi), 1);
%!  y = zeros (size (t));
%!  if (width != 0)
%!      y (t < width) = 2 * t (t < width) / width - 1;
%!  endif
%!  if (width != 1)
%!    y( t >= width) = -2 * (t (t >= width) - width) / (1 - width) + 1;
%!  endif
%!endfunction

%!demo
%! predictor = [0 1];
%! nbits = 4;
%! delta = 2^(-nbits+1);
%! codebook = [-1+delta/2 : delta : 1-delta/2];
%! partition = (codebook(1:end-1) + codebook(2:end))/2;
%! t = linspace (0, 2*pi, 128);
%! x = my_sawtooth (2*t, 0.25);
%! [idx, xq, distor] = dpcmenco (x, codebook, partition, predictor);
%! xr = dpcmdeco (idx, codebook, predictor);
%! plot (t, x, 'k--','linewidth',1, t, xr, 'b-','linewidth',1);
%! xlim ([0 2*pi]);
%! xlabel ('t'); ylabel ('x(t)');
%! legend ('original', 'dpcm');
%! title ( sprintf ('distortion = %.3f', distor) );

%% Test input validation
%!error dpcmenco ()
%!error dpcmenco (1)
%!error dpcmenco (1, 2)
%!error dpcmenco (1, 2, 3)
%!error dpcmenco (1, 2, 3, 4, 5)

%!test
%! predictor = [0 1];
%! partition = [-0.5, 0, 0.5];
%! codebook  = [-1, -0.25, 0.25, 1];
%! t = [0:pi/2:2*pi];
%! x = my_sawtooth(2*t);
%! [idx, qx, distor] = dpcmenco(x, codebook, partition, predictor);
%! assert (idx, [0, 3, 0, 3, 0])
%! assert (qx, [-1, 1, -1, 1, -1])
%! assert (distor, 0)
