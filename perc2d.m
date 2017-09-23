function [h0 d nd siz ier] = perc2d(u,hran,nh,o)
% PERC2D  find percolation threshold and domains for set u>h on 2D grid.
%
% function [h0 d nd siz ier] = perc2d(u,hran,nh,o)
% Inputs:
%  u (N*N) = square array of real numbers.
%  hran (length-2) = lower and upper limit of h range to explore
%  nh = number of uniformly spaced h values in range. (in empty, uses 1e5)
%  o sets options such as o.verb : verbosity (0,1,...)
% Outputs:
%  h0 = estimate of lowest h where the set u>h percolates (across z faces).
%  d (N*N) = integer nodal domain labeling for sets u>0 and u<0.
%            Labels are 1...nd.
%  nd = number of nodal domains = max(d(:))
%  siz (1*nd) = sizes of nodal domains
%  ier = status: 0 success, >0 failure
%
% (c) Alex Barnett 9/22/17

N = size(u,1);
if size(u,2)~=N, error('u must be square'); end
n=N*N;
if nargin<3 | isempty(nh), nh=1e5; end
if nargin<4, o = []; end
verb = 0; if isfield(o,'verb'); verb = o.verb; end

mex_id_ = 'o int = perc2d(i int, i double[], o int[xx], o int[x], o int[x], i int, i double[x], o double[x], i int)';
[ier, d, siz, nd, h0] = gateway(mex_id_, N, u, nh, hran, verb, N, N, n, 1, 2, 1);

siz = siz(1:nd)';        % truncate

% -------------------------------------------------------------------
