function [d nd siz ier] = nodal2dziff(u,o)
% NODAL2DZIFF  compute nodal domains (sets u>0, u<0) on 2D grid.
%
% function [d nd siz] = nodal2dziff(u)
% function [d nd siz ier] = nodal2dziff(u,o)
% Inputs:
%  u (N*N) = square array of real numbers.
%  o optional struct allowing one to set:
%     o.sign : if 0 count all domains (default); if >0 just +ve, <0 just -ve.
%     o.verb : verbosity (0,1,2,3)
% Outputs:
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
if nargin<2, o = []; end
verb = 0; if isfield(o,'verb'); verb = o.verb; end
sign = 0; if isfield(o,'sign'); sign = o.sign; end

mex_id_ = 'o int = nodal2dziff(i int, i double[], o int[xx], o int[x], o int[x], i int, i int)';
[ier, d, siz, nd] = gateway(mex_id_, N, u, sign, verb, N, N, n, 1);

siz = siz(1:nd)';        % truncate

% -------------------------------------------------------------------
