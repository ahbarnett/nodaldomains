function [d nd siz ier] = nodal3dziff(u,o)
% NODAL3DZIFF  compute nodal domains (sets u>0, u<0) on 3D grid.
%
% function [d nd siz] = nodal3dziff(u)
% function [d nd siz ier] = nodal3dziff(u,o)
% Inputs:
%  u (N*N*N) = cube of real numbers.
%  o sets options such as o.verb : verbosity (0,1,...)
% Outputs:
%  d (N*N*N) = integer nodal domain labeling for sets u>0 and u<0.
%               Labels are 1...nd.
%  nd = number of nodal domains = max(d(:))
%  siz (1*nd) = sizes of nodal domains
%  ier = status: 0 success, >0 failure
%
% (c) Alex Barnett 8/31/17

N = size(u,1);
if size(u,2)~=N | size(u,3)~=N, error('u must be cubical'); end
n=N*N*N;
if nargin<2, o = []; end
verb = 0; if isfield(o,'verb'); verb = o.verb; end

mex_id_ = 'o int = nodal3dziff(i int, i double[], o int[x], o int[x], o int[x], i int)';
[ier, d, siz, nd] = gateway(mex_id_, N, u, verb, n, n, 1);

d = reshape(d,[N N N]);  % mwrap can't handle 3d arrays
siz = siz(1:nd)';        % truncate

% -------------------------------------------------------------------
