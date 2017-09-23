% MATLAB/octave test driver for nodal2dziff. Barnett 9/22/17

clear;
N = 3000; u = rand(N,N)-1/2;    % iid site perc on Z^2.
o.verb = 1;
o.sign = 1;

tic; [d nd siz ier]  = nodal2dziff(u,o); toc  % runs for <1 s

fprintf('nd=%d domains found\nbiggest few sizes:\n',nd)
ssiz = sort(siz,'descend'); ssiz(1:10)

figure; imagesc(d); axis xy equal tight; caxis([1 nd]); colorbar
