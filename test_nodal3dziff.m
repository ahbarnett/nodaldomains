% MATLAB/octave test driver for nodal3dziff. Barnett 8/30/17

clear;
N = 200; u = rand(N,N,N)-1/2;    % iid site perc on Z^3.
o.verb = 1;

tic; [d nd siz ier]  = nodal3dziff(u,o); toc  % runs for <1 s

fprintf('nd=%d domains found\nbiggest few sizes:\n',nd)
ssiz = sort(siz,'descend'); ssiz(1:10)  % expect 2 macroscopic, others small

% animate a slice...
figure; for z=1:5:N, imagesc(squeeze(d(:,:,z))); axis equal xy tight;
  title(sprintf('z=%d',z)); caxis([1 nd]); drawnow; end, colorbar
