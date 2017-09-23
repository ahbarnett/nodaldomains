% MATLAB/octave test driver for perc2d. Barnett 9/22/17

clear;
N = 3000; u = rand(N,N)-1/2;    % iid site perc on Z^2.
o.verb = 1;

h0known = 1/2 - 0.592746;  % site perc thresh for Z^2, https://en.wikipedia.org/wiki/Percolation_threshold
hran = [-0.2,0]; nh = 1e5;
tic; [h0 d nd siz ier]  = perc2d(u,hran,nh,o); toc  % <1 s
fprintf('threshold h0=%g (known result for site perc on Z^2 = %g)\n',h0,h0known)
