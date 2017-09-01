% MATLAB/octave test driver for perc3d. Barnett 9/1/17

clear;
N = 200; u = rand(N,N,N)-1/2;    % iid site perc on Z^3.
o.verb = 1;

h0known = 1/2 - 0.3116076;  % site perc thresh for Z3 cubic, from wikipedia.
hran = [0,0.5]; nh = 1e5;
tic; [h0 d nd siz ier]  = perc3d(u,hran,nh,o); toc  % <1 s
fprintf('threshold h0=%g\n',h0)
