% GENUS3DJIN  compute genus of all 3d gridded domains and whether touch box.
%
% function [out g vol nei ins]=genus3djin(A)
% function [out ...]=genus3djin(A,numdoms)
%
% The input to this function is a 3D matrix A of nodal domain numbers, labeled
%  1 to numdoms. numdoms is optional input. A need not be zero-padded.
%  Outputs:
%  out - cell array of length numdoms.
%   The ith entry in the cell array contains the data for the ith nodal
%   domain in the following format: {genus, volume, neighbors}. Genus is the
%   calculated genus for the outer surface of nodal domain i, volume is the 
%   number of elements of A that have value i, and neighbors is a list of 
%   nodal domains that are orthogonally adjacent to at least one entry of
%   domain i.
%  g - list of genuses for all domains (1...numdoms)
%  vol - volumes (in voxels) for all domains
%  nei - cell array of neighbor labels for each domain (incl 0 if touches box)
%  ins - boolean list of whether each domain is inside the box (0 if touching
%        any box wall)
%  [Note: out is a different packaging of the other outputs]
%
% function [out ...]=genus3djin(A,numdoms,opts) also set the following options:
%  opts.corr=0 : doesn't correct for ambiguous 2x2x2 arrangements (default).
%            1 : handle ambiguous arrangements such as [0 1; 1 0] and a 2x2x2
%                block that has 6 1's in it, with the two 0's lying at
%                diagonally opposite corners in 3D space. 
%
% Without inputs, a self-test is done.
%
% See also: TOUCHBBOX3DJIN (was: findContained)
%
% Matthew Jin 8/26/14 ("genus3Dv4.m"). Tidy up & doc & test, Barnett 9/4/17

function [distribution g vol nei ins]=genus3djin(A,numdoms,opts)

if nargin==0, test_genus3djin; return, end
if nargin<2, numdoms = max(A(:)); end
if nargin<3, o=[]; end
if ~isfield(opts,'corr'), opts.corr=0; end   % default

[m,n,o]=size(A);
    
% pad with 0's for counting V,F,E later, and for touchbbox3djin...
tempmat=zeros(m+2,n+2,o+2);
tempmat(2:end-1,2:end-1,2:end-1)=A;
A=tempmat; clear tempmat;
    
% find bounding boxes and which domains touched by each domain...
% (slow: takes about 1/2 the time)
t=tic; database = cell(1,numdoms);
for l=1:numdoms
  database{l} = touchbbox3djin(A,l);
end
fprintf('genus3djin: touchbbox3djin (aka findContained) done, %g s\n',toc(t))

distribution=cell(1,numdoms);

% now must update each domain's contained list to include not just
% those immediately interior to domain, but also those which are within
% contained domains, and so on...    
updated=cell(1,numdoms); % updated is a dummy storage variable
for i=1:numdoms
  uncounted=database{i}{3};
  newcontained=database{i}{3}; % updated contained list
  
  while ~isempty(uncounted)
    additions=database{uncounted(1)}{3};  % contained in the contained
    newcontained=[newcontained,additions]; % add in such domains
    uncounted=[uncounted,additions]; % addition domains must also be counted
    uncounted(1)=[]; % entry is now counted 
  end    
  updated{i}=newcontained;
end
    
for i=1:numdoms
  database{i}{3}=[updated{i},i]; % domain i contains itself
end

% create the solid domain from data in database, and compute its genus:

for i=1:numdoms   % --------------- loop over domains
       
  nodaldata=cell(1,4);
  cutA=A(database{i}{4}:database{i}{5},database{i}{6}:database{i}{7},database{i}{8}:database{i}{9});
  
  sizemat=(cutA==i);
  volume=sum(sum(sum(sizemat)));
  
  soliddomain=ismember(cutA,database{i}{3});
  
  vertmat=soliddomain(1:end-1,1:end-1,1:end-1)+soliddomain(1:end-1,1:end-1,2:end)+soliddomain(1:end-1,2:end,1:end-1)+soliddomain(1:end-1,2:end,2:end)+soliddomain(2:end,1:end-1,1:end-1)+soliddomain(2:end,1:end-1,2:end)+soliddomain(2:end,2:end,1:end-1)+soliddomain(2:end,2:end,2:end);

  % compute V, F, E
  
  xedgemat=soliddomain(:,1:end-1,1:end-1)+soliddomain(:,1:end-1,2:end)+soliddomain(:,2:end,1:end-1)+soliddomain(:,2:end,2:end);
  yedgemat=soliddomain(1:end-1,:,1:end-1)+soliddomain(2:end,:,1:end-1)+soliddomain(1:end-1,:,2:end)+soliddomain(2:end,:,2:end);
  zedgemat=soliddomain(1:end-1,1:end-1,:)+soliddomain(2:end,1:end-1,:)+soliddomain(1:end-1,2:end,:)+soliddomain(2:end,2:end,:);
  
  xedgemat=(xedgemat<4 & xedgemat>0);
  yedgemat=(yedgemat<4 & yedgemat>0);
  zedgemat=(zedgemat<4 & zedgemat>0);
  
  E=sum(sum(sum(xedgemat)))+sum(sum(sum(yedgemat)))+sum(sum(sum(zedgemat)));
  
  if opts.corr       % Jin's v4 correction
    [x1,x2,x3]=ind2sub(size(xedgemat),find(xedgemat==2));
    
    ind1=find(soliddomain(sub2ind(size(soliddomain),x1,x2,x3))==soliddomain(sub2ind(size(soliddomain),x1,x2+1,x3+1)));
    
    f1=x1(ind1);
    g1=x2(ind1);
    h1=x3(ind1);
    
    xcorr=numel(f1);
    
    [y1,y2,y3]=ind2sub(size(yedgemat),find(yedgemat==2));
    
    ind2=find(soliddomain(sub2ind(size(soliddomain),y1,y2,y3))==soliddomain(sub2ind(size(soliddomain),y1+1,y2,y3+1)));
    
    f2=y1(ind2);
    g2=y2(ind2);
    h2=y3(ind2);
    
    ycorr=numel(f2);
    
    [z1,z2,z3]=ind2sub(size(zedgemat),find(zedgemat==2));
    
    ind3=find(soliddomain(sub2ind(size(soliddomain),z1,z2,z3))==soliddomain(sub2ind(size(soliddomain),z1+1,z2+1,z3)));
    
    f3=z1(ind3);
    g3=z2(ind3);
    h3=z3(ind3);
    
    zcorr=numel(f3);        
      
    E = E + xcorr+ycorr+zcorr;
                
    [i1,j1,k1]=ind2sub(size(vertmat),find(vertmat==6));
        
    f=(soliddomain==0);
        
    % vcorr 1-4 accounts for one type of problematic arrangment, where
    % opposite corners of a 2x2x2 block are removed (value 0)
        
    vcorr1=sum(f(sub2ind(size(f),i1,j1,k1)).*f(sub2ind(size(f),i1+1,j1+1,k1+1)));
    vcorr2=sum(f(sub2ind(size(f),i1,j1,k1+1)).*f(sub2ind(size(f),i1+1,j1+1,k1)));
    vcorr3=sum(f(sub2ind(size(f),i1,j1+1,k1)).*f(sub2ind(size(f),i1+1,j1,k1+1)));
    vcorr4=sum(f(sub2ind(size(f),i1,j1+1,k1+1)).*f(sub2ind(size(f),i1+1,j1,k1)));
    
    % this is accounting for a different problematic arrangement: count
    % extra vertices for the diagonal planar arrangement, i.e. [1 0; 0 1]
    
    vcorr5=sum(vertmat(sub2ind(size(vertmat),f1-1,g1,h1))==2);
    vcorr6=sum(vertmat(sub2ind(size(vertmat),f1,g1,h1))==2);
    vcorr7=sum(vertmat(sub2ind(size(vertmat),f2,g2-1,h2))==2);
    vcorr8=sum(vertmat(sub2ind(size(vertmat),f2,g2,h2))==2);
    vcorr9=sum(vertmat(sub2ind(size(vertmat),f3,g3,h3-1))==2);
    vcorr10=sum(vertmat(sub2ind(size(vertmat),f3,g3,h3))==2);
    
    vcorr=vcorr1+vcorr2+vcorr3+vcorr4+vcorr5+vcorr6+vcorr7+vcorr8+vcorr9+vcorr10;
  else, vcorr = 0;
  end
        
  vertmat=(vertmat<8 & vertmat>0);
  V=sum(sum(sum(vertmat)))+vcorr;
  
  % redo face sum for solid domain
  
  % down (y direction)
  nposy=numel(find(soliddomain(1:end-1,:,:).*~soliddomain(2:end,:,:)==1));
  
  % up (y direction)
  nnegy=numel(find(~soliddomain(1:end-1,:,:).*soliddomain(2:end,:,:)==1));
  
  % right (x direction)
  nposx=numel(find(soliddomain(:,1:end-1,:).*~soliddomain(:,2:end,:)==1));
  
  % left (x direction)
  nnegx=numel(find(~soliddomain(:,1:end-1,:).*soliddomain(:,2:end,:)==1));
  
  % above (z direction)
  nnegz=numel(find(~soliddomain(:,:,1:end-1).*soliddomain(:,:,2:end)==1));
  
  % below (z direction)
  nposz=numel(find(soliddomain(:,:,1:end-1).*~soliddomain(:,:,2:end)==1));
  
  F=sum(nposy+nposx+nposz+nnegy+nnegx+nnegz);
  
  genus=(2-V+E-F)/2;    % Euler
  
  % package for output...
  nodaldata{1}=genus;
  nodaldata{2}=volume;
  nodaldata{3}=database{i}{2}; % neighbor list
  nodaldata{4}=database{i}{1}; % -1 if touches box, 1 if not
  
  distribution{i}=nodaldata;
  
end    % ---------------------- end domains loop
fprintf('genus3djin finished, total time %g s\n',toc(t))

% Barnett repackage outputs easier...
g = []; vol = []; ins = []; nei = cell(1,numdoms);
for l=1:numdoms
  d = distribution{l};
  g = [g d{1}];
  vol = [vol d{2}];
  nei{l} = d{3};
  ins = [ins, d{4}==1];
end

%%%%%%%%%%%%%%%%%%%%%%%
function test_genus3djin
N=100; x=(-(N-1)/2:(N-1)/2)/(N/2);  % x grid in [-1,1]
[xx yy zz] = ndgrid(x,x,x);
o.corr = 0;

rr = xx.^2+yy.^2+zz.^2; r0 = 0.5; d = 1+(rr<r0^2); nd = 2;
'sphere of 2 surrounded by 1:'
[~,g,vol,nei,ins]=genus3djin(d,nd,o);%, 'nei:', nei{:}
if g==[0 0], 'genus ok', end

r1 = 0.25; d(rr<r1^2) = 3; nd = 3;
'...and with concentric smaller sphere of 3:'
[~,g,vol,nei,ins]=genus3djin(d,nd,o);%, 'nei:', nei{:}
if g==[0 0 0], 'genus ok', end

rho = sqrt(xx.^2+yy.^2); nd = 2;
'torus of 2 surrounded by 1:'
d = 1+((rho-0.5).^2+zz.^2<r1^2);
[~,g,vol,nei,ins]=genus3djin(d,nd,o);%, 'nei:', nei{:}
if g==[0 1], 'genus ok', end

sc = 12.5; x1=mod((1+xx)*sc/2,1)*2-1; y1=mod((1+yy)*sc/2,1)*2-1;
a = 0.5;
d = 1+(zz.^2<0.8-xx.^2-yy.^2 & abs(yy)<a & abs(xx)<a & x1.^2+y1.^2>.5^2);
nd = 2;
figure; show3ddomain(d,2,0);
'nasty shape of 2 with 36 holes in it surrounded by 1:'
[~,g,vol,nei,ins]=genus3djin(d,nd,o);%, 'nei:', nei{:}
if g==[0 36], 'genus ok', end
