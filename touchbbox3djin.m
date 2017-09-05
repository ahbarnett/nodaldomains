% TOUCHBBOX3DJIN   find touching domains and bounding box for 3D domain labels
%
% function [domaindata]=touchbbox3djin(A,l)
%
% Input A is a 3D matrix of nodal domain labels, which must be zero-padded on
%  all faces. l is the label of the domain of interest.
% Output is cell array containing entries:
%  1) whether or not touches the border (-1 if touches, 1 if not)
%  2) the neighbors of the nodal domain, ie domains that touch it.
%  3) the other domains immediately contained within this
%     nodal domain (does not account for domains within "contained" domains)
%  4-9) the indices defining the minimal-sized bounding box for 
%     the nodal domain (min row, max row, mincol, maxcol, minplane, maxplane),
%     1-indexed (including the zero-padding), and inclusive. Eg in 1d the domain
%     number 1 in the list [0 0 1 1 0] would have minrow = 2, maxrow = 5.
%
% Without inputs, a self-test is done.
%
% Matthew Jin, Aug 2014 ("findContained.m")
% Loop removed; old code removed, doc & self-test added. Barnett 9/4/17

function [data]=touchbbox3djin(A,l)

if nargin==0, test_touchbbox3djin; return, end

data=cell(1,9);
% matrix where 1s represent entry of domain in question, 0's otherwise
temp=(A==l);

[width,length,depth]=ind2sub(size(temp),find(temp==1));
% locate entries of a bounding box to reduce computation time
minrow=min(width)-1;       
maxrow=max(width)+1;
mincol=min(length)-1;      
maxcol=max(length)+1;
minplane=min(depth)-1;
maxplane=max(depth)+1;

% reducing temp to only necessary size reduces overall computation time
temp=temp(minrow:maxrow,mincol:maxcol,minplane:maxplane);

% temp2 stores the original domain numbers
temp2=A(minrow:maxrow,mincol:maxcol,minplane:maxplane);

% down (y direction), n for neighbor, b for boundary

[i,j,k]=ind2sub(size(temp(1:end-1,:,:)),find(temp(1:end-1,:,:).*~temp(2:end,:,:)==1));

nposy=[i+1,j,k,zeros(numel(i),1)];

% up (y direction)
[i,j,k]=ind2sub(size(temp(1:end-1,:,:)),find(~temp(1:end-1,:,:).*temp(2:end,:,:)==1));
nnegy=[i,j,k,zeros(numel(i),1)];

% right (x direction)
[i,j,k]=ind2sub(size(temp(:,1:end-1,:)),find(temp(:,1:end-1,:).*~temp(:,2:end,:)==1));
nposx=[i,j+1,k,zeros(numel(i),1)];

% left (x direction)
[i,j,k]=ind2sub(size(temp(:,1:end-1,:)),find(~temp(:,1:end-1,:).*temp(:,2:end,:)==1));
nnegx=[i,j,k,zeros(numel(i),1)];

% above (z direction)
[i,j,k]=ind2sub(size(temp(:,:,1:end-1)),find(~temp(:,:,1:end-1).*temp(:,:,2:end)==1));
nnegz=[i,j,k,zeros(numel(i),1)];

% below (z direction)
[i,j,k]=ind2sub(size(temp(:,:,1:end-1)),find(temp(:,:,1:end-1).*~temp(:,:,2:end)==1));
nposz=[i,j,k+1,zeros(numel(i),1)];
        
% append the neighboring nodal domain numbers to their corresponding
% indices
nposyind=sub2ind(size(temp),nposy(:,1),nposy(:,2),nposy(:,3));
nposy(:,4)=temp2(nposyind);

nnegyind=sub2ind(size(temp),nnegy(:,1),nnegy(:,2),nnegy(:,3));
nnegy(:,4)=temp2(nnegyind);

nposxind=sub2ind(size(temp),nposx(:,1),nposx(:,2),nposx(:,3));
nposx(:,4)=temp2(nposxind);

nnegxind=sub2ind(size(temp),nnegx(:,1),nnegx(:,2),nnegx(:,3));
nnegx(:,4)=temp2(nnegxind);

nposzind=sub2ind(size(temp),nposz(:,1),nposz(:,2),nposz(:,3));
nposz(:,4)=temp2(nposzind);

nnegzind=sub2ind(size(temp),nnegz(:,1),nnegz(:,2),nnegz(:,3));
nnegz(:,4)=temp2(nnegzind);

% list of all faces
facelist=[nposy; nnegy; nposx; nnegx; nposz; nnegz];

neighbors=unique(facelist(:,4))';

if sum(neighbors==0)~=0;
  data{1}=-1; % -1 for border-touching domains
else
  data{1}=1;
end
[boxrows,boxcols,boxplanes]=size(temp); % dimensions to be compared to

contained=NaN(numel(neighbors),1);
for c=1:numel(neighbors)
  [nrows,ncols,nplanes]=ind2sub(size(temp2),find(temp2==neighbors(c)));

  % neighbor bounding indices
  nminrow=min(nrows);
  nmaxrow=max(nrows);
  nmincol=min(ncols);
  nmaxcol=max(ncols);
  nminplane=min(nplanes);
  nmaxplane=max(nplanes);

  % in order for neighbor to be contained in the domain
  % all of its bounding indices must lie within those
  % of the domain. Note that we say >2 and < DIM-1 because the
  % domain has been "padded"
  cond=nminrow>2 && nmincol>2 && nminplane>2 && nmaxrow<boxrows-1 && nmaxcol<boxcols-1 && nmaxplane<boxplanes-1;
            
  if cond
    contained(c)=1; % neighbor is contained within domain
  else
    contained(c)=0;
  end
            
end
contained=neighbors(contained==1);
data{2}=neighbors;
data{3}=contained;
data{4}=minrow;
data{5}=maxrow;
data{6}=mincol;
data{7}=maxcol;
data{8}=minplane;
data{9}=maxplane;

%%%%%%%%%%%%%%%%%%%%%%%%
function test_touchbbox3djin
N=100; x=(-(N-1)/2:(N-1)/2)/(N/2);  % x grid in [-1,1]
[xx yy zz] = ndgrid(x,x,x);
rr = xx.^2+yy.^2+zz.^2;
r0 = 0.5;
u = 1+(rr<r0^2);   % sphere of 2 surrounded by 1.
A = zeros(N+2,N+2,N+2); A(2:N+1,2:N+1,2:N+1)=u;  % zero pad thickness 1
c = touchbbox3djin(A,1), c{2}
% should be:
%c =
%  1×9 cell array
%    [-1]    [1×2 double]    [2]    [1]    [102]    [1]    [102]    [1]    [102]
%ans =
%     0     2
c = touchbbox3djin(A,2)
% should be:
%c =
%  1×9 cell array
%    [1]    [1]    []    [26]    [77]    [26]    [77]    [26]    [77]
