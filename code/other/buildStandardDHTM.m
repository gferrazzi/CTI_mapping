function [F,FH]=buildStandardDHTM(N,fn,gpu)

%BUILDSTANDARDDHTM   Builds standard DHT matrices
%   [F,FH]=BUILDSTANDARDDHTM(N,{FN},{GPU})
%   * N are the dimensions of the space
%   * {FN} indicates whether to generate fully unitary Hartley matrices. It
%   defaults to 0
%   * {GPU} is a flag that determines whether to generate gpu (1) or cpu
%   (0) matrices (empty, default depending on machine)
%   ** F is a cell of discrete Hartley transform matrices along the 
%   different dimensions
%   ** FH is a cell of inverse discrete Hartley transform matrices along 
%   the different dimensions
%

if nargin<2 || isempty(fn);fn=0;end
if nargin<3 || isempty(gpu);gpu=useGPU;end

ND=length(N);
F=cell(1,ND);FH=cell(1,ND);
for m=1:ND
    if nargout>1;[F{m},FH{m}]=build1DHTM(N(m),fn,gpu);
    else F{m}=build1DHTM(N(m),fn,gpu);
    end
end
