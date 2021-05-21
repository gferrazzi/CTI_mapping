%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements: Electrical conductivity and permittivity maps of brain      %          
% tissues derived                                                         % 
% from water content based on T1-weighted acquisition. Michel E, Hernandez%
% D, Lee SY, MRM, 2017 77(3):1094-1103,                                   %
% https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.26193              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

thisFolder=pwd;
cd(thisFolder)
cd('../../')
thisFolder=pwd;

addpath(genpath('software/'))

niiFolder = [thisFolder '/data/nifties/T1/'];
filename  = 'TR700.nii';
hdr1= load_nii_hdr([niiFolder filename]);
V1=load_untouch_nii([niiFolder filename]);
Im1=double(V1.img)*(hdr1.dime.scl_slope);

filename  = 'TR3000.nii';
hdr2= load_nii_hdr([niiFolder filename]);
V2=load_untouch_nii([niiFolder filename]);
Im2=double(V2.img)*(hdr2.dime.scl_slope);

cd('processing/Michel_et_al/')
thisFolder = pwd;

[~,~,slices] = size(Im1);

%------------------------------------------------%
%                MODEL COEFFICIENTS              %
%------------------------------------------------%  

sliceToShow = 40; % change this parameter to visualize different slices
c1=0.286;
c2=1.526E-5;
c3=11.852;
w1=1.525;
w2=1.443;
         
CondRange = [0 2.5]; 
WaterRange = [60 100];
IntensRange = [100 1400];

%------------------------------------------------------%
%    Get ROI mask to identify the object of interest   %
%------------------------------------------------------%

porcenThld = 0.03;      
tmp1x=abs(Im2);
ROI=ones(size(tmp1x));
tmp2x = tmp1x;

 for sl=1:slices
        tmp1 = medfilt2(tmp1x(:,:,sl), [3 3]);
        H = fspecial('average', [3 3]);
        tmp1 = imfilter(tmp1,H,'replicate'); 
        tmp1x(:,:,sl) = tmp1x(:,:,sl)/max(max(tmp1x(:)));
        tmp2x(:,:,sl) = tmp1/max(tmp1(:));
 end
        
Thld= max(abs(tmp2x(:)))*porcenThld;                 
ROI(abs(tmp1x)<=Thld)=0;                            
obj_v=(ROI==1);                                     
clear tmp1; clear tmp2;  
clear tmp1x
clear tmp2x

%----------------------------%
%    Spin Echo image ratio   %
%----------------------------%

Im1(Im1==0)=eps;    
Im2(Im2==0)=eps;    
ImRatio = abs(Im1./Im2); 

ImRatio(isnan(ImRatio))=0; 
ImRatio(isinf(ImRatio))=0;
Ir = ImRatio;
IrC = abs((Ir).*ROI);
          
%----------------%
%    Water map   %
%----------------%

IW = w1.*exp(-w2.*IrC);  
IW = abs(IW.*ROI);
V1.img = abs(IW.*ROI*10000);
save_untouch_nii(V1,[thisFolder '/water_v1.nii']);
IW(IW > 1) = 1;
V1.img = abs(IW.*ROI*10000);
save_untouch_nii(V1,[thisFolder '/water_v2.nii']);

tmp = IW(:,end:-1:1,sliceToShow);
figure; set(gcf, 'defaultaxesfontsize', 14);
imagesc(abs(tmp'*100), [0 100]);
colormap gray;
axis image;
axis ij; 
msg =  strcat('Water Content Map (%)');
title(msg,'FontSize',16);
colorbar;
axis off                   
          
%-----------------------%
%    Conductivity map   %
%-----------------------%

ImCond = (c1+c2.*exp(c3.*IW));
ROIept = IW;
ROIept(ROIept~=0)=1;
V1.img = abs(ImCond.*ROIept*10000);
save_untouch_nii(V1,[thisFolder '/conductivity.nii']);

ImCond = ImCond.*ROIept;
tmp = abs(ImCond(:,end:-1:1,sliceToShow));
figure; set(gcf, 'defaultaxesfontsize', 14);
imagesc(tmp', CondRange);
axis image;
axis ij;
msg =  strcat('Conductivity Map (S/m)');
title(msg,'FontSize',16);
colorbar;
axis off         
