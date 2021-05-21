%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conductivity Tensor Imaging of the human brain using water mapping      %    
% techniques. Marino M, Cordero-Grande L, Mantini M, Ferrazzi G           %
% Frontiers of Neuroscince, In Press.                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

thisFolder=pwd;
cd(thisFolder)
cd('../../')
addpath(genpath('software/'))
addpath(genpath('code/'))
cd('processing/CTI/')
thisFolder = pwd;

%% Reconstruction parameters 
smoothing = 0.7;
sliceToShow = 40; 
zrange = 1:81; % modify this paramer to run on less slices (faster)
beta = 0.41;
randPerm = 16;

b3 = 1000;

system(['mkdir ' num2str(b3)])

%% DTI data
DTIi = load_untouch_nii('DTI.nii');
DTI  = double(DTIi.img);
bval = load('DTI.bval');
bval = bval';
GradientOrientations = load('DTI.bvec');

%% High frequency conductivity
CTIs = load_untouch_nii('conductivity_hf.nii');
conductivity  = (double((CTIs.img))/10000);

%% Mask
[XGO, YGO, ZGO, temp] = size(DTI);
mask = load_untouch_nii('MPRAGE_T1SPACE_BET_MASK.nii'); mask = double(mask.img);
mask = repmat(mask,[1 1 1 temp]);
DTI  = DTI.*mask;

%% Smoothing
for bi = 1 : length(bval)
   for Z = zrange
        vol = DTI(:,:,Z,bi);
        vol = imgaussfilt(vol,smoothing);
        DTI(:,:,Z,bi) = vol;
   end
end

%% Diffusion tensor estimation
DTIb3 = DTI(:,:,:,bval == b3);
bval_b3 = bval(bval == b3);
GradientOrientations_b3 = GradientOrientations(:,bval == b3);
GradientOrientations_b3 = GradientOrientations_b3';
DTIb3  = cat(4,DTI(:,:,:,1),DTIb3);

DTIi.hdr.dime.dim(5) = length(bval_b3)+1;
DTIi.img = DTIb3;
save_untouch_nii(DTIi,[num2str(b3) '/DTIb3.nii'])
bval_b3 = [0; bval_b3];
GradientOrientations_b3 = [0 0 0; GradientOrientations_b3];
dlmwrite([num2str(b3) '/bval_b3.txt'],bval_b3')
dlmwrite([num2str(b3) '/bvec_b3.txt'],GradientOrientations_b3');

% You may experience problems running these commands from Matlab. If so, copy and paste
% into terminal and execute this script in debug mode (folder 1000/)
system('mrconvert -fslgrad bvec_b3.txt bval_b3.txt DTIb3.nii DTIb3.mif -force');
system('dwi2tensor DTIb3.mif Db3.nii -force')

Db3 = load_untouch_nii([num2str(b3) '/Db3.nii']);
Db3 = Db3.img;
Db3(isnan(Db3)) = 0;
Db3 = double(Db3);

DTIi.hdr.dime.dim(5) = 6;
DTIi.img = Db3*10000000;
save_untouch_nii(DTIi,[num2str(b3) '/Db3.nii'])

%% Model fitting, eq 6, 7, 8, 9 and 10
[Xe,de,di,vic,viso,dest] = fittingFuntion_final(DTI,randPerm,bval,zrange,mask(:,:,:,1));

Xen = CTIs;
Xe(isinf(Xe)) = 0;
Xe(isnan(Xe)) = 0;
Xen.img = Xe*10000;
save_untouch_nii(Xen,[num2str(b3) '/Xe.nii'])

figure; set(gcf, 'defaultaxesfontsize', 14);
imagesc(Xe(:,end:-1:1,sliceToShow)',[0 1]);
axis image;
axis ij;
msg =  strcat('Xe [0 1]');
title(msg,'FontSize',16);
colorbar;
colormap(gray)
axis off       

den = CTIs;
de(isinf(de)) = 0;
de(isnan(de)) = 0;
den.img = de*100000;
save_untouch_nii(den,[num2str(b3) '/de.nii'])

figure; set(gcf, 'defaultaxesfontsize', 14);
imagesc(de(:,end:-1:1,sliceToShow)',[0 3*10^-3]);
axis image;
axis ij;
msg =  strcat('de (mm^2/s)');
title(msg,'FontSize',16);
colorbar;
colormap(gray)
axis off       

din = CTIs;
di(isinf(di)) = 0;
di(isnan(di)) = 0;
din.img = di*100000;
save_untouch_nii(din,[num2str(b3) '/di.nii'])

figure; set(gcf, 'defaultaxesfontsize', 14);
imagesc(di(:,end:-1:1,sliceToShow)',[0 3*10^-3]);
axis image;
axis ij;
msg =  strcat('di (mm^2/s)');
title(msg,'FontSize',16);
colorbar;
colormap(gray)
axis off 

Xen = CTIs;
vic(isinf(vic)) = 0;
vic(isnan(vic)) = 0;
Xen.img = vic*10000;
save_untouch_nii(Xen,[num2str(b3) '/vic.nii'])

figure; set(gcf, 'defaultaxesfontsize', 14);
imagesc(vic(:,end:-1:1,sliceToShow)',[0 1]);
axis image;
axis ij;
msg =  strcat('vic [0 1]');
title(msg,'FontSize',16);
colorbar;
colormap(gray)
axis off 

den = CTIs;
viso(isinf(viso)) = 0;
viso(isnan(viso)) = 0;
den.img = viso*10000;
save_untouch_nii(den,[num2str(b3) '/viso.nii'])

figure; set(gcf, 'defaultaxesfontsize', 14);
imagesc(viso(:,end:-1:1,sliceToShow)',[0 1]);
axis image;
axis ij;
msg =  strcat('viso [0 1]');
title(msg,'FontSize',16);
colorbar;
colormap(gray)
axis off 

din = CTIs;
dest(isinf(dest)) = 0;
dest(isnan(dest)) = 0;
din.img = dest*100000;
save_untouch_nii(din,[num2str(b3) '/dest.nii'])

figure; set(gcf, 'defaultaxesfontsize', 14);
imagesc(dest(:,end:-1:1,sliceToShow)',[0 1*10^-3]);
axis image;
axis ij;
msg =  strcat('d_e^* (mm^2/s)');
title(msg,'FontSize',16);
colorbar;
colormap(gray)
axis off 

%% Extracellular diffusion tensor calculation, eq 11 and 12
Db3SVD = zeros(XGO,YGO,ZGO,3);
for z = zrange
    for x = 1 : XGO
        for y = 1 : YGO
            if mask(x,y,z) == 1             
                
                Tensor = [Db3(x,y,z,1) Db3(x,y,z,4) Db3(x,y,z,5)
                    Db3(x,y,z,4) Db3(x,y,z,2) Db3(x,y,z,6)
                    Db3(x,y,z,5) Db3(x,y,z,6) Db3(x,y,z,3)];          
                [U,S,V] = svd(Tensor);
                S = diag(S);
                Db3SVD(x,y,z,1) = S(1);
                Db3SVD(x,y,z,2) = S(2);
                Db3SVD(x,y,z,3) = S(3);
                
            end
        end
    end
end
Db3SVD = sum(Db3SVD,4);
ni = (3.*de)./(Db3SVD);

DTIi.hdr.dime.dim(5) = 1;
DTIi.img = ni*10000;
save_untouch_nii(DTIi,[num2str(b3) '/ni.nii'])

%% CTI calculation, eq 1

ce = (conductivity)./(Xe.*de + (1-Xe).*di.*beta);
ce(isinf(ce)) = 0;
ce(isnan(ce)) = 0;

CTI = zeros(XGO,YGO,ZGO,6);
CTI(:,:,:,1) = Xe.*ni.*ce.*Db3(:,:,:,1);
CTI(:,:,:,2) = Xe.*ni.*ce.*Db3(:,:,:,2);
CTI(:,:,:,3) = Xe.*ni.*ce.*Db3(:,:,:,3);
CTI(:,:,:,4) = Xe.*ni.*ce.*Db3(:,:,:,4);
CTI(:,:,:,5) = Xe.*ni.*ce.*Db3(:,:,:,5);
CTI(:,:,:,6) = Xe.*ni.*ce.*Db3(:,:,:,6);

DTIi.hdr.dime.dim(5) = 6;
DTIi.img = CTI*10000;
save_untouch_nii(DTIi,[num2str(b3) '/CTItensor.nii']);

CTIcomp = permute(CTI(:,end:-1:1,sliceToShow,:),[2 1 3 4]);
CTIcomp = [CTIcomp(:,:,1,1), CTIcomp(:,:,1,4), CTIcomp(:,:,1,5)
    CTIcomp(:,:,1,4), CTIcomp(:,:,1,2), CTIcomp(:,:,1,6)
    CTIcomp(:,:,1,5), CTIcomp(:,:,1,4), CTIcomp(:,:,1,3)];

figure; set(gcf, 'defaultaxesfontsize', 14);
imagesc(CTIcomp,[0 1]);
axis image;
axis ij;
msg =  strcat('CTI (S/m)');
title(msg,'FontSize',16);
colorbar;
colormap(magma)
axis off 