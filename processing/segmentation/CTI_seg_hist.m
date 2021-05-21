%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     CREATE FIGURE 7 SUBPANEL of                         %
% Conductivity Tensor Imaging of the human brain using water mapping      %    
% techniques. Marino M, Cordero-Grande L, Mantini M, Ferrazzi G           %
% Frontiers of Neuroscince, In Press.                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all, clc

thisFolder=pwd;
cd(thisFolder)
cd('../../')
addpath(genpath('software/'))
addpath(genpath('code/'))
cd('processing/segmentation/')
thisFolder = pwd;

%% param
sliceToShow = 1:81;
color       = viridis;

%% load data
cti_filename = 'CTItensor.nii';
anatomy_CSF_filename = 'csf.nii';
MichelWaterMap    = '../Michel_et_al/water_v1.nii';
MichelCondMap     = '../Michel_et_al/conductivity.nii';
anatomy_GM_filename = 'MPRAGE_T1SPACE_brain_pve_1.nii';
anatomy_WM_filename = 'MPRAGE_T1SPACE_brain_pve_2.nii';
brain_mask_filename = 'mask.nii';
isotropy_filename = 'CTIisotropic.nii';
ratioCond         = 'ratioT.nii';

CMmap = load_untouch_nii(MichelCondMap);
CMmap = CMmap.img;
Wmap = load_untouch_nii(MichelWaterMap);
Wmap = Wmap.img;
Wmap(Wmap<=10000)=0;
Wmap(Wmap>0)=1;
CTI = load_untouch_nii(cti_filename);
CTI = CTI.img;
anatomy_CSF = load_untouch_nii(anatomy_CSF_filename);
CSF_mask = anatomy_CSF.img;
if not(strcmp(anatomy_CSF_filename,'csf.nii'))
    CSF_mask(CSF_mask<0.99) = 0;
end
GM_mask = load_untouch_nii(anatomy_GM_filename);
GM_mask = GM_mask.img;
GM_mask(GM_mask<0.99) = 0;
WM_mask = load_untouch_nii(anatomy_WM_filename);
WM_mask = WM_mask.img;
WM_mask(WM_mask<0.99) = 0;
brain  = load_untouch_nii(brain_mask_filename);
brain  = brain.img;

%% SVD on CTI map to calculate sigma_iso
sigma11 = CTI(:,:,:,1);
sigma12 = CTI(:,:,:,4);
sigma13 = CTI(:,:,:,5);
sigma22 = CTI(:,:,:,2);
sigma23 = CTI(:,:,:,6);
sigma33 = CTI(:,:,:,3);

xdim = size(CSF_mask,1);
ydim = size(CSF_mask,2);
zdim = size(CSF_mask,3);

CSF_VOL_tot_ani = mean(CTI(:,:,:,[1 2 3]),4);
VOL_tot_ani = zeros(size(CSF_VOL_tot_ani));
for ix = 1:xdim
    for iy = 1:ydim
        for iz = 1:zdim

            sig_xx = sigma11(ix,iy,iz);
            sig_xy = sigma12(ix,iy,iz);
            sig_yx = sig_xy;
            sig_xz = sigma13(ix,iy,iz);
            sig_zx = sig_xz;
            sig_yy = sigma22(ix,iy,iz);
            sig_yz = sigma23(ix,iy,iz);
            sig_zy = sig_yz;
            sig_zz = sigma33(ix,iy,iz);
            
            C = [sig_xx sig_xy sig_xz;
                sig_yx sig_yy sig_yz;
                sig_zx sig_zy sig_zz];
            
            C = double(C);
            
            [U,S,V] = svd(C);
            
            main_axes = diag(S);
            
            VOL_ani=(main_axes(1))*(main_axes(2))*(main_axes(3));

            if isinf(VOL_ani)
                VOL_ani = 0;
            end
            
            VOL_tot_ani(ix,iy,iz) = (VOL_ani)^(1/3);
                       
        end
    end
end

V = anatomy_CSF;
V.img = VOL_tot_ani;
save_untouch_nii(V,isotropy_filename);

%% GM, WM, CSF sigma_iso
Wmap = double(Wmap);
VOL_tot_ani_GM = (VOL_tot_ani.*GM_mask)/10000;
VOL_tot_ani_WM = (VOL_tot_ani.*WM_mask)/10000;
VOL_tot_ani_CSF = (VOL_tot_ani.*CSF_mask.*(1-Wmap))/10000;

VOL_tot_ani_GM  = VOL_tot_ani_GM(:,:,sliceToShow);
VOL_tot_ani_WM  = VOL_tot_ani_WM(:,:,sliceToShow);
VOL_tot_ani_CSF = VOL_tot_ani_CSF(:,:,sliceToShow);

X_GM = VOL_tot_ani_GM((VOL_tot_ani_GM>0));
X_WM = VOL_tot_ani_WM((VOL_tot_ani_WM>0));
X_CSF = VOL_tot_ani_CSF((VOL_tot_ani_CSF>0));

%% FIGURE 7
figure
h1 = histogram(X_GM,200,'Normalization','probability','FaceColor',color(1,:));
hold on
h2 = histogram(X_WM,200,'Normalization','probability','FaceColor',color(128,:));
hold on
h3 = histogram(X_CSF,100,'Normalization','probability','FaceColor',color(230,:));
set(gcf,'color','w');
legend('\phi_{GM}','\phi_{WM}','\phi_{CSF}','Location','northeast')
title('LF conductivity','Fontsize',12)
set(gca,'FontSize',15)
grid on
xlabel('S/m')
ylabel('Normalized frequency')
axis([0 2.5 0 0.2])
saveas(gcf,'distribution_all_range.pdf')