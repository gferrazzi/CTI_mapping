%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     CREATE FIGURE 5                                     %
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
cd('processing/segmentation/')
thisFolder = pwd;

condM = load_untouch_nii('../Michel_et_al/conductivity.nii');
condCTI = load_untouch_nii('CTIisotropic.nii');
 
slice = 40;

yimg1 = double(condM.img);
yimg2 = double(condCTI.img);

yimg1 = yimg1(:,:,slice);
yimg2 = yimg2(:,:,slice);

mask  = yimg2>0;
yimg1 = yimg1.*mask;
yimg2 = yimg2.*mask;

yimg1 = (yimg1./(mean(yimg1(:))))*1-1; % demean
yimg2 = (yimg2./(mean(yimg2(:))))*1-1;

mask  = not(yimg2==-1);
yimg1 = yimg1.*mask;
yimg2 = yimg2.*mask;

N=[size(yimg1), 1];%Size
y{1}=yimg1;
y{2}=yimg2;

figure
hold on
for s=2:2%Readout/Phase encode/Slice (use slice with discretion in MRI...)
    for t=1:2%First/Second image
        [x,k]=estimatePSD(y{t},s,1,'dB');
        x=ifftshift(x);k=ifftshift(k);
        x=x(1:end/2);k=k(1:end/2)*2;        
        x=x-x(1);
        plot(k,x,'LineWidth',2)  
        hold on
    end   
end

ylim([-40 5])
xl=xticklabels;
for l=1:length(xl);xl{l}=num2str(1/str2num(xl{l}));end
xticklabels(xl);   
grid on  
set(gcf,'color','w');
title('PSD','Fontsize',12)
set(gca,'FontSize',15)
grid on
xlabel('mm');
ylabel('dB');
legend('\sigma_{HF}','\sigma_{LF}^{iso}', '')    
saveas(gcf,'PSD.pdf')
