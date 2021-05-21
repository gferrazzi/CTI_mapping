N=[64 64 32];%Size
y{1}=randn(N);%+1i*randn(N);%Random image (best resolution, flat spectrum)
MS=[1 1 2];%Resolution RO/PE/SL
H=buildFilter(N(1:2),'tukey',[1 1],0,[0.2 0.4]);%We use a Gibbs ringing filter with stronger suppression in PE than in RO [0.2 0.4] parameter, we don't filter in the SL direction
y{2}=filtering(y{1},H);%We filter the original image

figure
hold on
for s=1:3%Readout/Phase encode/Slice (use slice with discretion in MRI...)
    for t=1:2%First/Second image
        [x,k]=estimatePSD(y{t},s,MS(s),'dB');
        x=ifftshift(x);k=ifftshift(k);
        x=x(1:end/2);k=k(1:end/2)*2;        
        x=x-x(1);%Normalization, a different one (or no normalization) may be used
        plot(k,x)           
    end
    ylim([-40 5])%Modify for improving plot
    xl=xticklabels;
    for l=1:length(xl);xl{l}=num2str(1/str2num(xl{l}));end
    xticklabels(xl);   
    xlabel('mm');
    ylabel('dB');
    grid on  
end
legend('Unfiltered RO','Filtered RO','Unfiltered PE','Filtered PE','Unfiltered SL','Filtered SL')    
title('PSD')
%Note how the slice only goes to 2mm, the readout starts to decrease at 0.8
%of the grid resolution (1.25mm) and the phase encode at 0.6 (1.666mm),
%these is one minus the parameters 0.2 and 0.4. 1.25mm and 1.666mm would
%be the solution to your question for this toy application. In real images
%it may not be so easy to see when the PSDs start to diverge. There you
%could use the PSD ratio, as xLowRes-xHighRes (actually in this case 
%xHighRes is 0, so we are seeing the ratio)