function [Xi,de,di,vic,viso,dest] = fittingFuntion_final(DTI,REP,bval,slice,mask)
    
    %% PARAMETERS
    bvalU = unique(bval); bvalU = bvalU';
    N    = 16;
    K    = 15;
    [XGO,YGO,ZGO,~] = size(DTI);
    BVAL  = length(bvalU);
    
    mask_final = mask;
    
    %% FIT
    x = bvalU';       
    
    Xi = zeros(XGO,YGO,ZGO);
    de = zeros(XGO,YGO,ZGO);    
    di = zeros(XGO,YGO,ZGO);
    
    vic_out = zeros(XGO,YGO,ZGO);
    viso_out = zeros(XGO,YGO,ZGO);    
    dest_out = zeros(XGO,YGO,ZGO);
    
    % preparing permutations
    for Z = slice
        
       vic  = zeros(XGO,YGO,1,REP);
       viso = zeros(XGO,YGO,1,REP);
       dest = zeros(XGO,YGO,1,REP);
       v0   = zeros(XGO,YGO,1,REP);

        Sb = zeros(XGO,YGO,1,BVAL,REP); % signal averaged across bdir at each bval  (eq 5)
        for rp = 1 : REP
            for bi = 1 : BVAL
                temp = DTI(:,:,Z,bval == bvalU(bi));
                if bi > 1
                    subset   = randperm(N);
                    subset   = subset(1:K);
                    %subset = [(1:REP-1) REP+1:N];
                    temp   = temp(:,:,1,subset);
                    Sb(:,:,1,bi,rp) = mean(temp,4);
                else
                    Sb(:,:,1,bi,rp) = mean(temp,4);
                end

            end
        end

        parfor R = 1 : REP      

            disp(['fitting slice = ' num2str(Z) ' repetition = ' num2str(R)])
            
            Sbm = Sb(:,:,1,:,R);

                for X = 1 : XGO
                    for Y = 1 : YGO      

                        if mask(X,Y,Z) == 1

                            yx     = squeeze(Sbm(X,Y,1,1:end,1)); 
                            yx     = yx./(yx(1));                   

                            y = @(b,x) (1-b(3)).*(b(1).*exp(-b(1).*(1.7*1E-03.*x)) + (1-b(1)).*exp(-(1-b(1)).*b(2).*x))+b(3).*exp(-3*1E-03.*x);
                            OLS = @(b) sum(sum((y(b,x)-yx).^2)); % Ordinary Least Squares cost function
                            opts = optimset('MaxFunEvals',50000, 'MaxIter',10000,'Display','off'); % options
                            B = fminsearchbnd(OLS, [0 0 0.5], [0 0 0], [1 Inf 1], opts);

                            vic(X,Y,1,R)  = B(1);
                            dest(X,Y,1,R) = B(2);
                            viso(X,Y,1,R) = B(3);

                        end

                    end
                end
        end
    
       % MIP
       [viso, ~] = max(viso,[],4);    
       [vic, ~]  = max(vic,[],4);  
       [dest, ~] = max(dest,[],4);
       
       vic_out(:,:,Z)  = vic;
       viso_out(:,:,Z) = viso;
       dest_out(:,:,Z) = dest;
       
    end 
    
    viso = viso_out;
    vic = vic_out;
    dest = dest_out;
    
    %% TAKING OUT DARK VOXELS AND INTERPOLATING WITH DELAUNAY 
    cutOFF = 0.15;
    brain = viso+vic+dest; brain(brain>0)=1; 
    se = strel('disk',10);
    brain = imclose(brain,se);
    
    mask = viso<=cutOFF;
    viso(mask == 1) = 0;
    
    viso_old = viso;
    vic_old = vic;
    dest_old = dest;
    
    for Z = slice
    
	    %%        
	    cutOFF = 0.15;
	    slice_toInterp = viso(:,:,Z);

	    if(sum(slice_toInterp(:))~=0)
	    
		    [x,y] = size(slice_toInterp);
		    [X,Y] = meshgrid(1:x,1:y);

		    lattice_tointerp = (slice_toInterp == 0);
		    lattice_tointerp = lattice_tointerp.*double(brain(:,:,Z));
		    lattice_signal   = not(lattice_tointerp);
		    lattice_signal   = lattice_signal.*double(brain(:,:,Z));

		    Xhave = X(lattice_signal == 1);
		    Yhave = Y(lattice_signal == 1);
		    Shave = slice_toInterp(lattice_signal == 1);

		    sliceInterp = griddata(Xhave,Yhave,Shave,X,Y);
		    sliceInterp(lattice_signal == 1) = slice_toInterp(lattice_signal == 1);
		    viso(:,:,Z) = sliceInterp;

            end
    
	    %%
	    cutOFF = 0.15;
	    csf = viso(:,:,Z)>0.9;
	    mask = vic(:,:,Z)<=cutOFF;
	    mask = mask.*not(csf);
	    slice_toInterp = vic(:,:,Z);
	    slice_toInterp(mask == 1) = 0;
	    
	    if(sum(slice_toInterp(:))~=0)
	    
		    [x,y] = size(slice_toInterp);
		    [X,Y] = meshgrid(1:x,1:y);

		    lattice_tointerp = (slice_toInterp == 0);
		    lattice_tointerp = lattice_tointerp.*double(brain(:,:,Z));
		    lattice_signal   = not(lattice_tointerp);
		    lattice_signal   = lattice_signal.*double(brain(:,:,Z));

		    Xhave = X(lattice_signal == 1);
		    Yhave = Y(lattice_signal == 1);
		    Shave = slice_toInterp(lattice_signal == 1);

		    sliceInterp = griddata(Xhave,Yhave,Shave,X,Y);
		    sliceInterp(lattice_signal == 1) = slice_toInterp(lattice_signal == 1);
		    vic(:,:,Z) = sliceInterp;

   	    end
    
    end
    
    viso = viso.*mask_final;
    vic  = vic.*mask_final;
    dest = dest.*mask_final;
    
    Xi  = (1-viso).*(1-vic)+viso;
    de  = ((1-viso).*((1-vic).^2).*dest)./((1-viso).*(1-vic)+viso)+(viso.*(3*1E-03))./((1-viso).*(1-vic)+viso);
    di  = vic.*1.7*1E-03; 
    
    Xi = Xi.*mask_final;
    de = de.*mask_final;
    di = di.*mask_final;
      
end