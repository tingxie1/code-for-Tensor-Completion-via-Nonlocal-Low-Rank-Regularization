function [Emsi] = NLRRTC(emsi,Omsi,Omega,mu,alpha)
Emsi=emsi;
sizeData        = size(Emsi);
 par.patsize       =   5;                            % Patch size
 par.patnum        =   10;                           % Initial Non-local Patch number
 par.step          =  2; 
%% main loop
tic
for  iter = 1 : 1
    Curmsi      	= Emsi;
    Curpatch        = Im2Patch3D(Curmsi, par); % image to FBPs    
    Omegapatch      = Im2Patch3D(Omega, par);
    sizePatch       = size(Curpatch);
    % block matching to find samilar FBP goups
    unfoldPatch     = Unfold(Curpatch,sizePatch,3)';
    [label, centroid,~] = fkmeans(unfoldPatch',par.patnum);
    tempPatch       = cell(par.patnum,1);
    omePatch       = cell(par.patnum,1);
     index= cell(par.patnum,1);
    for i = 1:par.patnum
        index{i}=find(label==i);
        tempPatch{i} = Curpatch(:,:,index{i}(:));
        omePatch{i} = Omegapatch(:,:,index{i}(:));
    end    
    
  
    Epatch          = zeros(sizePatch);
    W               = zeros(sizePatch(1),sizePatch(3));
    parfor i = 1:par.patnum  
        tempPatch{i} =LRRTC1(tempPatch{i},omePatch{i},mu,alpha); 
    end
    for i = 1:par.patnum  
        Epatch(:,:,index{i})  = Epatch(:,:,index{i}) + tempPatch{i};
        W(:,index{i})         = W(:,index{i})+ones(size(tempPatch{i},1),size(tempPatch{i},3));
    end 
     EmsiLast = Emsi;
    [Emsi, ~]  =  Patch2Im3D( Epatch, W, par, sizeData); % recconstruct the estimated HSI by aggregating all reconstructed FBP goups.  
    Emsi(logical(Omega)) =Omsi (logical(Omega)); 
end
toc;
end



