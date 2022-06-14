
function [Y, errList]  =  LRRTC(D,ome,mu,alpha)
sizeD          = size(D);
ndim           = length(sizeD);
Rank           = ones(1,ndim)*6;
Rank(1)        = sizeD(1)-20;
maxIter        = 400   ;  
errList = zeros(maxIter, 1);
%% initialization about M
M         = cell(ndim, 1);
Lam       = cell(ndim, 1);
tempX     = cell(ndim, 1);
sValue    = cell(ndim, 1);
Mnorm     = zeros(ndim, 1);
Msum      = zeros(sizeD);
for i = 1:ndim
    M{i}      = D;
    Lam{i}    = zeros(sizeD);
    tempX{i}  = Unfold(D, sizeD, i);
    sValue{i} = svd(tempX{i}, 'econ');
    Mnorm(i)  = min(sum(sValue{i}>5),Rank(i));
    Msum      = Msum + Fold(M{i},sizeD,i);
end
%% initialization about other parameters
LamSum    = zeros(sizeD);
temp_n    = zeros(1,ndim);
Y         = D;
tic;
%% main loop
for i = 1:300    
   if mod(i, 10) == 0
        fprintf('lrrtc: iterations = %d   difference=%f\n', i, errList(i-1));
    end
    %% updating M
    Msum   = 0*Msum;
    LamSum = 0*LamSum;
     mu    = mu*1.1;
    for k = 1:ndim
   [tempX{k}, temp_n(k), sValue{k}, Mnorm(k)] = Pro2WNNM(Unfold(double(Y + Lam{k}/mu), sizeD, k),alpha(k)/mu);%,Rank(k));
       M{k}      = Fold(tempX{k}, sizeD, k);     
        Msum      = Msum + M{k};
        LamSum    = LamSum + Lam{k}; % updating the multipliers       
    end   
    lastY = Y;
     Y = (mu* Msum- LamSum) / (ndim*mu); 
     Y(logical(ome)) =D(logical(ome)); 
     for k = 1:ndim
        Lam{k}    = Lam{k}-mu*(M{k}-Y);
     end
   errList(i) = norm(Y(:)-lastY(:)) / norm(Y(:));
   if errList(i) < 1e-4
  break;
    end
end
toc;
end



