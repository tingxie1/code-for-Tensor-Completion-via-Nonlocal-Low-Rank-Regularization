
% This is an implementation of the algorithm for tensor completion
% 
% Please cite the following paper if you use this code:
%T. Xie, S. Li, L. Fang and L. Liu
%''Tensor Completion via Nonlocal Low-Rank Regularization,''
%% in IEEE Transactions on Cybernetics, vol. 49, no. 6, pp. 2344-2354, June 2019.
%Email: xting@hnu.edu.cn
%-------------------------------------------------------------------

 clc;
 clear;
close all;
warning off
addpath('tensor_toolbox');
load j.mat   
T=z;
[n1,n2,n3]             =        size(T)                      ;
p                      =        0.4                       ;
Omega                  =        zeros(size(T))               ;
chosen                 =        randperm(n1*n2*n3,...
                                       round(p*n1*n2*n3))     ;
Omega(chosen)          =        1                             ;

A                      =        diag(sparse(double(Omega(:)))); % sampling operator
b                      =        A * T(:)                     ; % available data
bb                     =        reshape(b,[n1,n2,n3]);


%% the parameter setting for Natural HSIs
mu1=  1/120        ;     
mu2=  1/160        ;       
alpha1     = [1 1 25];   
alpha2     = [1 1.5 1.2]; 

%% the parameter setting for  Remote Sensing HSIs
% mu1=  1/100        ;     
% mu2=  1/160        ;       
% alpha1     = [1 1 2.25];   
% alpha2     = [1 6 5]; 

%% LRR-TC
 [Y, errList]  =  LRRTC(bb,Omega,mu1,alpha1);

%%NLRR-TC
 patsize       =   5;                            % Patch size
 step          =  2; 
 int_x        =   (n1-patsize)/step ;
 int_y        =   (n2-patsize)/step ;

Y1=[Y;Y(n1,:,:)];
ad=[Y(:,n2,:);Y(n1,n2,:)];
Yn=[Y1 ad];
T1=[T;T(n1,:,:)];
ad1=[T(:,n2,:);T(n1,n2,:)];
Tn=[T1 ad1];
Omegaa_n=[Omega;Omega(n1,:,:)];
ad2=[Omega(:,n2,:);Omega(n1,n2,:)];
Omega_n=[Omegaa_n ad2];

 if  (int_x     ==   fix(int_x))  &&  (int_y     ==   fix(int_y) )
        [YM]  =  NLRRTC(Yn(1:n1,1:n2,:),Tn(1:n1,1:n2,:),Omega_n(1:n1,1:n2,:),mu2,alpha2);
 else
        [YM]  =  NLRRTC(Yn,Tn,Omega_n,mu2,alpha2);
 end
 Yt=YM(1:n1,1:n2,:);
[psnr, ssim, fsim, ergas, msam] = MSIQA(255*T, 255*Yt)