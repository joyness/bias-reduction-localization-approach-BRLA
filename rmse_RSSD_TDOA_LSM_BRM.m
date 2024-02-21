%copyright@ Yuanyuan Zhang,School of Computer and Information Engineering,Changzhou Institute of Technology.
%This code is the reference of the paper:"An Efficient Estimator for Source Localization in WSNs
%using RSSD and TDOA Measurements"
% ---------------------------------------------------------
clear;clc;close all
%% Initialize
N_obj=1000;          % Number of iterations
N_MC=1000;           % Number of iterations
noi=9;
gamma =4.5;
L_noise=length(noi);
%% RSSD_TDOA_LSM_BRM
mse_RSSD_TDOA_LSM_BRM_x = zeros(L_noise,1);
mse_RSSD_TDOA_LSM_BRM_t = zeros(L_noise,1);
ACT_WIV_CDF=zeros(N_MC,1);
warning off
for nn2=1:L_noise 
    disp(['parameter #' num2str(nn2) '/' num2str(numel(noi))  ]);
    for mc=1:N_MC % loop for Monte Carlo trials
        disp(['... Monte Carlo simulation #' num2str(mc) '/' num2str(N_MC)]);
        x =[8;26];%
        ai1 =[10;0;0;5;20;15;25;5;30];%
        ai2 =[10;0;5;0;15;20;5;25;20];%
        ai=[ai1,ai2];
        ar1=[ai1(1),ai2(1)];
        rdi1=x(1)-ai1(1);
        rdi2=x(2)-ai2(1);
        x_r1=[rdi1,rdi2];
        r1 = sqrt(sum((x_r1).^2));
        aa = [ai1,ai2];        % sort the coordinates into a matrix.
        aa = aa';              % size(aa)==(2,N)
        aa(:,1)=[];
        m = length(aa(1,:));       % Number of sensors after deleting¡¡refer¡¡anchor.
        r = sqrt(sum((aa-repmat(x,1,m)).^2));
        r = r';                % size(r)==(m,1)
        noise1=sqrt(noi(nn2))*randn(m,1);
        %%RSSD
        P1r=10*gamma*log10(r/r1)+noise1;%the noise RSSD measurement with respect to the first sensor
        C1r =  10.^(P1r/(5*gamma));            % the corresponding C1r = [C1r_1,...,C1r_m]^T.
        V1r=C1r* log(10)/(5*gamma);
        A0=[2*(C1r*ar1-aa'),1-C1r];
        b0= C1r*sum(ar1.^2)-sum(aa.^2)';
        W1=diag(1- C1r/sum( C1r));
        W0=W1^(-1);
        W=[W0, zeros(size(W0,1),size(W0,1));
            zeros(size(W0,1),size(W0,1)),W0];
        %%TDOA
        noise2=sqrt(noi(nn2))*randn(m,1);
        RD=r-ones(m,1).*r1+noise2;
        G1=-[aa'-ones(m,1).*ar1,RD];
        h1=1/2*(RD.^2-sum(aa.^2)'+ones(m,1)*sum(ar1.^2));
        Q=diag(ones(m,1)*(noi(nn2)));
        %%X=[x' r1 ||x||^2]
        A0v=[2*(C1r*ar1-aa'),zeros(m,1),1-C1r];%%RSSD
        G1v=[G1, zeros(m,1)];%%TDOA
        mixA2=[A0v;G1v];
        mixb2=[b0;h1];
        tilde_x_mixLS2 =((mixA2'*mixA2)^-1)*mixA2'*mixb2;
        
        tic
        method = 'nnlsm_activeset';
        Result = solveNNLS(mixA2,mixb2,method);
        x_t =Result(1:2);
        %%BRM
        Aw=W*mixA2;
        rw=sqrt(sum((aa-repmat(x_t,1,m)).^2));
        rw=rw';
        r1_e=sqrt(sum(( x_t(1)-ai1(1)).^2)+(x_t(2)-ai2(1)).^2);
        P1r_e=10*gamma*log10(rw/r1_e)+noise1;
        C1r_e = 10.^(P1r_e/(5*gamma));
        Grssd_A0=[2*(C1r_e*ar1-aa'),zeros(m,1),1-C1r_e];
        %%
        RD_e=rw-ones(m,1).*r1_e+noise2;
        Gtdoa_G1=-[aa'-ones(m,1).*ar1,RD_e,zeros(m,1)];%%
        G=[Grssd_A0; Gtdoa_G1];
        Gw=W*G;
        Cin_spm=Aw'* Gw *(Gw'* Gw)^-1* Gw'* Aw;
        Fg=[eye(2);2*(x_t(1)-ai1(1))*((r1_e+Result(3)).^(-1)),2*(x_t(2)-ai2(1))*((r1_e+Result(3)).^(-1));2*x_t'];
        Ptau=[0; 0;r1_e-Result(3);sum(x_t.^2)- Result(4)];
        x_t2=x_t-(Fg'*Cin_spm*Fg)^-1*Fg'*Cin_spm*Ptau;
        mse_RSSD_TDOA_LSM_BRM_t(nn2) = mse_RSSD_TDOA_LSM_BRM_t(nn2)+toc;
        % Calculate mse of x
        mse_RSSD_TDOA_LSM_BRM_x(nn2)=mse_RSSD_TDOA_LSM_BRM_x(nn2)+sum((x_t2-x).^2);
        ACT_WIV_CDF(mc)=sqrt(sum(x_t2-x).^2);
    end
end
