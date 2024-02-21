clear;clc;close all
%% Initialize
N_obj=1000;          % Number of iterations
N_MC=1000;           % Number of iterations
noi=9;
gamma =4.5;
L_noise=length(noi);
%% RSSD_TDOA_WLS
mse_RSSD_TDOA_WLS_x = zeros(L_noise,1);
mse_RSSD_TDOA_WLS_t = zeros(L_noise,1);
WLS_CDF=zeros(N_MC,1);
warning off
for nn2=1:L_noise
    disp(['parameter #' num2str(nn2) '/' num2str(numel(noi))  ]);
    for mc=1:N_MC % loop for Monte Carlo trials
        disp(['... Monte Carlo simulation #' num2str(mc) '/' num2str(N_MC)]);
        x =[8;26];
        ai1 =[10;0;0;5;20;15;25;5;30];
        ai2 =[10;0;5;0;15;20;5;25;20];
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
        r = r';                % size(r)==(m,1), r=[r_1,r_2,...,r_m]^T
        
        noise1=sqrt(noi(nn2))*randn(m,1);
        %%RSSD
        P1r=10*gamma*log10(r/r1)+noise1;%the noise RSSD measurement with respect to the first sensor
        C1r =  10.^(P1r/(5*gamma));            % the corresponding C1r = [C1r_1,...,C1r_m]^T.
        V1r=C1r* log(10)/(5*gamma);
        A0=[2*(C1r*ar1-aa'),1-C1r];
        b0= C1r*sum(ar1.^2)-sum(aa.^2)';
        W1=diag(1- C1r/sum( C1r));
        W0=sqrt(W1);
        W0=W0^(-1);%
        %%TDOA
        noise2=sqrt(noi(nn2))*randn(m,1);
        RD=r-ones(m,1).*r1+noise2;%the noise TDOA measurement with respect to the first sensor
        G1=-[aa'-ones(m,1).*ar1,RD]; %
        h1=1/2*(RD.^2-sum(aa.^2)'+ones(m,1)*sum(ar1.^2)); %
        % WLS
        %%X=[x' r1 ||x||^2]
        A0v=[2*(C1r*ar1-aa'),zeros(m,1),1-C1r];%%RSSD
        G1v=[G1, zeros(m,1)];%%TDOA
        mixA2=[A0v;G1v];
        mixb2=[b0;h1];
        Wt =diag(noise2.^(-2)) ;
        W=[Wt, zeros(size(Wt,1),size(Wt,1));
            zeros(size(Wt,1),size(Wt,1)),Wt];
        tic
        tilde_x_mixLS2 =((mixA2'*W*mixA2)^-1)*mixA2'*W*mixb2;
        x_t= tilde_x_mixLS2(1:2);
        
        mse_RSSD_TDOA_WLS_t(nn2) = mse_RSSD_TDOA_WLS_t(nn2)+toc;
        % Calculate mse of x
        mse_RSSD_TDOA_WLS_x(nn2)=mse_RSSD_TDOA_WLS_x(nn2)+sum((x_t-x).^2);%%  2.59641 4.61268 5.34078 6.15298 6.52514
        WLS_CDF(mc)=sqrt(sum(x_t-x).^2);
    end
end