%% Used for plotting the CDF of various localization algorithms from their .mat files
 clear;clc;close all
load 'ACT_CDF.mat'%
load 'LSM_BRM_CDF.mat'%
load 'SRWLS_CDF.mat'
load 'LS_CDF.mat'
load 'WLS_CDF.mat'
load 'IPM_CDF.mat'
load 'RSSD_ACT_CDF.mat'
load 'TDOA_ACT_CDF.mat' 

%% Plot RMSE  
figure(1);box on;
% yline([0.9],'--k','90%','LineWidth',1.5,'LabelHorizontalAlignment','center');hold on;
warning off
h1=cdfplot(LS_CDF);
set(h1,'linestyle','-','Color','b','LineWidth',2.0,'DisplayName','LLS');hold on;
h2=cdfplot(WLS_CDF);
set(h2,'linestyle','-','Color','[0.4 0.6 1]','LineWidth',2.0,'DisplayName','WLS');
h3=cdfplot(SRWLS_CDF);
set(h3,'linestyle','-','Color','[0.7 0 1]','LineWidth',2.0,'DisplayName','SRWLS');
h4=cdfplot(IPM_CDF );
set(h4,'linestyle','-','Color','g','LineWidth',2.0,'DisplayName','IPM');
h5=cdfplot(ACT_CDF);
set(h5,'linestyle','- ','Color','m','LineWidth',2.0,'DisplayName','LSM');
h6=cdfplot(ACT_WIV_CDF);
set(h6,'linestyle','-','Color','r','LineWidth',2.0,'DisplayName','BRLA');
h7=cdfplot(RSSD_ACT_CDF);
set(h7,'linestyle','-','Color','c','LineWidth',2.0,'DisplayName','BRLA$_{RSSD}$');
h8=cdfplot(TDOA_ACT_CDF );
set(h8,'linestyle','-','Color','[1 0.5 0]','LineWidth',2.0,'DisplayName',' BRLA$_{TDOA}$ ');

legend_handle  = legend([h1 h2 h3 h4 h5 h6 h7 h8],...
{'LLS','WLS','SRWLS','IPM','LSM','BRLA','BRLA$_{RSSD}$','BRLA$_{TDOA}$'},'Interpreter','Latex');
 legend('show','Interpreter','latex','FontSize',10,'Location','Best');
yhandle = ylabel('CDF','Interpreter','Latex');
xhandle = xlabel('$||\widehat{\theta}-\theta||$ (m) ','Interpreter','latex');
title('')
set(gca,'xlim',[0,25]);
set(gca,'xtick',0:25/5:25);
grid on

