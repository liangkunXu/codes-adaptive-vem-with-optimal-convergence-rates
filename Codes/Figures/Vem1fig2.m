function Vem1fig2()
% 在SS型域域山求解前6个特征值
%NNdof=[551	2134	9893	38635	226221	988882];%DOF
NNdof=[1141	4345	19962	77621	453379	1979720];%NT
errorestimate=[0.000797433187417691	0.000212478700575169	5.37799435207698e-05	1.44707931566704e-05	3.37243158608173e-06	8.56231889701704e-07];

eigv1=[19.6383097714374	19.5967885095133	19.5699969817813	19.5641689888785	19.5625712363173	19.562527232517];
ref1=eigv1(end);

eigv2=[19.7734139760737	19.7140902924483	19.6863038774025	19.6801182143939	19.6784780083649	19.6784325742696];
ref2=eigv2(end);

eigv3=[50.2213241894106	49.2882552038109	48.9149311484497	48.8474937455224	48.8310118069735	48.8286388135706];
ref3=eigv3(end);

eigv4=[50.2213241894106	49.4744053830165	49.1798814885896	49.1083423626248	49.0862744538178	49.0809050370341];
ref4=eigv4(end);

eigv5=[50.7674223592945	49.6716150305379	49.2983224601317	49.2199488572367	49.200260487227	49.1970538042334];
ref5=eigv5(end);

eigv6=[50.7674223592945	49.6716150305379	49.2983224601317	49.24004060754	49.2235861911	49.2212087675303];
ref6=eigv6(end);

%运用loglog作图
%% 作误差曲线
loglog(NNdof(1:4),abs(eigv1(1:4)-ref1),'bs-.','MarkerSize',8,'LineWidth',2);
hold on 
loglog(NNdof(1:4),abs(eigv2(1:4)-ref2),'md-.','MarkerSize',8,'LineWidth',2);
loglog(NNdof(1:4),abs(eigv3(1:4)-ref3),'ro-.','MarkerSize',8,'LineWidth',2);
loglog(NNdof(1:4),abs(eigv4(1:4)-ref4),'c>-.','MarkerSize',8,'LineWidth',2);
loglog(NNdof(1:4),abs(eigv5(1:4)-ref5),'gp-.','MarkerSize',8,'LineWidth',2);
loglog(NNdof(1:4),abs(eigv6(1:4)-ref6),'ko-.','MarkerSize',8,'LineWidth',2);
%% 作误差指示子曲线
loglog(NNdof(1:4),errorestimate(1:4),'kX-','MarkerSize',8,'LineWidth',2);
%求回归系数
au=regress(log10(abs(eigv1(1:end-1)-ref1))',[ones(length(eigv1)-1,1) log10(NNdof(1:end-1))'])
%a=au(1),b=au(2),10
x=linspace(4000,30000);y=10.^(au(1)-2)*x.^(-1);hold on; loglog(x,y,'k-','LineWidth',1)
%添加主题和坐标标题
xlabel('$\#\mathcal{T}_\ell$','Interpreter','latex','fontsize',11);ylabel('Error','fontsize',12);

%添加图例
fg=legend('$|\lambda_{1}-\lambda_{1,\ell}|$','$|\lambda_{2}-\lambda_{2,\ell}|$', ...
    '$|\lambda_{3}-\lambda_{3,\ell}|$','$|\lambda_{4}-\lambda_{4,\ell}|$', ...
    '$|\lambda_{5}-\lambda_{5,\ell}|$','$|\lambda_{6}-\lambda_{6,\ell}|$', ...
    '$\eta_\ell^2$','location','best');
fg.Interpreter="latex";
tx=text(9000,10^-2.98,'slope$=-1$','fontsize',12);
%tx=text(26000,10^-3.4,'slope$=-1$','fontsize',12);
tx.Interpreter="latex";
axis tight;
hold off
end