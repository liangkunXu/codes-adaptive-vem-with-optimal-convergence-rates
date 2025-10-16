function Vem1fig1()
% 在L型域域山求解前5个特征值
%NNdof=[1093	3964	18416	62543	275290	1188054];
NNdof=[2258	8070	37170	125650	551778	2378878];%NT
errorestimate=[0.000757662717895104	0.000214676730177507	4.92237645358673e-05	1.40487966804863e-05	3.21126458255986e-06	7.69164453137186e-07];

eigv1=[9.68027955458956	9.6488482995807	9.64201987167898	9.64024331748183	9.63985813235456	9.63975593017892];
ref1=eigv1(end);

eigv2=[15.2347879244154	15.2021441788682	15.1995141024825	15.1976213637404	15.1973692310332	15.1972839925513];
ref2=eigv2(end);

eigv3=[19.8233828481596	19.7545894422882	19.7441107761932	19.7402369999171	19.7394698594893	19.7392790247862];
ref3=eigv3(end);

eigv4=[29.7339865156865	29.5702970727998	29.5341190065409	29.524605747607	29.5221912464946	29.5216579691475];
ref4=eigv4(end);

eigv5=[32.2645180987301	31.9962100204309	31.9336852228776	31.9179320300791	31.913937452121	31.9129439980999];
ref5=eigv5(end);

%运用loglog作图
%% 作误差曲线
loglog(NNdof,abs(eigv1-ref1),'bs-.','MarkerSize',8,'LineWidth',2);
hold on 
loglog(NNdof,abs(eigv2-ref2),'md-.','MarkerSize',8,'LineWidth',2);
loglog(NNdof,abs(eigv3-ref3),'ro-.','MarkerSize',8,'LineWidth',2);
loglog(NNdof,abs(eigv4-ref4),'c>-.','MarkerSize',8,'LineWidth',2);
loglog(NNdof,abs(eigv5-ref5),'gp-.','MarkerSize',8,'LineWidth',2);
%% 作误差指示子曲线
loglog(NNdof(1:end-1),errorestimate(1:end-1),'kX-','MarkerSize',8,'LineWidth',2);

%求回归系数
au=regress(log10(abs(eigv1(1:end-1)-ref1))',[ones(length(eigv1)-1,1) log10(NNdof(1:end-1))'])
%a=au(1),b=au(2),10
x=linspace(19000,200000);y=10.^(au(1)-1.3)*x.^(-1);hold on; loglog(x,y,'k-','LineWidth',1)
%添加主题和坐标标题
xlabel('$\#\mathcal{T}_\ell$','Interpreter','latex','fontsize',11);ylabel('Error','fontsize',12);

%添加图例
fg=legend('$|\lambda_{1}-\lambda_{1,\ell}|$','$|\lambda_{2}-\lambda_{2,\ell}|$', ...
    '$|\lambda_{3}-\lambda_{3,\ell}|$','$|\lambda_{4}-\lambda_{4,\ell}|$', ...
    '$|\lambda_{5}-\lambda_{5,\ell}|$','$\eta_\ell^2$','location','best');
fg.Interpreter="latex";
tx=text(8000,10^-3.3,'slope$=-1$','fontsize',12);
tx.Interpreter="latex";
axis tight;
hold off
end