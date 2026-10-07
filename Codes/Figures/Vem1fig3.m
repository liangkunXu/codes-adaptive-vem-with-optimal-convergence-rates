function Vem1fig3()
% 在SL型域域山求解前5个特征值
NNdof=[643	2885	14225	55276	220509	894696];
errorestimate=[0.000238043427925615	5.6466211808914e-05	1.27997616201333e-05	3.48631321705253e-06	8.92991065952726e-07	2.19929417070763e-07];

eigv1=[35.4526679917161	34.03626933268	33.5690605775501	33.5064855140301	33.4909313757645	33.4866342478284];
ref1=eigv1(end);

eigv2=[51.2634636049152	50.1703069597631	49.4680618274363	49.3825261541545	49.3572374364383	49.3502590526081];
ref2=eigv2(end);

eigv3=[69.491332138883	68.0282653457098	66.795263430072	66.6432269338101	66.5976837994584	66.5853488616109];
ref3=eigv3(end);

eigv4=[82.9280398731542	80.8535170931403	79.2582475230904	79.0474906837564	78.9804233107537	78.9625140327813];
ref4=eigv4(end);

eigv5=[118.503018600813	115.448136085416	112.479501773445	112.048473060999	111.948216973776	111.919782078715];
ref5=eigv5(end);

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
x=linspace(8000,200000);y=10.^(au(1)-2)*x.^(-1);hold on; loglog(x,y,'k-','LineWidth',1)
%添加主题和坐标标题
xlabel('$\#\mathcal{T}_\ell$','Interpreter','latex','fontsize',11);ylabel('Error','fontsize',12);
%添加图例
fg=legend('$|\lambda_{1}-\lambda_{1,\ell}|$','$|\lambda_{2}-\lambda_{2,\ell}|$', ...
    '$|\lambda_{3}-\lambda_{3,\ell}|$','$|\lambda_{4}-\lambda_{4,\ell}|$', ...
    '$|\lambda_{5}-\lambda_{5,\ell}|$','$\eta_\ell^2$','location','best');
fg.Interpreter="latex";
tx=text(30000,10^-3.9,'slope$=-1$','fontsize',12);
tx.Interpreter="latex";
axis tight;
hold off
end
