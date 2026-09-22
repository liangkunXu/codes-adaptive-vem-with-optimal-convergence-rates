function Vem1fig1()
NNdof=[1025	3911	15538	63346	261247	1080811];
errorestimate=[0.00073221293812672	0.000208843131140529	5.57034730497011e-05	1.37878746735335e-05	3.35131580464649e-06	8.16020351419103e-07];
eigv1=[9.8113736741961	9.68983823321175	9.64889113927064	9.64194029259666	9.64030496596059	9.63986834671297];
ref1=eigv1(end);
eigv2=[15.5679194202884	15.3122530330146	15.2205217578903	15.2024295616874	15.1985608958722	15.1976098150465];
ref2=eigv2(end);
eigv3=[20.3825541268633	19.9129886152585	19.7803470902706	19.7482005140098	19.7414476779286	19.7398148625007];
ref3=eigv3(end);
eigv4=[30.7533183620552	29.8849562078775	29.5948126744794	29.5381850001157	29.5257646991643	29.5226127675472];
ref4=eigv4(end);
eigv5=[33.149823495295	32.2542913408914	32.0020401252792	31.9316327297098	31.9177688498429	31.9138735969868];
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
x=linspace(3000,100000);y=10.^(au(1)-1.8)*x.^(-1);hold on; loglog(x,y,'k-','LineWidth',1)
%添加主题和坐标标题
xlabel('$\#\mathcal{T}_\ell$','Interpreter','latex','fontsize',11);ylabel('Error','fontsize',12);

%添加图例
fg=legend('$|\lambda_{1}-\lambda_{1,\ell}|$','$|\lambda_{2}-\lambda_{2,\ell}|$', ...
    '$|\lambda_{3}-\lambda_{3,\ell}|$','$|\lambda_{4}-\lambda_{4,\ell}|$', ...
    '$|\lambda_{5}-\lambda_{5,\ell}|$','$\eta_\ell^2$','location','best');
fg.Interpreter="latex";
tx=text(13000,10^-3.2,'slope$=-1$','fontsize',12);
tx.Interpreter="latex";
axis tight;
hold off
end
