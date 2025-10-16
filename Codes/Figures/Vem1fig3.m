function Vem1fig3()
% 在SL型域域山求解前5个特征值
%NNdof=[462	1733	7887	33722	140957	433906	1816514];
NNdof=[1002	3589	16026	67958	282970	869702	3636932];
errorestimate=[0.00017474502732705	5.30606036423727e-05	1.19370463444403e-05	2.79955462577445e-06	6.67810719946677e-07	2.36948494102752e-07	5.52854650764236e-08];

eigv1=[34.4311286191117	33.7488640638226	33.5410364257742	33.4969601188562	33.4877864599987	33.4860890098625	33.4854886304384];
ref1=eigv1(end);

eigv2=[49.9109432877075	49.4743820864767	49.378500844312	49.3547125706669	49.3495235416497	49.3485205133464	49.3481442837588];
ref2=eigv2(end);

eigv3=[67.6794044436472	66.794198917172	66.6390401907047	66.5939100846878	66.583991043978	66.5819304224197	66.5813512482329];
ref3=eigv3(end);

eigv4=[80.6951717315469	79.3374479721106	79.0512772174647	78.9778945231684	78.9615873934615	78.9581780104158	78.9571816814588];
ref4=eigv4(end);

eigv5=[115.332585018903	112.591778729354	112.099239144525	111.952910585423	111.919933943172	111.912676173525	111.910728402761];
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