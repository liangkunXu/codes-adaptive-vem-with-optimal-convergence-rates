function [etaN,dofN,NT,eigenN]=Main_Solver()
%% Adaptive Virtual Element Method
tic;
%% Parameters
theta = 0.5;  Tol = 1e-02; epsrong=1; % TOL=1e-03
%% Generate an initial mesh
data=load('meshs.mat');
node=data.meshs.node1;
elemt=data.meshs.elem1;
elem=elemt(:,1:3);
log13=false(size(elem,1),1);
log13(elemt(:,4)==1,:)=true;% 1 3
showmesh(node,elem);
pause(0.025);
etaN = [];
dofN=etaN;
NT=etaN;
eigenN=etaN;
k=1;
% solve
[Eigen,~,uh,info,item,~]=Directlysolveeigenvalues2(node,elem,log13);
while epsrong > 0.5*Tol
    % Step 1: SOLVE
    flagg=size(elem,1); 
    jishi=1;
    while  true
        if jishi==1
            % Step 4: ESTIMATE
            eta=0;
            for arg=1:item
                eta=eta+indicatorBVP2(node,elem,uh(:,arg),1/Eigen(arg)*uh(:,arg),info,log13);
            end
            eta1=eta./item;
        else
            % Step 4: Solve
            [uhb,uh,~,info,log13]=SolveBVP2(node,elem,uh,info,item,HB,belong,log13);
            % Step 4: ESTIMATE
            eta=0;
            for arg=1:item
                eta = eta+indicatorBVP2(node,elem,uh(:,arg),uhb(:,arg),info,log13);
            end
            eta1=eta./item;
        end

        if sum(eta1)^0.5 < epsrong
            if flagg==size(elem,1)
                epsrong=0.5*epsrong;
                continue;
            else
                % solve
                [Eigen,dof,uh,info,item,kNT]=Directlysolveeigenvalues2(node,elem,log13);
                dofN(k)=dof;
                NT(k)=kNT;
                eigenN(k,:)=Eigen';
                etaN(k) = sum(eta1);
                break;
            end
        else
            % Step 3: MARK
            elemMarked = mark(elem,eta1.^0.5,theta);
            % Step 4: REFINE
            bdFlag = SetBoundary(node,elem,'Dirichlet');
            [node,elem,~,~,HB,belong] = bisect(node,elem,elemMarked,bdFlag);
            belong(belong==0)=find(belong==0);
            log13=log13(belong);
            showmesh(node,elem);
            pause(0.025);
            jishi=jishi+1;
        end
    end
    epsrong=0.5*epsrong;
    k=k+1;
end
toc
end