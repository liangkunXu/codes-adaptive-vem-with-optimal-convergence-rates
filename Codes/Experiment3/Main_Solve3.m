function [etaN,dofN,eigenN,NT]=Main_Solve3()
%% Adaptive Virtual Element Method
%% Parameters
theta = 0.5;  Tol = 7e-03; epsrong=1; %Tol = 3e-04
%% S-shaped region with cracks
node = [0,1;0.5,1;1,1;0,0.5;0.5,0.5;1,0.5;0,0;0.5,0;1,0;1,0.5];
elem=[1 4 2; 5 2 4; 2 5 3; 6 3 5; 4 7 5; 8 5 7;5 8 10;9 10 8];
% Uniform encryption on a coarse grid
for i=1:3
    [node,elem] = uniformrefine(node,elem);
end
showmesh(node,elem);
pause(0.025);
etaN = [];
dofN=etaN;
NT=etaN;
eigenN=etaN;
k=1;
% solve
[Eigen,~,uh,info,item,~] = Directlysolveeigenvalues3(node,elem);
while epsrong > 0.5*Tol
    % Step 1: SOLVE
    flagg=size(elem,1); 
    jishi=1;
    while  true
        if jishi==1
            % Step 4: ESTIMATE
            eta=0;
            for arg=1:item
                eta=eta+indicatorBVP3(node,elem,uh(:,arg),1/Eigen(arg)*uh(:,arg),info);
            end
            eta1=eta./item;
            etaN(k) = sum(eta1);
        else
            % Step 4: Solve
            [uhb,uh,~,info]=SolveBVP3(node,elem,uh,info,item,HB,belong);
            % Step 4: ESTIMATE
            eta=0;
            for arg=1:item
                eta = eta+indicatorBVP3(node,elem,uh(:,arg),uhb(:,arg),info);
            end
            eta1=eta./item;
        end

        if sum(eta1)^0.5 < epsrong
            if flagg==size(elem,1)
                epsrong=0.5*epsrong;
                continue;
            else
                % solve
                [Eigen,dof,uh,info,item,kNT] = Directlysolveeigenvalues3(node,elem);
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
            [node,elem,~,~,HB,belong]=bisect(node,elem,elemMarked,bdFlag);
            showmesh(node,elem);
            pause(0.025);
            jishi=jishi+1;
        end
    end
    epsrong=0.5*epsrong;
    k=k+1;
end
%uhh=uh(:,2);
%showsolution(node,elem,uhh(1:size(node,1),1));
end