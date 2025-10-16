function [Eigen,dof,uh,info,item,NT] = Directlysolveeigenvalues3(node,elem)
%case 1
Axi=1;c=0;sagamaa=1;
%% Get auxiliary data
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid;  diameter = aux.diameter;  
N = size(node,1);  NT = size(elem,1);
Nm = 3;
%% Compute and assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
Ph0 = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);  Dm = cell(NT,1);
elemLen = cellfun('length',elem); 
nnz = sum(elemLen.^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); 
ssA = zeros(nnz,1);  ssB = zeros(nnz,1);
ia = 0; 
for iel = 1:NT
    % ------- element information --------
    index = elem{iel};  Nv = length(index);    
    xK = centroid(iel,1); yK = centroid(iel,2); 
    hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv;  v2 = [2:Nv,1]; 
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    nodeT = [node(index,:); centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % ------- scaled monomials ----------
    m1 = @(x,y) 1+0*x;                gradm1 = @(x,y) [0+0*x, 0+0*x];
    m2 = @(x,y) (x-xK)./hK;           gradm2 = @(x,y) [1+0*x, 0+0*x]./hK;
    m3 = @(x,y) (y-yK)./hK;           gradm3 = @(x,y) [0+0*x, 1+0*x]./hK;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];   
    mc = {m1,m2,m3};  % cell 
    Gradmc = {gradm1, gradm2, gradm3};    
    
   % -------- transition matrix ----------
    D = m(x,y);   Dm{iel} = D;
    D1=m1(x,y);

    % --------- elliptic projection -----------
    % first term  = 0
    B = zeros(Nm, Nv);
    % second term   
    elem1 = [v1(:), v2(:)];
    for im = 1:Nm
        gradmc = Gradmc{im};
        F1 = 0.5*sum(gradmc(x(v1), y(v1)).*Ne, 2);
        F2 = 0.5*sum(gradmc(x(v2), y(v2)).*Ne, 2);
        F = [F1, F2];
        B(im, :) = accumarray(elem1(:), F(:), [Nv 1]);
    end  
    % constraint
    Bs = B;  Bs(1,:) = 1/Nv; 
    % consistency relation
    G = B*D;  Gs = Bs*D;      
    % --------- L2 projection -----------   
    H = zeros(Nm,Nm);
    for i = 1:Nm
        fun = @(x,y) repmat(mc{i}(x,y),1,Nm).*m(x,y);
        H(i,:) = integralTri(fun,2,nodeT,elemT);
    end 

    % --------- Piecewise constant projection-----------   
    fun0 = @(x,y) (m1(x,y).*m1(x,y));
    G0 = integralTri(fun0,2,nodeT,elemT);
    % consistency relation
    B0 = G0/D1; 
    Pis0 = G0\B0;  
    Pi0 = D1*Pis0; 
    
    fun_0 = @(x,y) repmat(m1(x,y),1,Nm).*m(x,y);
    H_0 = integralTri(fun_0,2,nodeT,elemT);  

    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi = D*Pis;   I = eye(size(Pi)); 
    AK  = Axi*(Pis'*G*Pis) + sagamaa*(I-Pi)'*(I-Pi);
    BK = c*(Pis'*H*Pis);  % Pis = Pi0s
    BK0 = (Pis0'*H_0*Pis);   
    
    A = reshape(AK',1,[]); 
    B = reshape(BK',1,[]); 
    B0 = reshape(BK0',1,[]);  % New add 2

    % --------- assembly index for ellptic projection -----------
    indexDof = index;  Ndof = Nv;  % local to global index
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ssA(ia+1:ia+Ndof^2) = A(:);
    ssB(ia+1:ia+Ndof^2) = B(:);
    ssB0(ia+1:ia+Ndof^2) = B0(:); % New add 3

    ia = ia + Ndof^2;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    Ph{iel} = Pis;Ph0{iel} = Pis0;
    elem2dof{iel} = indexDof;
end
kkA = sparse(ii,jj,ssA,N,N);
kkB = sparse(ii,jj,ssB,N,N);
kkB0 = sparse(ii,jj,ssB0,N,N);% New add 4

%% Apply Dirichlet boundary conditions
% bdNodeIdx = bdStruct.bdNodeIdx;
% isBdNode = false(N,1); isBdNode(bdNodeIdx) = true;
% freeDof = (~isBdNode);
bdStruct = setboundary(node,elem);
bdNodeIdx = bdStruct.bdNodeIdx;
bdDof = false(N,1); bdDof(bdNodeIdx) = true;
inDof = (~bdDof);

%% New Set solver
Nbd=sum(bdDof);Nin=sum(inDof);

A=kkA(inDof,inDof)+kkB(inDof,inDof);
M=kkB0(inDof,inDof);

dof=Nin;
item=5;
[eigf,lam] = eigs(A,M,item,'sm');
lam = diag(lam);
uh=zeros(dof,item);
Eigen =lam;
for arg=1:item
    Eigf = eigf(:,arg);
    Eigf=Eigf/abs((Eigf'*(A)*Eigf)^0.5);
    uh(inDof,arg)=Eigf(1:Nin);uh(bdDof,arg)=zeros(Nbd,1);
end
%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; info.Ph0 = Ph0;
%info.kk = kk; info.DofI = freeDof; % d.o.f.s for computing errors in the energy norm
info.D = Dm;
end