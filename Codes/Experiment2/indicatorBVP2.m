function eta=indicatorBVP2(node,elem,uh0,uh,info,log13)
% This function returns the local error indicator of solving the laplace
% eigenproblem
% using virtual element method in V1
% Copyright (C) Xu.
%case 1
c=0;
A13=161.4476387975881;A24=1;
log24=~log13;
%% Pis and chi
% Pis
Ph = info.Ph; 
Ph0 = info.Ph0;
% elementwise numerical d.o.f.s
indexg = info.elem2dof; % elemenwise global index
chiu = cellfun(@(id) uh(id), indexg, 'UniformOutput', false); 
chiu0 = cellfun(@(id) uh0(id), indexg, 'UniformOutput', false); 
% coefficient of elliptic projection: Ph{iel}*chi{iel}
au = cellfun(@mtimes, Ph, chiu, 'UniformOutput', false);
au0 = cellfun(@mtimes, Ph0, chiu0, 'UniformOutput', false);
%% auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid;  diameter = aux.diameter;
% auxstructure
auxT = auxstructure(node,elem);
edge = auxT.edge; 
elem2edge = auxT.elem2edge; edge2elem = auxT.edge2elem;
% number
NT = size(elem,1); NE = size(edge,1);

%% elementwise residuals
eta1 = zeros(NT,1); 
% squared
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    % nodeT = [node(index,:);aux.centroid(iel,:)];
    % elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    % scaled monomials
    m1 = @(x,y)  1+0*x;
    m2 = @(x,y) (x-xK)./hK;
    m3 = @(x,y) (y-yK)./hK;
    
    Piuh = @(x,y) ( au0{iel}*1-c*(au{iel}(1)*m1(x,y)+au{iel}(2)*m2(x,y)+au{iel}(3)*m3(x,y)) ).^2;
    %eta1(iel) = hK^2*integralTri(Piuh,2,nodeT,elemT);
    eta1(iel) = hK^2*integralTri(Piuh,2,node(index,:),[1 2 3]);
end
elemRes = eta1;

%% elementwise edge jumps
% scaled norm vectors he*ne
e = node(edge(:,2),:)-node(edge(:,1),:);
Ne = [-e(:,2), e(:,1)];
% edgeJump
edgeJump = zeros(NE,1);
gradm = @(hK) [0 0; 1/hK 0; 0 1/hK];    
for s = 1:NE
    % element information
    k1 = edge2elem(s,1);  k2 = edge2elem(s,2);
    if k1==k2, continue; end  
    h1 = diameter(k1);   h2 = diameter(k2);
    % grad of Pi(uh)
    gradLu = au{k1}'*gradm(h1)*(A13*log13(k1)+A24*log24(k1)); 
    gradRu = au{k2}'*gradm(h2)*(A13*log13(k2)+A24*log24(k2));
    % jump of grad(Pi(uh))
    Jumpu = (gradLu-gradRu);
    % edgeJump
    edgeJump(s) = 0.5*(dot(Jumpu,Ne(s,:)))^2; % ce = 0.5
end
elemJump = cellfun(@(ie) sum(edgeJump(ie)), elem2edge);
%% Local error indicator  
%eta = sqrt(elemRes+elemJump);
eta = elemRes+elemJump;