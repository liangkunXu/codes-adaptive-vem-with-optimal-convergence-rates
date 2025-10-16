function [node,elem,bdEdge,brother,HB,belong] = bisect(node,elem,markedElem,bdEdge)
% BISECT refine the current triangulation by bisecting marked elements
% and minimal neighboring elements to get a conforming and shape regular
% triangulation. Newest vertex bisection is implemented.
% 
% [node,elem] = bisect(node,elem,markedElem);
% [node,elem,bdEdge] = bisect(node,elem,markedElem,bdEdge);
% [node,elem,bdEdge,HB,brother] = bisect(node,elem,markedElem,bdEdge);
%
% Input 
%    markedElem is a vector containing the index of elements to be bisected.
%       It could be a logical vector of length size(elem,1).
%    bdEdge stores boundary edges. This argument is optional.
%
% Output 
%  HB(:,1:3): hierarchical basis structure for new added nodes
%  	HB(:,1) the index of new added nodes
%  	HB(:,2:3) two paraent nodes of new added nodes
% HB is useful for the interpolation between two grids.
% See the function: nodeinterpolate
%
%  brother(:,1:2): two triangles coming from the same parent and
%      brother(t,1) < brother(t,2)
% brother is useful for the interpolation of elementwise function
% See function: eleminterpolate

%--------------------------------------------------------------------------
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details. 
%--------------------------------------------------------------------------
if iscell(elem) % transform to mat
    elem = cell2mat(elem);
end

HB = []; brother = []; 
if isempty(markedElem), return; end
%------------------- Construct auxiliary data structure -------------------
%[neighbor,elem2edge,edge]=auxstructure1(elem);
T=auxstructure1(elem);
neighbor=T.neighbor;
elem2edge=T.elem2edge;
edge=T.edge;
%[neighbor,elem2edge,edge] = auxstructurec(int32(elem));
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
%------------------- Add new nodes ----------------------------------------
isCutEdge = false(NE,1);
while sum(markedElem)>0
    isCutEdge(elem2edge(markedElem,1)) = true;
    refineNeighbor = neighbor(markedElem,1);
    markedElem = refineNeighbor(~isCutEdge(elem2edge(refineNeighbor,1)));
end
edge2newNode = zeros(NE,1,'uint32');
edge2newNode(isCutEdge) = N+1:N+sum(isCutEdge);%The number of the new node on the storage edge
HB = zeros(sum(isCutEdge),3,'uint32');
HB(:,1) = edge2newNode(isCutEdge);
HB(:,[2 3]) = edge(isCutEdge,[1 2]);
node(HB(:,1),:) = (node(HB(:,2),:) + node(HB(:,3),:))/2;
%------------------- Refine marked elements -------------------------------
Nb = 0; brother = zeros(3*NT,2,'uint32');belong=zeros(4*NT,1,'uint32');
%A unit is divided into two parts by at most three edges.
for k=1:2%refined mesh  after two bisects
    t = find(edge2newNode(elem2edge(:,1))>0);
    newNT = length(t);%The number of newly added elements = the number of edges cut.
    if (newNT == 0), break; end
    L = t; R = NT+1:NT+newNT;
    %The number of the ( cut ) element with the new node is assigned to L ( the element number on the left ), and the rest is assigned to R
    p1 = elem(t,1); p2 = elem(t,2); p3 = elem(t,3);
    p4 = edge2newNode(elem2edge(t,1));
    elem(L,:) = [p4, p1, p2];
    elem(R,:) = [p4, p3, p1];
	if nargin==4 %------------------ Refine boundary edges ----------------
   		bdEdge(R,[1 3]) = bdEdge(t,[2 1]);
   		bdEdge(L,[1 2]) = bdEdge(t,[3 1]);
        bdEdge(L,3) = 0;
        %The third edge of L must not be a boundary edge, so is the second edge of R. By default : bdEdge(L,2) = 0
	end
    brother(Nb+1:Nb+newNT,1) = L;%brother 
    brother(Nb+1:Nb+newNT,2) = R;
    elem2edge(L,1) = elem2edge(t,3);%
    elem2edge(R,1) = elem2edge(t,2);
    NT = NT + newNT; Nb = Nb + newNT;
    if(k==1); belong(L)=L;belong(R)=L;else belong(R)=belong(L);end% belong is only used for refined elements %added by me 
end
belong=belong(1:NT);
brother = brother(1:Nb,:);
if (nargin==3), bdEdge=[]; end