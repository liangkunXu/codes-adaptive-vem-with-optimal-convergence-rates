function [node,elem,bdEdge,brother,HB,belong,meshState] = bisect(node,elem,markedElem,bdEdge,meshState,Lambda) 
% BISECT refine marked elements by nonconforming newest-vertex bisection.
%
% until max(lambda(x)) <= Lambda.
%
% HB(:,1:3): [new node, first parent, second parent].
% belong(t): index of the old element containing new element t.
%
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details.

% begin Lambda-admissible nonconforming NVB
if nargin < 4
    bdEdge = [];
end
if nargin < 5
    meshState = [];
end
if nargin < 6 || isempty(Lambda)
    Lambda = 1;
end
if Lambda < 1 || Lambda ~= floor(Lambda)
    error('Lambda must be a positive integer.');
end

% VEM boundary nodes, including hanging nodes.
if isempty(meshState)
    tri = initialTriangles(elem);
    tri = orientTriangles(node,tri);
    meshState.tri = tri;
    meshState.nodeParent = zeros(size(node,1),2);
else
    if ~isfield(meshState,'tri') || ~isfield(meshState,'nodeParent')
        error('meshState must contain tri and nodeParent.');
    end
    tri = double(meshState.tri);
    if size(meshState.nodeParent,1) < size(node,1)
        meshState.nodeParent(end+1:size(node,1),1:2) = 0;
    elseif size(meshState.nodeParent,1) > size(node,1)
        error('meshState.nodeParent is inconsistent with node.');
    end
end

oldNT = size(tri,1);
if islogical(markedElem)
    if numel(markedElem) ~= oldNT
        error('A logical markedElem must have one entry per element.');
    end
    markedElem = find(markedElem);
else
    markedElem = markedElem(:);
end
markedElem = unique(double(markedElem),'stable');
if any(markedElem < 1) || any(markedElem > oldNT) || ...
        any(markedElem ~= floor(markedElem))
    error('markedElem contains an invalid element index.');
end

HB = zeros(0,3);
brother = zeros(0,2);
belong = (1:oldNT)';

% edgeMidpoint maps an edge to the node created at its midpoint.
edgeMidpoint = containers.Map('KeyType','char','ValueType','double');
for z = 1:size(meshState.nodeParent,1)
    parents = meshState.nodeParent(z,:);
    if all(parents > 0)
        key = edgeKey(parents(1),parents(2));
        if isKey(edgeMidpoint,key) && edgeMidpoint(key) ~= z
            error('Two midpoint nodes are associated with the same edge.');
        end
        edgeMidpoint(key) = z;
    end
end

% REFINE
for k = 1:numel(markedElem)
    bisectOne(markedElem(k));
end

% REFINE
numberOfCompletionSteps = 0;
while true
    elem = polygonalConnectivity();
    [isHanging,lambda,maxLambda,xHat] = globalIndices(elem);
    if maxLambda <= Lambda
        break;
    end

    numberOfCompletionSteps = numberOfCompletionSteps+1;
    if numberOfCompletionSteps > 100000
        error('MAKE_ADMISSIBLE did not terminate.');
    end

    elementHat = hangingElement(xHat,1:size(tri,1),elem);
    if isempty(elementHat)
        error('The maximal-index hanging node has no containing element.');
    end

    oppositeNodes = edgeInteriorNodes(tri(elementHat,2),tri(elementHat,3));
    if any(oppositeNodes == xHat)
        % Case A: one bisection is sufficient.
        bisectOne(elementHat);
    else
        % Case B: bisect E and then the child that still contains xHat as
        % a hanging node.
        children = bisectOne(elementHat);
        elemAfterFirstCut = polygonalConnectivity();
        childHat = hangingElement(xHat,children,elemAfterFirstCut);
        if isempty(childHat)
            error('Case B could not locate the required child element.');
        end
        bisectOne(childHat);
    end
end

meshState.tri = tri;
meshState.isHanging = isHanging;
meshState.lambda = lambda;
meshState.maxLambda = maxLambda;

bdEdge = [];

    function children = bisectOne(t)
        % One newest-vertex bisection, without refining a neighbour.
        p1 = tri(t,1);
        p2 = tri(t,2);
        p3 = tri(t,3);
        key = edgeKey(p2,p3);

        if isKey(edgeMidpoint,key)
            p4 = edgeMidpoint(key);
        else
            p4 = size(node,1)+1;
            node(p4,:) = 0.5*(node(p2,:)+node(p3,:));
            meshState.nodeParent(p4,:) = sort([p2,p3]);
            edgeMidpoint(key) = p4;
            HB(end+1,:) = [p4,p2,p3];
        end

        left = t;
        right = size(tri,1)+1;
        oldOwner = belong(t);
        tri(left,:) = [p4,p1,p2];
        tri(right,:) = [p4,p3,p1];
        belong(right,1) = oldOwner;
        brother(end+1,:) = [left,right];
        children = [left,right];
    end

    function polygons = polygonalConnectivity()
        % Insert all descendant midpoints on the three geometric sides.
        polygons = cell(size(tri,1),1);
        for t = 1:size(tri,1)
            corners = tri(t,:);
            polygon = zeros(1,0);
            for j = 1:3
                a = corners(j);
                b = corners(mod(j,3)+1);
                polygon = [polygon,a,edgeInteriorNodes(a,b)]; %#ok<AGROW>
            end
            polygons{t} = polygon;
        end
    end

    function interior = edgeInteriorNodes(a,b)
        % Nodes on [a,b], ordered from a to b, without the endpoints.
        key = edgeKey(a,b);
        if ~isKey(edgeMidpoint,key)
            interior = zeros(1,0);
            return;
        end
        midpoint = edgeMidpoint(key);
        interior = [edgeInteriorNodes(a,midpoint),midpoint, ...
                    edgeInteriorNodes(midpoint,b)];
    end

    function [hanging,index,largest,xLargest] = globalIndices(polygons)
        % Definition 2.1.
        hanging = false(size(node,1),1);
        for t = 1:size(tri,1)
            localHanging = setdiff(polygons{t},tri(t,:));
            hanging(localHanging) = true;
        end

        index = zeros(size(node,1),1);
        for z = 1:size(node,1)
            if ~hanging(z)
                continue;
            end
            parents = meshState.nodeParent(z,:);
            if any(parents == 0) || any(parents >= z)
                error('Invalid parent data for hanging node %d.',z);
            end
            index(z) = max(index(parents))+1;
        end
        [largest,xLargest] = max(index);
    end

    function tFound = hangingElement(z,candidates,polygons)
        tFound = [];
        for q = 1:numel(candidates)
            t = candidates(q);
            if any(polygons{t} == z) && ~any(tri(t,:) == z)
                tFound = t;
                return;
            end
        end
    end
% end Lambda-admissible nonconforming NVB
end

function key = edgeKey(a,b) 
edge = sort([a,b]);
key = sprintf('%d_%d',edge(1),edge(2));
end

function tri = initialTriangles(elem) 
if iscell(elem)
    if any(cellfun(@numel,elem) ~= 3)
        error(['A polygonal mesh with hanging nodes requires meshState ', ...
               'from the preceding call.']);
    end
    tri = zeros(numel(elem),3);
    for t = 1:numel(elem)
        tri(t,:) = elem{t}(:)';
    end
else
    if size(elem,2) ~= 3
        error('The initial mesh must consist of triangles.');
    end
    tri = double(elem);
end
end

function tri = orientTriangles(node,tri)
for t = 1:size(tri,1)
    p = node(tri(t,:),:);
    signedDoubleArea = det([p(2,:)-p(1,:);p(3,:)-p(1,:)]);
    if signedDoubleArea == 0
        error('Initial element %d is degenerate.',t);
    elseif signedDoubleArea < 0
        tri(t,[2,3]) = tri(t,[3,2]);
    end
end
end
