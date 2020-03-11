%% Initialisation
clear all;  close all;  clc
% load the data
T=readtable('PositionsObjects.txt'); PositionsObjets=T.x+1i*T.y;  % positions defined by complex numbers
T=readtable('PositionsBoxes.txt');   PositionsBoxes=T.x+1i*T.y;

%% Plot example with object i in the box i (1 in 1, 2 in 2, etc.)
n = length (PositionsObjets);   % number of objects 

l_boxes = length(PositionsBoxes);
l_objets = length(PositionsObjets);

dist_mat = zeros(l_boxes, l_objets);
for i=1:l_boxes
    for j=1:l_objets
        dist_mat(i,j) = norm(PositionsBoxes(i) - PositionsObjets(j));
    end
end

dist_vec = reshape(dist_mat.',1,[]);
size_vec = length(dist_vec);
Aeq = zeros(size_vec, 1);
beq = zeros(1);
Aineq = zeros(size_vec, 1);
bineq = zeros(1);


i = 0;
while i < size_vec
    row = zeros(size_vec, 1);
    for j= 1:l_objets
        row(i+j) = 1;
    end
    %fprintf("%d\n", i);
    Aeq = [Aeq, row];
    beq = [beq; 1];
    i = i + l_objets;
end

i = 1;
while i <= l_boxes
    row = zeros(size_vec, 1);
    for j= 0:l_boxes-1
        row(i + l_boxes * j) = 1;
    end
    Aeq = [Aeq, row];
    beq = [beq; 1];
    i = i + 1;
end

% object 1 can't be in the last box
row = zeros(size_vec, 1);
ind = to_ravel(l_boxes, 1, l_objets);
row(ind) = 1;

Aeq = [Aeq, row];
beq = [beq; 0];

% object 2 can't be on the first box 
row = zeros(size_vec, 1);
ind = to_ravel(1, 2, l_objets);
row(ind) = 1;
Aeq = [Aeq, row];
beq = [beq; 0];

% object 1 is on the left of the object 2
for i = 1:l_boxes - 1   
   row = zeros(size_vec, 1);
   ind = to_ravel(i, 1, l_objets);
   row(ind) = 1;
   ind = to_ravel(i+1, 2 , l_objets);
   row(ind) = -1;
   
   Aeq = [Aeq, row];
   beq = [beq; 0];
end

% % object 4 isn't on the right of the object 3
% for l = l_boxes : -1 : 2
%     row = zeros(size_vec, 1);
%     ind = to_ravel(l, 3, l_objets);
%     row(ind) = 1;
%     for k = l-1: -1 : 1
%         ind = to_ravel(k, 4, l_objets);
%         row(ind) = -1;
%     end
%    
%     Aeq = [Aeq, row];
%     beq = [beq; 0];
% end

for i = 2:l_boxes-1
    row = zeros(size_vec, 1);
    ind = to_ravel(i, 3, l_objets);
    row(ind) = 1;
    for k = 1: l_boxes
        ind = to_ravel(i+k, 4, l_objets);
        if ind < size_vec
            row(ind) = 1;    
        end
    end
    
    Aineq = [Aineq, row];
    bineq = [bineq; 1];
end
% object 3 can't be on the first box 
row = zeros(size_vec, 1);
ind = to_ravel(1, 3, l_objets);
row(ind) = 1;
Aeq = [Aeq, row];
beq = [beq; 0];

% object 4 can't be on the last box 
row = zeros(size_vec, 1);
ind = to_ravel(l_objets, 4, l_objets);
row(ind) = 1;
Aeq = [Aeq, row];
beq = [beq; 0];

% if object 9 is in the first box object 7 is in the second box
rows = zeros(size_vec, 1);
ind = to_ravel(1, 9, l_objets);
row(ind) = 1;
ind = to_ravel(2, 7, l_objets);
row(ind) = -1;

Aeq = [Aeq, row];
beq = [beq; 0];

% if object 9 is in the last box then object 7 is in the penultimate box
rows = zeros(size_vec, 1);
ind = to_ravel(l_objets, 9, l_objets);
row(ind) = 1;
ind = to_ravel(l_objets-1, 7, l_objets);
row(ind) = -1;

Aeq = [Aeq, row];
beq = [beq; 0];

% boxe containing 7 must be close to boxe containing 9 for i = 2,...,n-1
for i = 2:l_objets - 1
    rows = zeros(size_vec, 1);
    ind = to_ravel(i, 9, l_objets);
    row(ind) = 1;
    ind = to_ravel(i-1, 7, l_objets);
    row(ind) = -1;
    ind = to_ravel(i+1, 7, l_objets);
    row(ind) = -1;

    Aeq = [Aeq, row];
    beq = [beq; 0];
 
end

%Aeq(1, :) =[];
Aeq(:, 1) = [];
beq = beq(2:end);

Aineq(:, 1) = [];
bineq = bineq(2:end);

Aineq = Aineq';
Aeq = Aeq';

intcon = 1:size_vec;
intcon = intcon(:);

x = intlinprog(dist_vec, intcon, Aineq, bineq, Aeq, beq, zeros(size_vec, 1), ones(size_vec, 1));

x = reshape(x, l_objets, l_boxes);

ans_vec = zeros(l_boxes, 1);
for i = 1:l_boxes
    ans_vec(i) = find(x(i, :));
end
ans_vec'
PlotSolution (ans_vec, PositionsObjets, PositionsBoxes)


function x = to_ravel(i, j, l)
    x = (i-1) * l + j;
end
