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

i = 0;
while i < size_vec
    row = zeros(size_vec, 1);
    for j= 1:l_objets
        row(i+j) = 1;
    end
    %fprintf("%d\n", i);
    Aeq = [Aeq, row];
    i = i + l_objets;
end

i = 1;
while i <= l_boxes
    row = zeros(size_vec, 1);
    for j= 0:l_boxes-1
        row(i + l_boxes * j) = 1;
    end
    %fprintf("%d\n", i);
    Aeq = [Aeq, row];
    i = i + 1;
end

%Aeq(1, :) =[];
Aeq(:, 1) =[];
Aeq = Aeq';

[rows, cols] = size(Aeq);
beq = ones(rows, 1);


x = intlinprog(dist_vec, length(dist_vec), [], [], Aeq, beq, zeros(size_vec, 1), ones(size_vec, 1));

x = reshape(x, l_objets, l_boxes);

ans_vec = zeros(l_boxes, 1);
for i = 1:l_boxes
    ans_vec(i) = find(x(i, :));
end
ans_vec'
PlotSolution (ans_vec, PositionsObjets, PositionsBoxes)
