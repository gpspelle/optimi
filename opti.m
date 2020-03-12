%% Initialisation
clear all;  close all;  clc
% load the data
T=readtable('PositionsObjects.txt'); PositionsObjets=T.x+1i*T.y;  % positions defined by complex numbers
T=readtable('PositionsBoxes.txt');   PositionsBoxes=T.x+1i*T.y;

%% Plot example with object i in the box i (1 in 1, 2 in 2, etc.)
n = length (PositionsObjets);   % number of objects 

dist = zeros(n, n); % matrix of distances
for i=1:n
    for j=1:n
        dist(i,j) = norm(PositionsBoxes(i) - PositionsObjets(j));
    end
end

% unravelling a matrix to a vector
dist = reshape(dist.',1,[]);

% number of elements of the vector
size = length(dist);

% create the ineq constraints
Aineq = zeros(size, 1);
bineq = zeros(1);

% create the eq constraints
Aeq = zeros(size, 1);
beq = zeros(1);

% only one value in each row
i = 0;
while i < size
    row = zeros(size, 1);
    for j= 1:n
        row(i+j) = 1;
    end
    Aeq = [Aeq, row];
    beq = [beq; 1];
    i = i + n;
end

% only one value in each line
i = 1;
while i <= n
    row = zeros(size, 1);
    for j= 0:n-1
        row(i + n * j) = 1;
    end
    Aeq = [Aeq, row];
    beq = [beq; 1];
    i = i + 1;
end

% object 1 can't be in the last box
row = zeros(size, 1);
ind = to_ravel(n, 1, n);
row(ind) = 1;
Aeq = [Aeq, row];
beq = [beq; 0];

% object 2 can't be on the first box 
row = zeros(size, 1);
ind = to_ravel(1, 2, n);
row(ind) = 1;
Aeq = [Aeq, row];
beq = [beq; 0];

% object 1 is on the left of the object 2
for i = 1:n - 1   
   row = zeros(size, 1);
   ind = to_ravel(i, 1, n);
   row(ind) = 1;
   ind = to_ravel(i+1, 2 , n);
   row(ind) = -1;
   
   Aeq = [Aeq, row];
   beq = [beq; 0];
end

%  object 4 isn't on the right of the object 3
for i = 2:n-1
    row = zeros(size, 1);
    ind = to_ravel(i, 3, n);
    row(ind) = 1;
    for k = i+1: n-1
        ind = to_ravel(k, 4, n);
        row(ind) = 1;    
    end
    
    Aineq = [Aineq, row];
    bineq = [bineq; 1];
end

% object 3 can't be on the first box 
row = zeros(size, 1);
ind = to_ravel(1, 3, n);
row(ind) = 1;
Aeq = [Aeq, row];
beq = [beq; 0];

% object 4 can't be on the last box 
row = zeros(size, 1);
ind = to_ravel(n, 4, n);
row(ind) = 1;
Aeq = [Aeq, row];
beq = [beq; 0];

% if object 9 is in the first box object 7 is in the second box
rows = zeros(size, 1);
ind = to_ravel(1, 9, n);
row(ind) = 1;
ind = to_ravel(2, 7, n);
row(ind) = -1;
Aineq = [Aineq, row];
bineq = [bineq; 0];

% if object 9 is in the second box object 7 is in the first box
rows = zeros(size, 1);
ind = to_ravel(2, 9, n);
row(ind) = 1;
ind = to_ravel(1, 7, n);
row(ind) = -1;
Aineq = [Aineq, row];
bineq = [bineq; 0];

% if object 9 is in the last box then object 7 is in the penultimate box
rows = zeros(size, 1);
ind = to_ravel(n, 9, n);
row(ind) = 1;
ind = to_ravel(n-1, 7, n);
row(ind) = -1;
Aineq = [Aineq, row];
bineq = [bineq; 0];

% if object 9 is in the penultimate box then object 7 is in the last box
rows = zeros(size, 1);
ind = to_ravel(n, 7, n);
row(ind) = 1;
ind = to_ravel(n-1, 9, n);
row(ind) = -1;
Aineq = [Aineq, row];
bineq = [bineq; 14];

% boxe containing 7 must be close to boxe containing 9 for i = 2,...,n-1
for i = 2:n - 1
    row = zeros(size, 1);
    ind = to_ravel(i, 9, n);
    row(ind) = 1;
    ind = to_ravel(i-1, 7, n);
    row(ind) = -1;
    ind = to_ravel(i+1, 7, n);
    row(ind) = -1;

    Aineq = [Aineq, row];
    bineq = [bineq; 0];
 
end

% block the last solution we found
rows = zeros(size, 1);
row(to_ravel(1, 12, n)) = 1;
row(to_ravel(2, 13, n)) = 1;
row(to_ravel(3, 10, n)) = 1;
row(to_ravel(4, 9, n)) = 1;
row(to_ravel(5, 8, n)) = 1;
row(to_ravel(6, 7, n)) = 1;
row(to_ravel(7, 5, n)) = 1;
row(to_ravel(8, 14, n)) = 1;
row(to_ravel(9, 6, n)) = 1;
row(to_ravel(10, 11, n)) = 1;
row(to_ravel(11, 1, n)) = 1;
row(to_ravel(12, 2, n)) = 1;
row(to_ravel(13, 3, n)) = 1;
row(to_ravel(14, 4, n)) = 1;
row(to_ravel(15, 15, n)) = 1;
Aineq = [Aineq, row];
bineq = [bineq; n-1];

% remove the first line, we don't need them, it was created to allocate
% space in the begin
Aeq(:, 1) = [];
beq = beq(2:end);

Aineq(:, 1) = [];
bineq = bineq(2:end);

% transpose the matrix
Aineq = Aineq';
Aeq = Aeq';

% all variables are integers
intcon = 1:size;
intcon = intcon(:);


%options = optimoptions('intlinprog', 'AbsoluteGapTolerance', 0, 'CutGeneration', 'none', 'IntegerTolerance', 1e-6, 'CutMaxIterations', 50, 'Display', 'final', 'ConstraintTolerance', 1e-9);
%options = options + optimoptions('intlinprog', 'CutGeneration', 'none');
% solve the problem
x = intlinprog(dist, intcon, Aineq, bineq, Aeq, beq, zeros(size, 1), ones(size, 1));

% ravel the answer (vector -> matrix)
x = reshape(x, n, n);

x = round(x);

% get the indexes that are 1
ans_vec = zeros(n, 1);
for i = 1:n
    ans_vec(i) = find(x(i, :));
end

D = 0;

for j = 1:n
    D = D + abs(PositionsObjets(j) - PositionsBoxes(ans_vec(j)));
end

D

% transpose and plot
ans_vec'
PlotSolution (ans_vec, PositionsObjets, PositionsBoxes)

function x = to_ravel(i, j, l)
    x = (i-1) * l + j;
end
