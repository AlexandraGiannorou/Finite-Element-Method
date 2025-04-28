clear all; 

%dimensions
w = 1e-2;
d = 0.5e-2;
h = 1e-3; 
A = 5*w;
B = 5*w;

%defining the geometry 
gd = [  3          3           3;
        4          4           4;
      -A/2        w/2         w/2;
       A/2        w/2         w/2;
       A/2       -w/2        -w/2;
      -A/2       -w/2        -w/2;
      -B/2        d/2        -d/2;
      -B/2      d/2+h      -d/2-h;
       B/2      d/2+h      -d/2-h;
       B/2        d/2        -d/2
       ];

sf = 'R1-R2-R3'; %removes regions R2 and R3 from the mesh grid
ns = [82 82 82;
      49 50 51]; % every column corresponds to one variable of the set formula sf, 82, 49, 50, 51 are the ASCII codes for 'R', '1','2','3'
d1 = decsg(gd, sf, ns);
[p,e,t] = initmesh(d1); %Delaunay traingulation algorithm, creates the mesh. Returns the 3 matrices p,e, t that contain all the information for the mesh grid

% ie p is a 2xN array
% 2 rows for the 2 coordinates of each point
% N columns for the number of the mesh points

% t is the array that includes the triangle information, what points create
% a triangle

%[p,e,t] = refinemesh(d1,p,e,t); % optimizes the mesh by dividing an element 
pdeplot(p,e,t); %mesh visualization 
axis equal;
axis tight; 

%on the capacitor plates Dirichlet BC are applied, thus the value is known

Nn = size(p,2); %number of nodes
Nel = size(t,2); %number of elements
Ned = size(e,2); %number of edges

% e array includes information about the edges on the outter boundaries, or
% the boundaries between reagions
% e(1, id), e(2, id) are the two points of the edge
% e(6, id), e(7, id) are the two regions on the two sides of the edge


node_id = ones (Nn,1);

for id = 1:Ned  %runs a loop scanning all the edges 
    n1 = e(1,id); %gets the two nodes of each edge
    n2 = e(2,id);
    r1 = e(6,id); %gets the two regions of each edge
    r2 = e(7,id);

    if (r1==0 || r2 == 0) 
        node_id(n1) = 0;
        node_id(n2) = 0; 
    end
end 

for in = 1: Nn
    x = p(1, in); 
    y = p(2, in);
    text(x,y, num2str(node_id(in))); %in the coordinates of each point, it goes and writes as a string the node id, if it is zero or one, to check if the dirichlet bc are defined correctly
end

