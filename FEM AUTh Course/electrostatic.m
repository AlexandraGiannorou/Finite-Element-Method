clear all; 

%% dimensions
w = 1e-2;   %[m]
d = 0.5e-2; %[m]
h = 1e-3;   %[m]
A = 5*w;    %[m]
B = 5*w;    %[m]
er = 1;     %dielectric constant
%% defining the geometry and meshing
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
[p,e,t] = refinemesh(d1,p,e,t);

% ie p is a 2xN array
% 2 rows for the 2 coordinates of each point
% N columns for the number of the mesh points

% t is the array that includes the triangle information, what points create
% a triangle

%[p,e,t] = refinemesh(d1,p,e,t); % optimizes the mesh by dividing an element 
%pdeplot(p,e,t); %mesh visualization 
axis equal;
axis tight; 

%% Applying Boundary Conditions
%on the capacitor plates Dirichlet BC are applied, thus the value is known

Nn = size(p,2); %number of nodes
Nel = size(t,2); %number of elements
Ned = size(e,2); %number of edges

% e array includes information about the edges on the outter boundaries, or
% the boundaries between reagions
% e(1, id), e(2, id) are the two points of the edge
% e(6, id), e(7, id) are the two regions on the two sides of the edge


node_id = ones (Nn,1);
X0 = zeros(Nn, 1);
V0 = 100;  %[Volt]

for id = 1:Ned  %runs a loop scanning all the edges 
    n1 = e(1,id); %gets the two nodes of each edge
    n2 = e(2,id);
    r1 = e(6,id); %gets the two regions of each edge
    r2 = e(7,id);

    % gets the cordinates of the two points of the edges 
    x1 = p(1, n1);
    y1 = p(2, n1);
    x2 = p(1, n2);
    y2 = p(2, n2);
   
    if ((r1==0 || r2 == 0) && (y1>0 && y1<B/4 && abs(x1)<1.2*w/2 && y2>0 && y2<B/4 && abs(x2) < 1.2 *w/2 ) ) %sets the dirichlet boundaries for the capacitor plates, inside a region that inlcudes the plates, forces the voltage values on the points of the edges between two regions
        node_id(n1) = 0;
        node_id(n2) = 0; 
        X0(n1) = V0/2;
        X0(n2) = V0/2;

    end

     if ((r1==0 || r2 == 0) && (y1<0 && y1>-B/4 && abs(x1)<1.2*w/2 && y2<0 && y2>-B/4 && abs(x2) < 1.2 *w/2 ) )
        node_id(n1) = 0;
        node_id(n2) = 0; 
        X0(n1) = -V0/2;
        X0(n2) = -V0/2;
     end

end 

% for in = 1: Nn
%     x = p(1, in); 
%     y = p(2, in);
%     text(x,y, num2str(in));    %in the coordinates of each point, it goes and writes as a string the node id, if it is zero or one, to check if the dirichlet bc are defined correctly
% end


%% seperating the known from the unknown 
ic = 0; % counter
index = zeros(Nn,1);
for in = 1: Nn                 %scans all the nodes
    if (node_id(in)==1)        % if the node id is 1, meaning the node is unknown
                               % increment the counter by 1, else increment
                               % enter the loop again
        ic = ic + 1;             % the counter counts the number of the unknowns
        index(in) = ic;        % at the node index of the node with id==1 

    end    
end

Nunknowns = ic; 

% for ie = 1: Nel
%     n(1:3) = t(1:3, ie);  % n includes the three nodes of each triangular element
%     x(1:3) = p(1, n(1:3)); % x get the coordinates of each node of each  triangular element
%     y(1:3) = p(2, n(1:3));
%     region = t(4, ie);
%     xg = sum(x)/3; 
%     yg = sum(y)/3; 
%     text (xg, yg, num2str(region));
% end

%% Finite - Element - Methos 

% Initializing Stiffness Matrix S 
S = spalloc (Nunknowns, Nunknowns, 7*Nunknowns);
B = zeros (Nunknowns, 1);

for ie = 1: Nel
     
    n(1:3) = t(1:3, ie);  % n includes the three nodes of each triangular element
    x(1:3) = p(1, n(1:3)); % x get the coordinates of each node of each  triangular element
    y(1:3) = p(2, n(1:3));
    region = t(4, ie);
    %Simplex Coordinates Coefficients
    
    D = det( [1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3) ] );
    
    b(1) = (y(2)-y(3))/D;
    c(1) = (x(3)-x(2))/D;

    b(2) = (y(3)-y(1))/D;
    c(2) = (x(1)-x(3))/D;

    b(3) = (y(1)-y(2))/D;
    c(3) = (x(2)-x(1))/D;

    Ae = abs(D)/2;

    for i = 1:3
        for j = 1:3

            Se(i,j) = er * (b(i)*b(j) + c(i)*c(j)) *Ae; %Stiffness matrix for each triangular element

            %% Assembly
            % The total unknown function is a sum of φΝ for all the nodes of the grid
            % Np total is a sum of all the basis functions the p node belongs to 
            % Thus the Galerkin Formulation for the whole grid will be the sum of all
            % elements
            % The Stiffness Matrix of this formulation will be the total Stifness
            % Matrix,for all the nodes of the grid
            % The total Stiffness Matrix will be sparse, with non zero terms only when
            % two nodes are close to each other: Only from the elements in which both
            % nodes belong
            % Then the total Stiffness Matrix will be calculated by the Element Stiffness 
            % of those two elements.

            if (node_id(n(i))==1 )
                if (node_id(n(j))==1 )
                    S( index(n(i)) , index(n(j))) =  S(index(n(i)) , index(n(j))) + Se(i,j); 
                else
                    B(index (n(i))) = B(index(n(i))) - Se(i,j)*X0(n(j));
                end
            end

        end
    end

end

X = S\B; %solving Sx=B with Gauss-Jordan method
for in = 1:Nn
    if (node_id(in) == 1) 
        X0(in) = X( index(in) );
    end
end

% pdeplot(p,e,t,'xydata',X0, 'contour', 'on',  'mesh', 'on');
% axis equal
% axis tight;
% hold on;
% colormap jet;

figure
pdeplot(p,e,t, 'XYData', X0, 'mesh', 'off', 'contour', 'on');
colormap jet
hold on

% contour lines dashed & thicker
contourLines = findobj(gca, 'Type','Line');
set(contourLines, 'LineStyle','--', 'LineWidth', 1.5)

% overlay the mesh
pdeplot(p,e,t, 'mesh', 'on', 'contour', 'off');

axis equal
axis tight

