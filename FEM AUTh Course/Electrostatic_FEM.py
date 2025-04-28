from fenics import *
from mshr   import *
import matplotlib.pyplot as plt

# — your MATLAB parameters —
w = 1e-2
d = 0.5e-2
h = 1e-3
A = 5*w
B = 5*w

# — 1) main domain R1 (A×B) centered at (0,0) —
R1 = Rectangle(Point(-A/2, -B/2),
               Point( A/2,  B/2))

R2 = Rectangle (Point(-w/2, d/2),
                Point(w/2, d/2+h ) ) # top plate

R3 = Rectangle (Point(-w/2, -d/2-h),
                
                Point(w/2, -d/2) ) # bottom plate

# — 3) CSG: subtract plates from the dielectric domain —
domain = R1 -R2 -R3

# — 4) mesh generation (initmesh) —
mesh = generate_mesh(domain, 64)

# (optional) single refinement, if you want exactly one call to refinemesh:
# mesh = refine(mesh)

# — 5) plot just the mesh to compare with pdeplot —
plt.figure(figsize=(6,6))
plot(mesh, linewidth=0.5)
plt.axis('equal')
plt.title("Mesh of A×B domain with two w×h plates removed")
plt.tight_layout()
plt.show()

# — 6) counts to match Nn, Nel, Ned —
Nn  = mesh.num_vertices()
Nel = mesh.num_cells()
mesh.init(1)              # build edge (facet) connectivity
Ned = mesh.num_entities(1)

#print(f"Nn  = {Nn}")
#print(f"Nel = {Nel}")
#print(f"Ned = {Ned}")
