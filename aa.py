from __future__ import print_function
from dolfin import*
import numpy as np
#creat mesh

width = 1.0
height = 0.5
mesh = RectangleMesh(Point(0, 0), Point(width, height), 8, 4)

#define the function space
P1 = VectorFunctionSpace(mesh, "CG", 1)
PN = FunctionSpace(mesh, "Nedelec 1st kind H(curl)", 1)

V = FunctionSpace(mesh, "Nedelec 1st kind H(curl)", 1)
#define boundary conditons

#u0=np.array([0.0,0.0])

#def u0_boundary(x,on_boundary):
   # return on_boundary
#bc=DirichletBC(V,u0,u0_boundary)

def curl_t(w):
    return Dx(w[1], 0) - Dx(w[0], 1)

zero = Constant((0.0, 0.0))
#define boundary conditons

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Boundary condition
bc = DirichletBC(V, zero, DirichletBoundary())

#define f
f = Expression(("-1.0",  "1.0"), degree=1)
#f=Constant((-4, 4))
#g1=Expression("-4*x[0]",degree=2)
#g2=Expression("4*x[1]",degree=2)
#f=np.array([-4,4])

#define the test and trial functions
u=TrialFunction(V)
v=TestFunction(V)
#define variational problem 
a=inner(curl(u),curl(v))-inner(u,v)
L=inner(f,v)
T=Function(V)
solve(a==L,T,bc)
#A, b = assemble_system(a, L, bc)

#Compute solution
#T=Function(V)
#solve(inner(curl(v), curl(u))*dx -inner(u,v)==f*v*dx, T, bc)
#solve(curl(v)*curl(u)*dx+inner(u,v)*dx==f*v*dx,u,bc)


#plot(T)


plot(mesh)
interactive()
