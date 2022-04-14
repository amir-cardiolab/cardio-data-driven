from fenics import *
from mshr import *
import numpy 

time_input = numpy.array([0,
    0.0096,
    0.0192,
    0.0288,
    0.0384,
    0.0480,
    0.0576,
    0.0672,
    0.0768,
    0.0864,
    0.0960,
    0.1056,
    0.1152,
    0.1247,
    0.1343,
    0.1439,
    0.1535,
    0.1631,
    0.1727,
    0.1823,
    0.1919,
    0.2015,
    0.2111,
    0.2207,
    0.2303,
    0.2399,
    0.2495,
    0.2591,
    0.2687,
    0.2783,
    0.2879,
    0.2975,
    0.3071,
    0.3167,
    0.3263,
    0.3359,
    0.3455,
    0.3551,
    0.3646,
    0.3742,
    0.3838,
    0.3934,
    0.4030,
    0.4126,
    0.4222,
    0.4318,
    0.4414,
    0.4510,
    0.4606,
    0.4702,
    0.4798,
    0.4894,
    0.4990,
    0.5086,
    0.5182,
    0.5278,
    0.5374,
    0.5470,
    0.5566,
    0.5662,
    0.5758,
    0.5854,
    0.5949,
    0.6045,
    0.6141,
    0.6237,
    0.6333,
    0.6429,
    0.6525,
    0.6621,
    0.6717,
    0.6813,
    0.6909,
    0.7005,
    0.7101,
    0.7197,
    0.7293,
    0.7389,
    0.7485,
    0.7581,
    0.7677,
    0.7773,
    0.7869,
    0.7965,
    0.8061,
    0.8157,
    0.8253,
    0.8348,
    0.8444,
    0.8540,
    0.8636,
    0.8732,
    0.8828,
    0.8924,
    0.9020,
    0.9116,
    0.9212,
    0.9308,
    0.9404,
    0.9500,
    0.9596,
    0.9692,
    0.9788,
    0.9884,
    0.9980,
    1.0076,
    1.0172,
    1.0268,
    1.0364,
    1.0460,
    1.0556,
    1.0652,
    1.0747,
    1.0843,
    1.0939,
    1.1035,
    1.1131,
    1.1227,
    1.1323,
    1.1419,
    1.1515,
    1.1611,
    1.1707,
    1.1803,
    1.1899,
    1.1995,
    1.2091,
    1.2187,
    1.2283,
    1.2379,
    1.2475,
    1.2571,
    1.2667,
    1.2763,
    1.2859,
    1.2955,
    1.3051,
    1.3146,
    1.3242,
    1.3338,
    1.3434,
    1.3530,
    1.3626,
    1.3722,
    1.3818,
    1.3914,
    1.4010,
    1.4106,
    1.4202,
    1.4298,
    1.4394,
    1.4490,
    1.4586,
    1.4682,
    1.4778,
    1.4874,
    1.4970,
    1.5066,
    1.5162,
    1.5258,
    1.5354,
    1.5449,
    1.5545,
    1.5641,
    1.5737,
    1.5833,
    1.5929,
    1.6025,
    1.6121,
    1.6217,
    1.6313,
    1.6409,
    1.6505,
    1.6601,
    1.6697,
    1.6793,
    1.6889,
    1.6985,
    1.7081,
    1.7177,
    1.7273,
    1.7369,
    1.7465,
    1.7561,
    1.7657,
    1.7753,
    1.7848,
    1.7944,
    1.8040,
    1.8136,
    1.8232,
    1.8328,
    1.8424,
    1.8520,
    1.8616,
    1.8712,
    1.8808,
    1.8904,
    1.9000,
    1.9096,
    1.9192,
    1.9288,
    1.9384,
    1.9480,
    1.9576,
    1.9672,
    1.9768,
    1.9864,
    1.9960,
    2.0056,
    2.0152,
    2.0247,
    2.0343,
    2.0439,
    2.0535,
    2.0631,
    2.0727,
    2.0823,
    2.0919,
    2.1015,
    2.1111,
    2.1207,
    2.1303,
    2.1399,
    2.1495,
    2.1591,
    2.1687,
    2.1783,
    2.1879,
    2.1975,
    2.2071,
    2.2167,
    2.2263,
    2.2359,
    2.2455,
    2.2551,
    2.2646,
    2.2742,
    2.2838,
    2.2934,
    2.3030,
    2.3126,
    2.3222,
    2.3318,
    2.3414,
    2.3510,
    2.3606,
    2.3702,
    2.3798,
    2.3894,
    2.3990,
    2.4086,
    2.4182,
    2.4278,
    2.4374,
    2.4470,
    2.4566,
    2.4662,
    2.4758,
    2.4854,
    2.4949,
    2.5045,
    2.5141,
    2.5237,
    2.5333,
    2.5429,
    2.5525,
    2.5621,
    2.5717,
    2.5813,
    2.5909,
    2.6005,
    2.6101,
    2.6197,
    2.6293,
    2.6389,
    2.6485,
    2.6581,
    2.6677,
    2.6773,
    2.6869,
    2.6965,
    2.7061,
    2.7157,
    2.7253,
    2.7348,
    2.7444,
    2.7540,
    2.7636,
    2.7732,
    2.7828,
    2.7924,
    2.8020,
    2.8116,
    2.8212,
    2.8308,
    2.8404,
    2.8500,          
       ])

V_vel = numpy.array([17.8696,
   18.6069,
   20.8597,
   23.5008,
   26.5045,
   30.9845,
   36.2607,
   40.3874,
   44.6993,
   48.5726,
   51.5080,
   53.9016,
   55.0639,
   54.4802,
   53.1004,
   51.7180,
   50.4252,
   49.5198,
   49.4959,
   50.0309,
   50.6385,
   51.4014,
   52.4245,
   53.1627,
   53.2386,
   52.8589,
   52.2846,
   51.5976,
   50.7110,
   49.5496,
   48.1604,
   46.6457,
   45.0022,
   43.1104,
   40.8619,
   38.4350,
   36.2394,
   34.4397,
   32.9736,
   31.8165,
   31.1723,
   31.2849,
   31.9983,
   32.8593,
   33.5215,
   33.9174,
   34.0250,
   33.8807,
   33.5701,
   33.1802,
   32.7637,
   32.2714,
   31.6305,
   30.8181,
   29.9247,
   29.0551,
   28.2650,
   27.5183,
   26.7750,
   26.0506,
   25.3978,
   24.8311,
   24.3567,
   23.9744,
   23.6484,
   23.3344,
   22.9896,
   22.6543,
   22.4350,
   22.3966,
   22.4102,
   22.2976,
   22.0186,
   21.7907,
   21.6576,
   21.5586,
   21.4383,
   21.2787,
   21.0961,
   20.9084,
   20.7309,
   20.5773,
   20.4186,
   20.1984,
   19.8613,
   19.3997,
   18.9730,
   18.6786,
   18.4986,
   18.4115,
   18.3970,
   18.4311,
   18.4934,
   18.5608,
   18.6129,
   18.6274,
   18.5822,
   18.4550,
   18.2255,
   17.8696,
   18.6069,
   20.8597,
   23.5008,
   26.5045,
   30.9845,
   36.2607,
   40.3874,
   44.6993,
   48.5726,
   51.5080,
   53.9016,
   55.0639,
   54.4802,
   53.1004,
   51.7180,
   50.4252,
   49.5198,
   49.4959,
   50.0309,
   50.6385,
   51.4014,
   52.4245,
   53.1627,
   53.2386,
   52.8589,
   52.2846,
   51.5976,
   50.7110,
   49.5496,
   48.1604,
   46.6457,
   45.0022,
   43.1104,
   40.8619,
   38.4350,
   36.2394,
   34.4397,
   32.9736,
   31.8165,
   31.1723,
   31.2849,
   31.9983,
   32.8593,
   33.5215,
   33.9174,
   34.0250,
   33.8807,
   33.5701,
   33.1802,
   32.7637,
   32.2714,
   31.6305,
   30.8181,
   29.9247,
   29.0551,
   28.2650,
   27.5183,
   26.7750,
   26.0506,
   25.3978,
   24.8311,
   24.3567,
   23.9744,
   23.6484,
   23.3344,
   22.9896,
   22.6543,
   22.4350,
   22.3966,
   22.4102,
   22.2976,
   22.0186,
   21.7907,
   21.6576,
   21.5586,
   21.4383,
   21.2787,
   21.0961,
   20.9084,
   20.7309,
   20.5773,
   20.4186,
   20.1984,
   19.8613,
   19.3997,
   18.9730,
   18.6786,
   18.4986,
   18.4115,
   18.3970,
   18.4311,
   18.4934,
   18.5608,
   18.6129,
   18.6274,
   18.5822,
   18.4550,
   18.2255,
   17.8696,
   18.6069,
   20.8597,
   23.5008,
   26.5045,
   30.9845,
   36.2607,
   40.3874,
   44.6993,
   48.5726,
   51.5080,
   53.9016,
   55.0639,
   54.4802,
   53.1004,
   51.7180,
   50.4252,
   49.5198,
   49.4959,
   50.0309,
   50.6385,
   51.4014,
   52.4245,
   53.1627,
   53.2386,
   52.8589,
   52.2846,
   51.5976,
   50.7110,
   49.5496,
   48.1604,
   46.6457,
   45.0022,
   43.1104,
   40.8619,
   38.4350,
   36.2394,
   34.4397,
   32.9736,
   31.8165,
   31.1723,
   31.2849,
   31.9983,
   32.8593,
   33.5215,
   33.9174,
   34.0250,
   33.8807,
   33.5701,
   33.1802,
   32.7637,
   32.2714,
   31.6305,
   30.8181,
   29.9247,
   29.0551,
   28.2650,
   27.5183,
   26.7750,
   26.0506,
   25.3978,
   24.8311,
   24.3567,
   23.9744,
   23.6484,
   23.3344,
   22.9896,
   22.6543,
   22.4350,
   22.3966,
   22.4102,
   22.2976,
   22.0186,
   21.7907,
   21.6576,
   21.5586,
   21.4383,
   21.2787,
   21.0961,
   20.9084,
   20.7309,
   20.5773,
   20.4186,
   20.1984,
   19.8613,
   19.3997,
   18.9730,
   18.6786,
   18.4986,
   18.4115,
   18.3970,
   18.4311,
   18.4934,
   18.5608,
   18.6129,
   18.6274,
   18.5822,
   18.4550,
   18.2255,
   17.8696,])



output_dir = '/scratch/aa3878/2d_aneurysm/'
T = 2.85  #3.0            	# final time
num_steps = 10000*3 * 5 #6000   #60000   # number of time steps
dt = T / num_steps 	# time step size
mu = 0.04       		# dynamic viscosity
rho = 1.06           	# density

# Create mesh
tube = Rectangle(Point(0, 0), Point(3.8, 0.4))
aneurysm = Circle(Point(2.4, 0.7), 0.5)
domain = tube + aneurysm
mesh = generate_mesh(domain, 128)

# File(output_dir + 'mesh.pvd') << mesh

# Define function spaces
basis_order = 2 #Ali used 2
V = VectorFunctionSpace(mesh, 'P', basis_order)
Q = FunctionSpace(mesh, 'P', 1)

# Define boundaries
inflow   = 'near(x[0], 0)'
outflow  = 'near(x[0], 3.8)'
walls    = 'on_boundary && x[1]>=0.4 || near(x[1], 0)'

# class Inflow(SubDomain):
# 	def inside(self, x, on_boundary):
# 		return near(x[0], 0)
# class Outflow(SubDomain):
# 	def inside(self, x, on_boundary):
# 		return near(x[0], 3.8)
# class Walls(SubDomain):
# 	def inside(self, x, on_boundary):
# 		return on_boundary and x[1] >= 0.4 or near(x[1], 0)
# 
# inflow = Inflow()
# outflow = Outflow()
# walls = Walls()
# 				
# boundaries = MeshFunction('size_t', mesh, 1)
# boundaries.set_all(0)
# inflow.mark(boundaries, 1)
# outflow.mark(boundaries, 2)
# walls.mark(boundaries, 3)

# File(output_dir + 'boundaries.pvd') << boundaries					

# Define inflow profile
# inflow_profile = ('4.0*1.5*x[1]*(0.81 - x[1]) / pow(0.81, 2)', '0')

#inflow_profile = ('32', '0')

# Define boundary conditions
#bcu_inflow = DirichletBC(V, Expression(inflow_profile, degree = 2), inflow)
bcu_walls = DirichletBC(V, Constant((0., 0.)), walls)
#bcu = [bcu_inflow, bcu_walls]

bcp_outflow = DirichletBC(Q, Constant(0.), outflow)

bcp = [bcp_outflow]

# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = Function(V)
u_  = Function(V)
p_n = Function(Q)
p_  = Function(Q)

# Define expressions used in variational forms
U  = 0.5*(u_n + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0))
k  = Constant(dt)
mu = Constant(mu)
rho = Constant(rho)

# Define symmetric gradient
def epsilon(u):
	 return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
	 return 2*mu*epsilon(u) - p*Identity(len(u))

# Define variational problem for step 1
F1 = rho*dot((u - u_n) / k, v)*dx \
	+ rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
	+ inner(sigma(U, p_n), epsilon(v))*dx \
	+ dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
	- dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# Define variational problem for step 3
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
#[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# Create XDMF files for visualization output
xdmffile_u = XDMFFile(output_dir + 'velocity.xdmf')
xdmffile_p = XDMFFile(output_dir  + 'pressure.xdmf')

# Create time series (for use in reaction_system.py)
timeseries_u = TimeSeries(output_dir + 'velocity_series')
timeseries_p = TimeSeries(output_dir  + 'pressure_series')

# Save mesh to file (for use in reaction_system.py)
#File(output_dir + 'aneurysm.xml.gz') << mesh


# Create progress bar
progress = Progress('Looping', num_steps)
#set_log_level(LogLevel.PROGRESS)

# Time-stepping
t = 0
print('t', t)
for n in range(num_steps):

	v_inlet_BC = numpy.interp(t, time_input, V_vel)
	bcu_inflow = DirichletBC(V, Constant((v_inlet_BC, 0.)), inflow)
	bcu = [bcu_inflow, bcu_walls]

	# Apply boundary conditions to matrices
	[bc.apply(A1) for bc in bcu]

	# Step 1: Tentative velocity step
	b1 = assemble(L1)
	[bc.apply(b1) for bc in bcu]
	solve(A1, u_.vector(), b1, 'gmres', 'ilu')

	# Step 2: Pressure correction step
	b2 = assemble(L2)
	[bc.apply(b2) for bc in bcp]
	solve(A2, p_.vector(), b2, 'gmres', 'ilu')

	# Step 3: Velocity correction step
	b3 = assemble(L3)
	solve(A3, u_.vector(), b3, 'cg', 'sor')

	# Save solution to file (XDMF/HDF5)
	if n % 1000 == 0 and t > 1.9 :
		xdmffile_u.write(u_, t)
		# 	xdmffile_p.write(p_, t)
   		# Save nodal values to file
		#timeseries_u.store(u_.vector(), t)
    	#timeseries_p.store(p_.vector(), t)

	# Update previous solution
	u_n.assign(u_)
	p_n.assign(p_)

	# Update progress bar
	progress += 1

	# Update current time
	t += dt

	#print 'u max:', u_.vector().get_local().max()
	if n % 50 * 5 == 0:
	 print('u max:', u_.vector().get_local().max())
	 print('t', t)




