#========================================#
# Name: Arturo Burgos                    #
#                                        #
#                                        #
#           Spring Mass System           #
#                                        #
#========================================#

# Second order ODE: Harmonic Oscillation - nm spring and nm mass


using LinearAlgebra, LinearMaps, Plots
#============#
# USER INPUT #
#============#

T = 12
dt = 0.005 
# Here define the number of pair spring-mass: 
ncells = 5 

nm = 2*ncells

# Assuming all the masses are equal:
m1 = 0.01

m2 = 0.02

m = [m1, m2]
# Assuming all springs are equal:
k1 = 20

k2 = 25 

k = [k1, k2]

# γ and β parameters:
γ = 1/2
β = 1/4



# Initial Quantities
# Note that we should assign nm quatities for position and velocity

# Initial Position
#x0 = [0.1 0.0 0.0] 
# Initial Velocity
#ẋ0 = [0.0 0 0]

# Initial Position
#x0 = [0.0 0.2 0.4] 
# Initial Velocity
#ẋ0 = [0.1 0.5 1.1]

#x0 = [0.0 0.2] 
# Initial Velocity
#ẋ0 = [0.1 0.5]

# Initial displacement
x0 = zeros(nm)
#x0[1] = 1.0
# Initial Velocity
ẋ0 = zeros(nm)


#===========================#
# Procedure Numerical Model #
#===========================#

# Time grid:
N = Int(T/dt) + 1 # Number of time-steps
t = range(0, T, N)



# Define M̃ matrix
m_vec = repeat(m, ncells)
M̃ = diagm(m_vec)


# Define A and At matrices: We either can create the actual A matrix
# Or we can apply a function Atimes! that maps a matrix. 

#========#
# Matrix #
#========#

A = zeros(nm,nm)

for i in 1:nm
    for j in 1:nm
        if j == i
            A[i,j] = 1
        end
        if j == i - 1
            A[i,j] = -1
        end
    end 
end

At = A'

#===========#
# LinearMap #
#===========#

#= """Implements the matrix multiplication y = Ax.

For example, if n = 3:

```
    | 1       |
y = |-1  1    | x
    |   -1  1 |
```
"""
function A_times!(y::AbstractVector, x::AbstractVector)

    n = length(x)
    @assert n == length(y) && n > 0

    y[1] = x[1]
    for i in 2:n
        y[i] = x[i] - x[i-1]
    end

    return y
end


"""Implements the matrix multiplication y = Ax.

For example, if n = 3:

```
    | 1 -1    |
y = |    1 -1 | x
    |       1 |
```
"""
function A_times_T!(y::AbstractVector, x::AbstractVector)

    n = length(x)
    @assert n == length(y) && n > 0

    y[n] = x[n]
    for i in 1:n-1
        y[i] = x[i] - x[i+1]
    end

    return y
end
 
A = LinearMap(A_times!,A_times_T!, nm) =#

# Define K̃ matrix
k_vec = repeat(k, ncells)
K̃ = diagm(k_vec)

# Define K̂
K̂ = A'*K̃*A

# Define big K matrix
#K = [O -I; K̂ O]

# Define big C matrix
#C = (-inv(M) * K)

#M̃ = [10 0 0; 0 20 0; 0 0 30]
#K̂ = 1e3*[45 -20 -15;-20 45 -25;-15 -25 40]

# Define P matrix
P = M̃/(β*dt^2) + K̂
P_inv = inv(P)

# Define u, ud, udd:
u = zeros(nm,N)
ud = zeros(nm,N)
udd = zeros(nm,N)

u[:,1] = x0
ud[:,1] = ẋ0
udd[:,1]= -inv(M̃)*(K̂*u[:,1])


# External force acting only (in my physical model) in the last mass
function f(t)
    a = zeros(nm)
    a[end] = 0.01
    return a
end

# Defining solution array:
#x = zeros(N,2*nm)
# Assingning the initial conditions to solution array
#x[1,:] = [x0 ẋ0]


#================#
# Newmark Scheme #
#================#
function newmark(u, ud, udd, N, β, γ, M̃, P)
    for i in 1:N-1 
        fn = f(t[i])
        Q = fn + M̃*((1/(β*dt^2))*u[:,i] + (1/(β*dt))*ud[:,i] + (1/(2*β) - 1)*udd[:,i])
        #= if i < 5 
            @show Q
        end =# # Just to check Q
        u[:, i+1] = inv(P)*Q 
        udd[:, i+1] = (1/(β*dt^2))*(u[:,i+1]-u[:,i]) - (1/(β*dt))*ud[:,i]-((1/(2*β))-1)*udd[:,i]
        ud[:, i+1] = ud[:,i] + (1-γ)*dt*udd[:,i] + γ*dt*udd[:, i+1]
    end

    return u
end


u = newmark(u, ud, udd, N, β, γ, M̃, P)

#=================#
# Post Processing #
#=================#

display(plot(t[:],u[end,1:end], ylim = (0, 0.01)))