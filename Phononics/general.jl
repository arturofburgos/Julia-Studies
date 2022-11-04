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

T = 125
dt = 0.02 
# Here define the number of pair spring-mass: 
nm = 3
# Assuming all the masses are equal:
m = 3
# Assuming all springs are equal:
k = 8

# Initial Quantities
# Note that we should assign nm quatities for position and velocity

# Initial Position
#= x0 = [0.0 0.2 0.4] 
# Initial Velocity
ẋ0 = [0.1 0.5 1.1] =#

x0 = [0.1 0.0 0.0] 
# Initial Velocity
ẋ0 = [0.0 0 0]


#===========================#
# Procedure Numerical Model #
#===========================#

# Time grid:
N = Int(T/dt) + 1 # Number of time-steps
t = range(0, T, N)

# Define the zero matrix:
O = zeros(nm, nm)

# Define M̃ matrix
m_vec = fill(m, nm)
M̃ = diagm(m_vec)

# Define big M matrix 
M = [I O; O M̃]

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
k_vec = fill(k,nm)
K̃ = diagm(k_vec)

# Define K̂
K̂ = A'*K̃*A

# Define big K matrix
K = [O -I; K̂ O]

# Define big C matrix
C = (-inv(M) * K)

# Define State Space representation u:
u = [x0 ẋ0]'

# Defining solution array:
x = zeros(N,2*nm)
# Assingning the initial conditions to solution array
x[1,:] = [x0 ẋ0]

# External FORCE
f = (t -> sin.(t)*ones(nm))


#======================#
# Forward Euler Scheme #
#======================#

#= for i in 2:N
    u[:] = u[:] + dt *(C*u[:])
    x[i,:] = u[:]
end =#

function forward_euler(y, h, N, A, save_var, M, nm)
    for i in 2:N
        fn = f(t[i-1])
        fbig = [zeros(nm); fn]
        y[:] = y[:] + h * (A * y[:] + inv(M)*fbig)
        save_var[i,:] = y[:]
    end

    return save_var
end


#=======================#
# Backward Euler Scheme #
#=======================#

#= for i in 2:N     
    u[:] = inv(I - dt* C)*u[:] 
    x[i,:] = u[:]
end =#

function backward_euler(y, h, N, A, save_var, M, nm)
    for i in 2:N
        fn = f(t[i-1])
        fbig = [zeros(nm); fn]
        y[:] = inv(I - h * A) * y[:]
        save_var[i,:] = y[:]
    end

    return save_var
end


#====================#
# Trapezoidal Scheme #
#====================#

#= for i in 2:N 
    u[:] = inv(I - (dt/2)* C)*(I + (dt/2)* C)*u[:] 
    x[i,:] = u[:]
end =#

function trapezoidal(y, h, N, A, save_var)
    for i in 2:N
        y[:] = inv(I - (h/2) * A) * (I + (h/2) * A) * y[:]
        save_var[i,:] = y[:]
    end

    return save_var
end

x = forward_euler(u, dt, N, C, x, M, nm)
#x = backward_euler(u, dt, N, C, x)
#x = trapezoidal(u, dt, N, C, x)


#=================#
# Post Processing #
#=================#

plot(t,x[:,1:nm])