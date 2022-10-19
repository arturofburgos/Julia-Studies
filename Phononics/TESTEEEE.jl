#========================================#
# Name: Arturo Burgos                    #
#                                        #
#                                        #
#           Spring Mass System           #
#                                        #
#========================================#

# Second order ODE: Harmonic Oscillation - two spring and two mass
using LinearAlgebra


# TIME GRID

T = 2
dt = 0.0002 
N = Int(T/dt) + 1 # Number of time-steps
t = range(0, T, N)


# MATRIX O

O = zeros(2,2)

# MATRIX A AND At

A = [1 0 ; -1 1]
At = A'


# MATRIX K
k1 = 8
k2 = 8
k = [k1, k2]
K = diagm(k)

# MATRIX K̂

K̂ = At*K*A

# MATRIX M

m1 = 3
m2 = 3
m = [m1, m2]
M = diagm(m)


M̃ = [I O ; O M]
K̃ = [O I ; K̂ O]

x0_1 = 0.0 # Initial Position
x0_2 = 2.0
x0 = [x0_1 x0_2]

ẋ0_1 = 1.0 # Initial Velocity
ẋ0_2 = 1.0
ẋ0 = [ẋ0_1 ẋ0_2]

u = [x0 ẋ0]'

x = zeros(N,4)
x[1,:] = [x0 ẋ0]
################################################
# Inverse is: 

#= A^-1

# Division is multiplying by inverse:

A/K
K^-1*A
A/K == A*(K^-1) =#

##############################################
c = inv(M̃) * K̃
c[1:2, 3:4] .*= -1

println(u)
for i in 2:N 
    
    
    u[:] = u[:] + dt *(-c*u[:])
    println(u)
    x[i,:] = u[1:4]
    
end

using Plots

plot(t,x)



#, xlim = (0,0.5),ylim = (-5,5)