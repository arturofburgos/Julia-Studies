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
dt = 0.002 
N = Int(T/dt) + 1 # Number of time-steps
t = range(0, T, N)


# External FORCE
f = (t -> sin(t)*ones(2))

# Matrix O
O = zeros(2,2)

A = [1 0 ; -1 1]
At = A'
k1 = 8
k2 = 8
k = [k1, k2]
K = diagm(k)

m1 = 3
m2 = 3
m = [m1, m2]
M = diagm(m)


M̃ = [I O ; O M]
K̃ = [O I ; K O]

x0_1 = 0.0 # Initial Position
x0_2 = 2.0
x0 = [x0_1 x0_2]

ẋ0_1 = 1.0 # Initial Velocity
ẋ0_2 = 0.5
ẋ0 = [ẋ0_1 ẋ0_2]

u = [x0 ẋ0]'
ut = [0 0 ẋ0]'
x = zeros(N,2)
x[1,:] = x0 
################################################
# Inverse is: 

A^-1

# Division is multiplying by inverse:

A/K
A/K == A*(K^-1)

##############################################
c = K̃/M̃


for i in 2:N 
    fn = f(t[i-1])
    fbig = [zeros(2); fn]
    u[:] = u[:] + dt*(-c*u[:] + M̃^-1*fbig)
    #= u[1:2]  = u[1:2] + dt * u[3:4]
    u[3:4] = u[3:4] + dt * c*u[3:4] =#
    x[i,:] = u[1:2]
end

using Plots

plot(t,x[:,1])
plot!(t,x[:,2])