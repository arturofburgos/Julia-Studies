#========================================#
# Name: Arturo Burgos                    #
#                                        #
#                                        #
#           Spring Mass System           #
#                                        #
#========================================#

# Second order ODE: Harmonic Oscillation - two spring and two mass
using LinearAlgebra


# Time grid 

T = 15
dt = 0.02 
N = Int(T/dt) + 1 # Number of time-steps
t = range(0, T, N)


O = zeros(2,2)

# Matrix M̃

m1 = 3
m2 = 3 
m = [m1, m2]
M̃ = diagm(m)

M = [I O; O M̃]

# Matrix K̃

A = [1 0; -1 1]
At = A'

k1 = 8
k2 = 8
k = [k1, k2]
K̃ = diagm(k)

K̂ = At*K̃*A

K = [O -I; K̂ O] # CHECK

# State space representation u

x0_1 = 0.0 
x0_2 = 2.0
x0 = [x0_1 x0_2] # Initial Position

ẋ0_1 = 1.0
ẋ0_2 = 0.0
ẋ0 = [ẋ0_1 ẋ0_2]

u = [x0 ẋ0]'
un = zeros(4,1)

x = zeros(N,4)
x[1,:] = [x0 ẋ0]

C = (inv(M) * -K)



for i in 2:N 
    
    
    u[:] = u[:] + dt *(C*u[:])
    #println(u)
    x[i,:] = u[1:4]
    
end

using Plots

plot(t,x[:,1:2])

