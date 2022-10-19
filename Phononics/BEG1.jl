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

T = 1 
dt = 0.0002 
N = Int(T/dt) + 1 # Number of time-steps
t = range(0, T, N)


O = zeros(2,2)

# Matrix M̃

m1 = 3
m2 = 3 
m = [m1, m2]
M = diagm(m)

M̃ = [I O; O M]

# Matrix K̃

A = [1 0; -1 1]
At = A'

k1 = 8
k2 = 8
k = [k1, k2]
K = diagm(k)

K̂ = At*K*A

K̃ = [O I; K̂ O]

# State space representation u

x0_1 = 0.0 
x0_2 = 2.0
x0 = [x0_1 x0_2] # Initial Position

ẋ0_1 = 1.0
ẋ0_2 = 1.0
ẋ0 = [ẋ0_1 ẋ0_2]

u = [x0 ẋ0]'
un = zeros(4,1)

x = zeros(N,2)
x[1,:] = x0

C = inv(M̃) * K̃
println("for 1")
println(un)
println(u)
println(x[1:5,:])
println("\n")

for i in 2:N
    println("for $i")
    
    u[:] = u[:] + dt *(-C*un[:])
    println(un)
    x[i,:] = un[1:2]
    println(x[1:5,:])
    u[:] = un[:]
    println(u)
    println("\n")
    
end



using Plots

plot(t,x[:,1])
plot!(t,x[:,2])
