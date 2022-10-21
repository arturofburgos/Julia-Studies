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

T = 20 
dt = 0.0002 
N = Int(T/dt) + 1 # Number of time-steps
t = range(0, T, N)




# Matrix M̃

m1 = 1
m2 = 5 



# Matrix K̃



k1 = 8
k2 = 8



# State space representation u

x0_1 = 0.0 
x0_2 = 2.0


ẋ0_1 = 0.1
ẋ0_2 = 1.1


u1 = [x0_1 ẋ0_1]
u2 = [x0_2 ẋ0_2]


x1 = zeros(N)
x1[1] = x0_1

x2 = zeros(N)
x2[1] = x0_2


for i in 2:N
    
    u1[1] = u1[1] + dt * u1[2]
    u2[1] = u2[1] + dt * u2[2]

    u1[2] = u1[2] + dt * ((-k1/m1) * (u1[1]+1) + (k2/m1) * (u2[1]-u1[1]-2))
    u2[2] = u2[2] + dt * (-k2/m2) * (u2[1] - u1[1]-2)
    
    #Discuss this offset with Prof Andres.

    #= u1[2] = u1[2] + dt * ((-k1/m1) * (u1[1]+0) + (k2/m1) * (u2[1]-u1[1]-0))
    u2[2] = u2[2] + dt * (-k2/m2) * (u2[1] - u1[1]-0) =# # In this case I have the same results as the other schemes.
    
    
    x1[i] = u1[1]
    x2[i] = u2[1]
    
end



using Plots


plot(t, x1)
plot!(t,x2)