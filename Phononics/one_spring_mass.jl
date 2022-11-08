#========================================#
# Name: Arturo Burgos                    #
#                                        #
#                                        #
#           Spring Mass System           #
#                                        #
#========================================#

# Second order ODE: Harmonic Oscillation - one spring and one mass

# Create the time grid: 

T = 20 # Lenght of the time interval
dt = 0.002 # Time-step size
N = Int(T/dt) + 1 # Number of time-steps
t = range(0, T, N)

# Constants: 

k = 8
m = 3

c = sqrt(k/m)

# Using Eulers Method: 

# Set initial conditions:

x0 = 0.0 # Initial Position 
ẋ0 = 1.0 # Initial Velocity 

u = [x0 ẋ0] # Initial condition matrix 

# Array to store the position at each time-step: 

x = zeros(N)
x[1] = x0

for i in 2:N # Option 1
    rhs = [u[2], (-c^2)*u[1]]
    u[:] = u[:] + dt * rhs
    x[i] = u[1]
end

using Plots, LaTeXStrings, Measures # Possible to use PyPlot too.
gr()


@userplot springmass
@recipe function f!(var::springmass)
        
    size --> (950,400)
    margin --> 5mm
    top_margin --> 3mm
    framestyle --> :box
    grid --> :true
    gridalpha --> 5
    gridwidth --> 0.3
    minorgrid --> :true
    minorgridalpha --> 1
    minorgridwidth --> 0.05
    fontfamily --> "Computer Modern"

    t, x = var.args
    title --> "Mass position over time"
    legend --> :topright
    
    
    xaxis --> ("Time [s]")
    yaxis -->("Position [m]")
    #yguidefontrotation --> -90
    seriestype --> :line 
    linewidth --> 1.5
    #color --> :red 
    
    return t, x

end


figure = springmass(t,x, label = "Numerical")
display(figure)

# Compare with the Analytical Solution: 

x_exact = x0 * cos.(c*t) + (ẋ0/c) * sin.(c*t)
springmass!(t,x_exact, label = "Analytical")

# Check Error and Convergence

dt_values = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001] # Time step sizes in order to check convergence and error

# Same as below, but less efficient.
#= for dt_sizes in 1:size(dt_values,1)
 
    N = Int(T/dt_values[dt_sizes]) + 1 # Number of time-steps
    println(N)

end =#


for i in eachindex(dt_values)
 
    N = Int(T/dt_values[i]) + 1 # Number of time-steps
    t = range(0.0, T, length = N) # Time grid
    # Set the initial conditions: 
    z = [x0, ẋ0]
    X = []
    push!(X,x0)

    for n in 2:N
        push!(X,0.0)
        rhs = [z[2], (-c^2)*z[1]]
        z[:] = z[:] + dt_values[i] * rhs
        X[i] = z[1] 
    end

    


end