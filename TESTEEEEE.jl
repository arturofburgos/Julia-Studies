f = t-> sin.(t)*ones(3)

T = 25
dt = 0.02 
N = Int(T/dt) + 1 # Number of time-steps
t = range(0, T, N)

for i in 1:5
    global fn = f(t[i])
    println(fn)
end
