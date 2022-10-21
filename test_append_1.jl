f = Array{Float64}(undef,0)
g = [1, 2, 3, 4, 5, 6]
h = [1,2,3]

for i in 1:3
    append!(f,[g[j] for j in h[i]])
end


println(f)