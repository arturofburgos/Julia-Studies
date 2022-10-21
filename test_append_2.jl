f = Array{Float64}(undef,0)
g = [1, 2, 3, 4, 5, 6]
# = [1,2,3]

for i in 1:3
    append!(f,fill(i))
end


println(f)