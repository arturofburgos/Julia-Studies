include("a.jl")
include("b.jl")

using .Prism_mod, .Circle_mod

function main()

    p = Prisma()

    c = circ(5)

    println(c.radius)

    println("p.length = $(p.length)")

end

main()