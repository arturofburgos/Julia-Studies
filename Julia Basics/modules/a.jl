module Prism_mod

    struct Prisma
        length::Float64
        width::Float64
        height::Float64
        
        function Prisma()
            new(1,1,1) # This new() function is used to construct a new struct 
        end
        
        function Prisma(l::Float64, w::Float64, h::Float64)
            if l < 0 || w < 0 || h < 0
                error("Cannot have negative number!")
            elseif w < l
                error("Cannot have width shorter than length!")
            else
                new(l,w,h)
            end
        end
    end

    export Prisma

end # end of module