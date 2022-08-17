module Circle_mod

    mutable struct Circle
        radius
    end


    function circ(r::Real)
        Circle(r)
    end  

    export circ

end # end of module