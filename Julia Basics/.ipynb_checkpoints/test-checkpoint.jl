


mutable struct Student{T}
    age::Int64
    label::T
    function Student(x::Int64, label::Any)
        new{typeof(label)}(x,label)
    end
    function Student()
        Student(0, 0)
    end
end