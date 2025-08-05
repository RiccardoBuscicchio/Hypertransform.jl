module Hypertransform

export hypertriangulate, hypercubify

"""
    hypertriangulate(x::AbstractArray; bounds::Tuple{<:Real, <:Real} = (0.0, 1.0))
    Transform a vector or matrix from the hypercube space to the hypertriangle space.
    Each input vector is mapped such that the output is in the range defined by `bounds`.
    # Arguments
    - `x`: A vector of length `n` or a matrix of size `(m, n)` where each row is a point in the hypercube space.
    - `bounds`: A tuple `(lower, upper)` defining the target parameter bounds.
    # Returns
    - Transformed vector or matrix of the same shape as input, mapped to the hypertriangle space.
"""
function hypertriangulate(x::AbstractArray; bounds::Tuple{<:Real, <:Real} = (0.0, 1.0))
    lower, upper = bounds
    unit_x = (x .- lower) ./ (upper - lower)

    is_vec = ndims(unit_x) == 1
    unit_x = is_vec ? reshape(unit_x, 1, :) : unit_x
    
    n = size(unit_x, 2)
    m = size(unit_x, 1)
    idx = collect(0:(n - 1))
    exponents = 1.0 ./ (n .- idx)

    try
        # Apply to each row 
        inner_term = (1 .- unit_x) .^ reshape(exponents, 1, :)
        result = similar(inner_term)
        for i in 1:m
            result[i, :] = cumprod(inner_term[i, :])
        end
        unit_y = 1 .- result
        y = map(y -> y .* (upper-lower) .+ lower, unit_y)
        return is_vec ? vec(y[1, :]) : reshape(y, size(unit_x))
    catch e
        throw(ArgumentError("Values outside bounds passed to hypertriangulate: $e"))
    end
end

"""
    hypercubify(y::AbstractArray; bounds::Tuple{<:Real, <:Real} = (0.0, 1.0))
    Transform a vector or matrix from the hypertriangle space back to the hypercube.
    Each input vector is mapped such that the output is in the range defined by `bounds`.
    # Arguments
    - `y`: A vector of length `n` or a matrix of size `(m, n)` where each row is a point in the hypertriangle space.
    - `bounds`: A tuple `(lower, upper)` defining the target parameter bounds.
    # Returns
    - Transformed vector or matrix of the same shape as input, mapped back to the hypercube space.
"""
function hypercubify(y::AbstractArray; bounds::Tuple{<:Real, <:Real} = (0.0, 1.0))    
    lower, upper = bounds
    unit_y = (y .- lower) ./ (upper - lower)
    is_vec = ndims(unit_y) == 1
    unit_y = is_vec ? reshape(unit_y, 1, :) : unit_y
    n = size(unit_y, 2)
    idx = collect(0:(n - 1))
    powers = Float64.(n .- idx)
    unit_y_shifted = hcat(zeros(size(unit_y, 1)), unit_y[:, 1:end-1])
    frac = (1 .- unit_y) ./ (1 .- unit_y_shifted)
    #frac = [(1 .- row) ./ (1 .- row_shifted) for (row, row_shifted) in zip(eachrow(unit_y), eachrow(unit_y_shifted))]
    unit_x = 1 .- (frac .^ powers') #[1 .- (f .^ powers) for f in eachrow(frac)]
    x = lower .+ unit_x .* (upper-lower)
    #x = map(x -> x .* (upper-lower) .+ lower, unit_x)
    return is_vec ? vec(x) : reshape(x, size(unit_y))
    end

end # module
