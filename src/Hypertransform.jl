module Hypertransform

export hypertriangulate, hypercubify

"""
    hypertriangulate(x::AbstractArray; bounds::Tuple{<:Real, <:Real} = (0.0, 1.0))

Maps a vector or matrix from hypercube space to hypertriangle space, scaling each element to the specified bounds.

The hypercube is the space the samplers usually work in; the 
    components of x are in no particular order.
    
    The hypertriangle is the space where the components are sorted into
    ascenting order, \( y0 < y1 < ... < yn \). 
    
    The (unit) transformation is defined by:

    .. math::
        y_j = 1 - \\prod_{i=0}^{j} (1 - x_i)^{1/(n-i)}

    Example application. If we are analysing a number num_dim of DWD 
    sources, all with identical priors. Then this function would be
    called on the array `np.array([f_1, f_2, ..., f_num_sources])` with
    `bounds=(f_min, f_max)`.
    
# Arguments
- `x`: An array, either a vector of length `n` or a matrix of size `(m, n)`, where each row represents a point in hypercube space.
- `bounds`: A tuple `(lower, upper)` specifying the minimum and maximum values for the output range. Defaults to `(0.0, 1.0)`.

# Returns
- An array of the same shape as `x`, with each element transformed to lie within the specified bounds in hypertriangle space.

# Example
```julia
using Hypertransform
x = [0.5, 0.75]
bounds = (0.0, 1.0)
y = hypertriangulate(x, bounds=bounds)
```
# Raises
- `ArgumentError`: If any element of `x` is outside the specified bounds.


"""
function hypertriangulate(x::AbstractArray; bounds::Tuple{<:Real, <:Real} = (0.0, 1.0))
    lower, upper = bounds
    # Explicit bounds check
    if any(x .< lower) || any(x .> upper)
        throw(ArgumentError("Values outside bounds ($lower, $upper) passed to hypertriangulate"))
    end

    unit_x = (x .- lower) ./ (upper - lower)

    is_vec = ndims(unit_x) == 1
    unit_x = is_vec ? reshape(unit_x, 1, :) : unit_x
    
    n = size(unit_x, 2)
    m = size(unit_x, 1)
    idx = collect(0:(n - 1))
    exponents = 1.0 ./ (n .- idx)

    # Apply to each row 
    inner_term = (1 .- unit_x) .^ reshape(exponents, 1, :)
    result = similar(inner_term)
    for i in 1:m
        result[i, :] = cumprod(inner_term[i, :])
    end
    unit_y = 1 .- result
    y = map(y -> y .* (upper-lower) .+ lower, unit_y)
    return is_vec ? vec(y[1, :]) : reshape(y, size(unit_x))
end

"""
    hypercubify(y::AbstractArray; bounds::Tuple{<:Real, <:Real} = (0.0, 1.0))

Transform a vector or matrix from the hypertriangle space back to the hypercube.

# Arguments
- `y`: A vector of length `n` or a matrix of size `(m, n)` where each row is a point in the hypertriangle space.
- `bounds`: A tuple `(lower, upper)` defining the target parameter bounds. Defaults to `(0.0, 1.0)`.

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
