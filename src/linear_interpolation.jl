using LinearAlgebra, SpecialFunctions, Interpolations

export interp1_linear, linear_interpolation_1D


"""
    linear_interpolation_1D(x, y, x_interp)

Perform linear interpolation in 1D.

# Arguments
- `x::AbstractArray`: The x coordinates of data points.
- `y::AbstractArray`: The y coordinates of data points.
- `x_interp::Real`: The x coordinate at which interpolation is desired.

# Returns
- `y_interp::Real`: The interpolated y value at `x_interp`.

"""
function linear_interpolation_1D(x, y, x_interp)
    if length(x) != length(y)
        #Check if input arrays are of the same length
        error("x data and y data must be of the same length. Please check your inputs.")
    elseif x_interp < minimum(x) || x_interp > maximum(x)
        #Check if x_interp is within the bounds of x
        error("x_interp is out of bounds of x.")
    end
    # Find the interval in which x_interp lies and interpolate
    for (i,_) in enumerate(1:length(x)-1)
        if x_interp >= x[i] && x_interp <= x[i + 1]
            # Linear interpolation formula for a unequally spaced grid
            x0, y0   = x[i], y[i]
            x1, y1   = x[i + 1], y[i + 1]
            N1       = 1 - (x_interp - x0) * inv(x1 - x0)
            N2       =     (x_interp - x0) * inv(x1 - x0)
            y_interp =  N1 * y0 + N2 * y1
            return y_interp
        end
    end
end
"""
    interp1_linear(x, y, xq, extrapval)

Performs one-dimensional linear interpolation.

# Arguments
- `x::AbstractVector`: Vector of known x-coordinates (must be sorted in ascending order).
- `y::AbstractVector`: Vector of known y-values corresponding to `x`.
- `xq::AbstractVector`: Vector of query points where interpolation is desired.
- `extrapval`: Value to use for points in `xq` that are outside the range of `x`.

# Returns
- `yq::Vector`: Interpolated values at the query points `xq`. For points outside the range of `x`, returns `extrapval`.

# Notes
- Assumes `x` is strictly increasing.
- If `xq` contains values outside the range of `x`, those values are assigned `extrapval`.

"""

function interp1_linear(x, y, xq, extrapval)
    itp = LinearInterpolation(x, y, extrapolation_bc=Flat())
    v = itp(xq)
    if xq < minimum(x) || xq > maximum(x)
        return extrapval
    else
        return v
    end
end