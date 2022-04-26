
module BlockAverage

using Statistics: mean

export block_average

"""

```
block_average(x::AbstractVector{T}, by = mean) where T<:Real
```

Function that computes the block average and standard deviation as a function of block size for a time-dependent, possibly correlated, set of data.

Returns three vectors `avg`, `err`, and `sizes`.  The output `avg` and `err` vectors will contain the estimate of the mean value of the property and the *standard deviation of the mean,* for each possible block size (stored in `sizes`). Thus, `avg[i]` for `sizes[i] = N` will contain the mean of `by(x)` computed for blocks of `x` of length `N` (thus, with `length(x)÷N` samples).

The function defined by `by=` must be able to operate on a `view` of a subvector of `x`.  By default, `x` is assumed to contain the value of the property to be considered, and then `by = mean`. Nevertheless, one could compute other property from `x`.

Call with complete set of parameters:

```julia
avg, err, sizes = block_average(x::AbstractVector{T}; by = mean, 
                                max_block_size::Int = length(x)÷10, 
                                min_block_size::Int = 1 ) where T<:Real
```

## Examples

With `x = rand(1000)`:

1. `x` contains directly the property to be estimated:
```julia
avg, err, sizes = block_average(x)
```

2. The property to be estimated is the average of the square of each `x`:
```julia
func(x) = sum(x.^2)/length(x)
avg, err, sizes = block_average(x,by=func)
```

3. The property is the probabilty of finding an element of `x` greater than `0.5`:
```julia
func(x) = count( x .> 0.5 ) / length(x)
avg, err, sizes = block_average(x,by=func)
```

"""
function block_average(
    x::AbstractVector{T};
    by = mean,
    min_block_size::Int = 1,
    max_block_size::Int = length(x) ÷ 10,
) where {T<:Real}

    n = length(x)

    xmean = T[]
    xerr = T[]
    sizes = Int[]

    for block_size in min_block_size:max_block_size
        nblocks = n ÷ block_size
        remaining = n % block_size

        # If the remaining is greater than the number of blocks, add
        # some points to each block to maximize the use of the data
        if remaining >= nblocks
            block_size += remaining ÷ nblocks
        end
        nblocks = n ÷ block_size

        # If the number of blocks didn't change, continue
        if length(sizes) > 0 && nblocks == n ÷ sizes[end]
            continue
        end

        # Add new point to vectors
        push!(sizes, block_size)
        push!(xmean, zero(T))
        push!(xerr, zero(T))

        # Compute the property in each block
        for i in 1:nblocks
            xblock = @view x[brange(i, block_size)]
            xmean[end] += by(xblock)
        end
        # Averaging over the number of blocks
        xmean[end] /= nblocks

        # Compute the standard deviation of the property estimate (σ²/√N)
        if nblocks > 1
            for i in 1:nblocks
                xblock = @view x[brange(i, block_size)]
                val = by(xblock)
                xerr[end] += (val - xmean[end])^2
            end
            xerr[end] = sqrt(xerr[end] / ((nblocks - 1) * nblocks))
        else
            xerr[end] = 0.
        end
    end

    return xmean, xerr, sizes
end

# Range of indices for block i of size block_size
function brange(i, block_size)
    ifirst = (i - 1) * block_size + 1
    ilast = ifirst + block_size - 1
    return ifirst:ilast
end

# Generate correlated data to test
function test_data(n)
    temperature = 1.0
    x = Vector{Float64}(undef, n)
    x[1] = 0.0
    u = 0.0
    i = 1
    while i < n
        x_trial = x[i] - 0.1 + 0.2 * rand()
        u_trial = x_trial^2
        if u_trial < u || exp((u - u_trial) / temperature) > 0.5
            i += 1
            x[i] = x_trial
            u = u_trial
        end
    end
    x
end

end # module
