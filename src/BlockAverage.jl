module BlockAverage

using DocStringExtensions
using LaTeXStrings
using Statistics: mean
using StatsBase
using EasyFit

export block_average, BlockAverageData

"""

$(TYPEDEF)

Structure that contains the result of the block-average analysis of the sequence. 

`xmean` is the property value computed for all data (usually the mean, but not necessarily)

`blocksize` is an array of block sizes, in which the data was split. 
By default it goes from `1` to `length(x)`, with a number of points corresponding to the
number of integer divisions of `lenth(x)`.

`xmean_maxerr`: The property is computed for each block, and the maximum error (difference
between the property in the block and the average property) is stored in this array.

`xmean_stderr`: The standard error of the estimates of the property, meaning the standar 
deviation of the estimates divided by the square root of the number of blocks. 

`autocor`: Is the autocorrelation function of the data, as a function of the lag. 

"""
struct BlockAverageData{T}
    xmean::T
    blocksize::Vector{Int}
    xmean_maxerr::Vector{T}
    xmean_stderr::Vector{T}
    autocor::Vector{T}
    tau::T
end

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
    max_block_size::Int = length(x),
    lags::Union{Nothing,AbstractVector{Int}} = nothing,
) where {T<:Real}

    n = length(x)

    xmean = mean(x)
    xmean_maxerr = T[]
    xmean_stderr = T[]
    blocksize = Int[]

    for block_size in min_block_size:max_block_size
        nblocks = n ÷ block_size
        remaining = n % block_size
        if remaining != 0
           continue
        end

        # Add new point to vectors
        push!(blocksize, block_size)
        push!(xmean_maxerr, zero(T))
        push!(xmean_stderr, zero(T))

        # Compute the property in each block, and keep the maximum error
        diff_max = -Inf 
        for i in 1:nblocks
            xblock = @view x[brange(i, block_size)]
            this_block_mean = by(xblock)
            diff = abs(this_block_mean - xmean)
            if diff > diff_max
                diff_max = diff
                xmean_maxerr[end] = this_block_mean
            end
            # Compute the standard deviation of the property estimate (σ²/√N)
            xmean_stderr[end] += (this_block_mean - xmean)^2
        end
        if nblocks > 1
            xmean_stderr[end] = sqrt(xmean_stderr[end] / (nblocks - 1))
            # We want the standard error, so
            xmean_stderr[end] /= sqrt(nblocks)
        end
    end

    # Compute auto-correlation function of the data
    if isnothing(lags)
        auto_cor = autocor(x)
    else
        auto_cor = autocor(x, lags)
    end

    tau = fitexp(1:length(auto_cor), auto_cor, c=0., u=upper(a=1.1), l=lower(a=0.9)).b

    return BlockAverageData{T}(
        xmean,
        blocksize,
        xmean_maxerr,
        xmean_stderr,
        auto_cor,
        tau 
    )

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

function plot(
    data::BlockAverageData; 
    xlims=nothing, 
    ylims=nothing,
    xscale=:identity,
    title="",
)
    Plots = Main.Plots
    p = Plots.plot(layout=(3,1))
    Plots.hline!(
        [data.xmean],
        linestyle=:dash,
        alpha=0.5,
        color=:black,
        label=:none,
        title=title,
        subplot=1,
    )
    Plots.plot!(
        data.blocksize, data.xmean_maxerr, 
        ylabel="worst block value",
        xlabel=L"\textrm{block~size~}(N)",
        label=nothing,
        linewidth=2,
        marker=:circle,
        color=:black,
        xscale=xscale,
        subplot=1
    )
    Plots.annotate!(
        maximum(data.blocksize) - 0.1*maximum(data.blocksize),
        maximum(data.xmean_maxerr) - 0.1*(maximum(data.xmean_maxerr)-minimum(data.xmean_maxerr)),
        Plots.text("mean = $(round(data.xmean, digits=2))", "Computer Modern", 12, :right),
        subplot=1,
    )
    Plots.plot!(data.blocksize, data.xmean_stderr, 
        ylabel=L"\sigma^2 / \sqrt{N}",
        xlabel=L"\textrm{block~size~}(N)",
        label=nothing,
        linewidth=2,
        marker=:circle,
        color=:black,
        xscale=xscale,
        subplot=2
    )
    # Auto correlation function
    Plots.plot!(data.autocor, 
        ylabel=L"c(\Delta t)",
        xlabel=L"\Delta t",
        label=nothing,
        linewidth=2,
        color=:black,
        subplot=3
    )
    Plots.plot!(
        exp.(-inv(data.tau) * (0:length(data.autocor))),
        label=nothing,
        linewidth=2,
        color=:black,
        alpha=0.5,
        subplot=3,
    )
    Plots.annotate!(
        length(data.autocor) - 0.2*length(data.autocor),
        0.80,
        Plots.text("τ = $(round(data.tau, digits=2))", "Computer Modern", 12, :right),
        subplot=3,
    )
    Main.plot!(
        p,
        size=(400,600),
        framestyle=:box,
        fontfamily="Computer Modern",
        xlims = xlims,
        ylims = ylims,
        leftmargin=0.5Plots.Measures.cm,
        rightmargin=0.5Plots.Measures.cm,
    )
    return p
end

end # module
