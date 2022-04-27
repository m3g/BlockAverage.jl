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

`x` is the original data set. 

`xmean` is the property value computed for all data (usually the mean, but not necessarily)

`blocksize` is an array of block sizes, in which the data was split. 
By default it goes from `1` to `length(x)`, with a number of points corresponding to the
number of integer divisions of `length(x)`.

`xmean_maxerr`: The property is computed for each block, and the maximum error (difference
between the property in the block and the average property) is stored in this array, for
each blocks size. 

`xmean_stderr`: The standard error of the estimates of the property, meaning the standar 
deviation of the estimates divided by the square root of the number of blocks. 

`autocor`: Is the autocorrelation function of the data, as a function of the lag. 

`lags`: Is the set of "time" lags for which the autocorrelation will be computed. Defined by
the `lags` parameter of the `block_average` function, as a range. 

`tau`: The characteristic decay time of the autocorrelation function, as obtained by fitting
of a single exponential, of the form `exp(-t/tau)` to the data. 

"""
struct BlockAverageData{T}
    x::Vector{T}
    xmean::T
    blocksize::Vector{Int}
    xmean_maxerr::Vector{T}
    xmean_stderr::Vector{T}
    lags::AbstractVector{Int}
    autocor::Vector{T}
    tau::T
end

function Base.show(io::IO,::MIME"text/plain",b::BlockAverageData)
    merr = findmax(b.xmean_stderr)
    print(io, """
    -------------------------------------------------------------------
    $(typeof(b))
    -------------------------------------------------------------------
    Estimated value (mean by default) = $(b.xmean)
    Length of data series: $(length(b.x))

    Block size ranges: $(extrema(b.blocksize))

    Maximum standard error (error, block size): $((merr[1], b.blocksize[merr[2]]))

    Deviations in last 3 blocks:
             percentual: $((100/b.xmean)*(b.xmean_maxerr[end-2:end] .- b.xmean))  
               absolute: $((b.xmean_maxerr[end-2:end] .- b.xmean))  

    Autocorrelation is first zero with lag: $(b.lags[findfirst(x -> x <=0, b.autocor)])
    Characteristic time of autocorrelation decay: 
            as fraction of series length: $(b.tau / length(b.x))
                                absolute: $(b.tau)
    -------------------------------------------------------------------
    """)
end

"""

```
block_average(
    x::AbstractVector{T};
    by = mean,
    min_block_size::Int = 1,
    max_block_size::Int = length(x),
    lags::Union{Nothing,AbstractVector{Int}} = nothing,
) where {T<:Real}
```

This function peforms some convergence analysis for a property computed from a series of data, typically a time-series. 
The data is given in vector `x`, and `by` defines the property to be estimated, typically, and by default, the mean value.

Two analyses are performed: a block averaging, in which the data is split in to blocks, and the mean value (or `by` value)
in each block is computed idependently. The output will contain the worst estimate obtained for all blocks, and the
standard error of the estimates, as a function of the block size. 

Finally, the autocorrelation function of the data is computed, and a single exponential is fitted, to obtain the
characteristic time of the decay of the correlation. 

The output will be a structure of type `BlockAverageData{T}`. See the corresponding help entry for more information.

All results can be plot with a convenience function `BlockAverage.plot`

The `lags` keyword can be tuned to define the range of intervals and length of the autocorrelation calculation, with
important implications to the exponential fit and correlation curve shape. See the `StatsBase.autocor` help for 
further information.

## Example

```julia-repl
julia> x = BlockAverage.test_data(10^6); # example data generator

julia> b = block_average(x, lags=0:100:10^5)
-------------------------------------------------------------------
BlockAverageData{Float64}
-------------------------------------------------------------------
Estimated value (mean by default) = -0.07941666750957592
Length of data series: 1000000

Block size ranges: (1, 1000000)

Maximum standard error (error, block size): (0.22783752512097308, 40000)

Deviations in last 3 blocks:
         percentual: [188.6496035297814, -22.9503589818173, -0.0]  
           absolute: [-0.1498192283933797, 0.01822641028484394, 0.0]  

Autocorrelation is first zero at step: 145
Characteristic time of autocorrelation decay: 
        as fraction of series length: 4.537259628690915e-5
                            absolute: 45.37259628690915
-------------------------------------------------------------------

julia> using Plots

julia> BlockAverage.plot(b) # creates a plot with the results

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

    xmean = by(x)
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
        lags =  0:round(Int,min(size(x,1)-1, 10*log10(size(x,1))))
        auto_cor = autocor(x)
    else
        auto_cor = autocor(x, lags)
    end

    tau = fitexp(lags, auto_cor, c=0., u=upper(a=1.1), l=lower(a=0.9)).b

    return BlockAverageData{T}(
        x,
        xmean,
        blocksize,
        xmean_maxerr,
        xmean_stderr,
        lags,
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

"""

```
plot(
    data::BlockAverageData; 
    xlims=nothing, 
    ylims=nothing,
    xscale=:identity,
    title="",
)
```

Function that creates a plot of output data from `block_average`. The optional 
`xlims`, `ylims`, `xscale`, `title`, can be set to adjust the aparency of the plot.    

The function requires `Plots` to be loaded first, as it is not an explict dependency
of the `BlockAverage` package.

## Example

```julia-repl
julia> x = BlockAverage.test_data(10^6);

julia> b = block_average(x, lags=0:100:10^5);

julia> using Plots

julia> BlockAverage.plot(b)
```

"""
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
    Plots.plot!(
        data.lags,
        data.autocor, 
        ylabel=L"c(\Delta t)",
        xlabel=L"\Delta t",
        label=nothing,
        linewidth=2,
        color=:black,
        subplot=3
    )
    Plots.plot!(
        data.lags,
        exp.(-inv(data.tau) .* data.lags),
        label=nothing,
        linewidth=2,
        color=:black,
        alpha=0.5,
        subplot=3,
    )
    Plots.annotate!(
        data.lags[end] - 0.2*data.lags[end],
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
