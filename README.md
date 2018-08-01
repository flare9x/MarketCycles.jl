# MarketCycles

This package provides digital signal processing indicators developed by John F. Ehlers.

Currently the original indicators shared in his book: Cycle Analytics for Traders, Advanced Technical Trading Concepts are provided with intent to explore the DSP space and provide new intuitions based on the Ehlers framework. 

Outside contributions are welcome! 

## Available Indicators

*   Available Indicators
    *   Supersmoother
    *   Decycler
    *   Decycler Oscillator
    *   Band Pass Filter
    *   Hurst Coefficient 
    *   HP-LP Roofing Filter 
    *   Zero Mean Roofing Filter 
    *   Roofing Filter 
    *   Modified Stochastic 
    *   Modified RSI 
    *   Autocorrelation (Multiple Lag Matrix) 
    *   Autocorrelation (Single Lag) 
    *   Autocorrelation Periodogram 
    *   Autocorrelation Reversals 
    
*   TO DO
    *   Dominant Cycle
    *   DFT Spectral Estimate
    *   Comb Filter Spectral Estimate
    *   Adaptive RSI
    *   Adaptive Stochastic Indicator
    *   Adaptive CCI 
    *   Adaptive Band Pass Filter 
    *   Even Better SineWave Indicator 
    *   Compute and Display Convolution
    *   Classic Hilbert Transformer
    *   Hilbert Transformer Indicator
    *   Dominant Cycle Using the Dual Differentiator Method
    *   Dominant Cycle Using the Phase Accumulation Method
    *   Dominant Cycle Using the Homodyne Method
    *   Fisher Transform to the Adaptive RSI Indicator
    *   SwamiCharts RSI
    *   SwamiCharts Stochastic
    
 
## Usage
   ```julia
Pkg.add("MarketCycles.jl")
``` 

Each indicator function requires an input of a single dimension array of Float64 type. Each indicator function has its own set of arguments which change the level of some parameter. This may be the length of a look back period, the high or low pass periods or bandwidth value. For example we may call ```@doc function_name``` to see the associated documentation for the specific indicator function: 

   ```julia
julia> @doc AutoCorrelationReversals
  Autocorrelation Reversals - Equation 8-3

     The indicated reversals are very sensitive to the smoothing of the price data.
    ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

     LPLength is made available as an indicator input to decrease or increase the number of
 indicated reversals as desired.
    ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

     The AvgLength parameter is also made available as an indicator because this averaging also
 impacts the number of indicated reversals.
    ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

     Care should be taken when increasing the value of this input because the lag of the
 indicator increases in direct proportion to the increase of the value of the AvgLength.
    ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

     Typical delay of the indicator will be about three bars when the AvgLength parameter is set
 to a value of 3.
    ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

  AutoCorrelationReversals(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48,
  LPLength::Int64=10, HPLength::Int64=48, AvgLength::Int64=3)::Array{Float64}
```
Call the function as below for lags 1 to 48: 

 ```julia
AutoCorrelationReversals(your_data,min_lag=1,max_lag=48,LPLength=10,HPLength=48,AvgLength=3)
 ```
We may see what this looks like over dummy data: 

   ```julia
using MarketCycles
using Gadfly

# Generate dummy data
srand(121)
n = 1000
op = 100.0 + cumsum(randn(n))
hi = op + rand(n)
lo = op - rand(n)
cl = 100.0 + cumsum(randn(n))
index = collect(1:1:length(cl))
for i = 1:n
	if cl[i] > hi[i]
		cl[i] = hi[i]
	elseif cl[i] < lo[i]
		cl[i] = lo[i]
	end
end

# Apply autocorrelation reversals function
auto_cor_reversals = AutoCorrelationReversals(cl,min_lag=1,max_lag=48,LPLength=10,HPLength=48,AvgLength=3)

# Plot
white_panel = Theme(
    panel_fill="white",
    default_color="blue",
    background_color="white"
)
p1 = plot(x=index,y=cl,Geom.line,
Guide.xlabel(nothing), Guide.ylabel("Price"), Guide.title("Dummy Data"),white_panel)
p2 = plot(x=index,y=auto_cor_reversals,Geom.line,Guide.xlabel("Time Index"),Guide.title("Autocorrelation Reversals"), Guide.ylabel("Autocorrelation Reversals"),white_panel)
out = vstack(p1,p2)

# Save Plot
draw(PNG("C:/Users/Andrew.Bannerman/Desktop/Julia/auto_correlation_reversals.png", 1500px, 800px), out)
```

For the output: 

![John Ehlers Autocorrelation Reversals](https://github.com/flare9x/MarketCycles.jl/blob/master/examples/auto_correlation_reversals_readme.png)

One may line up the reversals to the turning points within the data!

## Available Function Names 
```julia
SuperSmoother, Decycler, Decycle_OSC, BandPassFilter, DominantCycle, HurstCoefficient, HPLPRoofingFilter,
    ZeroMeanRoofingFilterP0, ZeroMeanRoofingFilterP1, RoofingFilterIndicator,
    ModifiedStochastic, ModifiedRSI, AutoCorrelationIndicator, SingleLagAutoCorrelationIndicator, 
    AutoCorrelationPeriodogram, AutoCorrelationReversals
```
Feel free to explore any of the functions with:

```julia
@doc SuperSmoother
```


