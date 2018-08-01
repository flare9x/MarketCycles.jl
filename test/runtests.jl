using MarketCycles
using CSV
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# count 0
count_nans(x) = sum(isnan.(x))

# Validate against John Ehlers Original TradeStation Easylanguage code

# SuperSmoother Test
test = CSV.read("C:/Users/Andrew.Bannerman/.julia/v0.6/MarketCycles/test/test_3-3_Supersmoother.csv")
close = Float64.(test[:Close])
super_smoother_benchmark = Float64.(test[:Ten_Period_Supersmoother])
Super_Smoother = SuperSmoother(close,n=10)
Super_Smoother = round.(Super_Smoother,2) # round same as tradestation output
valid = ifelse.(Super_Smoother .== super_smoother_benchmark,1,0)
valid = valid[28:length(valid)] # remove indicator lead in period
@test sum(valid) == length(valid)
