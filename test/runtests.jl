using MarketCycles
using CSV
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

######################################################################
# Validate against John Ehlers Original TradeStation Easylanguage code
######################################################################

# 3-3 SuperSmoother Test
@testset "Ehlers Tests" begin
test = CSV.read("C:/Users/Andrew.Bannerman/.julia/v0.6/MarketCycles/test/test_3-3_Supersmoother.csv")
close = Float64.(test[:Close])
super_smoother_benchmark = Float64.(test[:Ten_Period_Supersmoother])
Super_Smoother = SuperSmoother(close,n=10)
Super_Smoother = round.(Super_Smoother,2) # round same as tradestation output
valid = ifelse.(Super_Smoother .== super_smoother_benchmark,1,0)
valid = valid[28:length(valid)] # remove indicator lead in period
@test sum(valid) == length(valid)

# 4-1 Decycler Test
test = CSV.read("C:/Users/Andrew.Bannerman/.julia/v0.6/MarketCycles/test/test_4-1_Decycler.csv")
decycler_benchmark = Float64.(test[:Sixty_Period_Decycler])
decycler = Decycler(close, n=60)
decycler = round.(decycler,2) # round same as tradestation output
valid = ifelse.(decycler .== decycler_benchmark,1,0)
@test sum(valid) == length(valid)
end

# 4-2 - Decycle Oscillator test
test = CSV.read("C:/Users/Andrew.Bannerman/.julia/v0.6/MarketCycles/test/test_4-2_Decycle_Oscillator.csv")
decycler_osc_benchmark = Float64.(test[:thirty_sixty_decycle_osc])
decycler_osc = Decycle_OSC(close, n1=30,n2=60)
decycler_osc = round.(decycler,2) # round same as tradestation output
valid = ifelse.(decycler_osc .== decycler_osc_benchmark,1,0)
@test sum(valid) == length(valid)
