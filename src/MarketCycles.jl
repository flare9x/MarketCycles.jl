VERSION >= v"0.0.2" && __precompile__(true)

module MarketCycles

export
    SuperSmoother, Decycler, Decycle_OSC, BandPassFilter, DominantCycle, HurstCoefficient, HPLPRoofingFilter,
    ZeroMeanRoofingFilterK0, ZeroMeanRoofingFilterK1, RoofingFilterIndicator,
    ModifiedStochastic, ModifiedRSI, AutoCorrelationIndicator, SingleLagAutoCorrelationIndicator, AutoCorrelationPeriodogram,
    AutoCorrelationReversals, DFTS, AdaptiveRSI, AdaptiveStochastic, PFish

    include("ehlers_cycles.jl")

end
