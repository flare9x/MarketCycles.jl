VERSION >= v"0.0.1" && __precompile__(true)

module MarketCycles

export
    SuperSmoother, Decycler, Decycle_OSC, BandPassFilter, DominantCycle, HurstCoefficient, HPLPRoofingFilter,
    ZeroMeanRoofingFilterK0, ZeroMeanRoofingFilterK1, RoofingFilterIndicator,
    ModifiedStochastic, ModifiedRSI, AutoCorrelationIndicator, SingleLagAutoCorrelationIndicator, AutoCorrelationPeriodogram,
    AutoCorrelationReversals, DFT, AdaptiveRSI, AdaptiveStochastic

    include("ehlers_cycles.jl")

end
