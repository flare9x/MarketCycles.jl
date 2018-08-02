VERSION >= v"0.0.1" && __precompile__(true)

module MarketCycles

export
    SuperSmoother, Decycler, Decycle_OSC, BandPassFilter, DominantCycle, HurstCoefficient, HPLPRoofingFilter,
    ZeroMeanRoofingFilterP0, ZeroMeanRoofingFilterP1, RoofingFilterIndicator,
    ModifiedStochastic, ModifiedRSI, AutoCorrelationIndicator, SingleLagAutoCorrelationIndicator, AutoCorrelationPeriodogram,
    AutoCorrelationReversals

    include("ehlers_cycles.jl")

end
