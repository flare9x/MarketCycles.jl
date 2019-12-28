using MarketCycles
using CSV
using Statistics

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

######################################################################
# Validate against John Ehlers Original TradeStation Easylanguage code
# Run against dummy data loaded in Tradestation
######################################################################

# Note - for some tests irratic results for the lead in period #

# 3-3 SuperSmoother Test
@testset "Ehlers" begin
    @testset "Super Smoother - Equation 3-3" begin
        filename = joinpath(dirname(@__FILE__), "test_3-3_Supersmoother.csv")
        test = CSV.read( filename)
        dat = Float64.(test[:, :x])
        super_smoother_benchmark = Float64.(test[:Ten_Period_Supersmoother])
        Super_Smoother = SuperSmoother(dat,n=10)
        Super_Smoother = round.(Super_Smoother; digits=2) # round same as tradestation output
        valid = ifelse.(Super_Smoother .== super_smoother_benchmark,1,0)
        valid = valid[28:length(valid)] # remove indicator lead in period
        @test sum(valid) == length(valid)
        # Passed but takes 28 lead in bars to do so
    end

    @testset "Decycler - Equation 4-1" begin
        filename = joinpath(dirname(@__FILE__), "test_4-1_Decycler.csv")
        test = CSV.read( filename)
        dat = Float64.(test[:, :x])
        decycler_benchmark = Float64.(test[:Sixty_Period_Decycler])
        decycler = Decycler(dat, n=60)
        decycler = round.(decycler; digits=2) # round same as tradestation output
        valid = ifelse.(decycler .== decycler_benchmark,1,0)
        @test sum(valid) == length(valid)
    end

    @testset "Decycle Oscillator - Equation 4-2" begin
        filename = joinpath(dirname(@__FILE__), "test_4-2_Decycle_Oscillator.csv")
        test = CSV.read( filename)
        dat = Float64.(test[:, :x])
        decycler_osc_benchmark = Float64.(test[:thirty_sixty_decycle_osc])
        decycler_osc = Decycle_OSC(dat, LPPeriod=30,HPPeriod=60)
        decycler_osc = round.(decycler_osc; digits=2) # round same as tradestation output
        valid = ifelse.(decycler_osc .== decycler_osc_benchmark,1,0)
        @test sum(valid) == length(valid)-148 # minus lead in
    end

    @testset "Band Pass Filter - Equation 5-1" begin
        filename = joinpath(dirname(@__FILE__), "test_5-1_Band_Pass_Filter.csv")
        test = CSV.read( filename)
        dat = Float64.(test[:, :x])
        BP_filt_benchmark = Float64.(test[:Trigger])
        BP_filt = BandPassFilter(dat,n=20,bandwidth=.3)
        BP_filt = round.(BP_filt;digits=2) # round same as tradestation output
        valid = ifelse.(BP_filt .== BP_filt_benchmark,1,0)
        @test sum(valid) == length(valid)-92  # 92 bar lead in period
    end

    @testset "Hurst Coefficient - Equation 6-1" begin
        filename = joinpath(dirname(@__FILE__), "test_6-1_Hurst_Coefficient.csv")
        test = CSV.read( filename)
        dat = Float64.(test[:, :x])
        Hurst_benchmark = Float64.(test[:Hurst_Coefficient])
        Hurst = HurstCoefficient(dat,n=30)
        Hurst = round.(Hurst;digits=4) # round same as tradestation output
        valid = ifelse.(Hurst .== Hurst_benchmark,1,0)
        @test sum(valid) == length(valid)-37 # 37 bar lead in period
    end

    @testset "HP LP Roofing Filter - Equation 7-1" begin
        filename = joinpath(dirname(@__FILE__), "test_7-1_HP_LP_Roofing_Filter.csv")
        test = CSV.read( filename)
        dat = Float64.(test[:, :x])
        HP_LP_Roof_benchmark = Float64.(test[:HP_LP_Roofing_Filter])
        HP_LP_Roof = HPLPRoofingFilter(dat)
        HP_LP_Roof = round.(HP_LP_Roof;digits=2) # round same as tradestation output
        valid = ifelse.(HP_LP_Roof .== HP_LP_Roof_benchmark,1,0)
        @test sum(valid) == length(valid)-48 # 48 bar lead in period
    end

    @testset "Zero Mean Roofing Filter - Lag 0 - Equation 7-2" begin
        filename = joinpath(dirname(@__FILE__), "test_7-2_Zero_Mean_Roofing_Filter.csv")
        test = CSV.read( filename)
        dat = Float64.(test[:, :x])
        Zero_Mean_Roofing_Filter_lag_0_benchmark = Float64.(test[:Filt])
        Zero_Mean_Roofing_Filter_lag_0 = HPLPRoofingFilter(dat)
        Zero_Mean_Roofing_Filter_lag_0 = round.(Zero_Mean_Roofing_Filter_lag_0;digits=2) # round same as tradestation output
        valid = ifelse.(Zero_Mean_Roofing_Filter_lag_0 .== Zero_Mean_Roofing_Filter_lag_0_benchmark,1,0)
        @test sum(valid) == length(valid)-48 # 48 bar lead in period
    end

    @testset "Zero Mean Roofing Filter - Lag 1 - Equation 7-2" begin
        filename = joinpath(dirname(@__FILE__), "test_7-2_Zero_Mean_Roofing_Filter.csv")
        test = CSV.read( filename)
        dat = Float64.(test[:, :x])
        Zero_Mean_Roofing_Filter_lag_1_benchmark = Float64.(test[:Filt2])
        Zero_Mean_Roofing_Filter_lag_1 = ZeroMeanRoofingFilterK1(dat)
        Zero_Mean_Roofing_Filter_lag_1 = round.(Zero_Mean_Roofing_Filter_lag_1;digits=2) # round same as tradestation output
        valid = ifelse.(Zero_Mean_Roofing_Filter_lag_1 .== Zero_Mean_Roofing_Filter_lag_1_benchmark,1,0)
        @test sum(valid) == length(valid)-60 # 48 bar lead in period
    end

    @testset "Roofing Filter Indicator Equation 7-3" begin
        filename = joinpath(dirname(@__FILE__), "test_7-3_Roofing_Filter_Indicator.csv")
        test = CSV.read(filename)
        dat = Float64.(test[:, :x])
        Roofing_Filter_Indicator_benchmark = Float64.(test[:Filt])
        Roofing_Filter_Indicator = RoofingFilterIndicator(dat)
        Roofing_Filter_Indicator = round.(Roofing_Filter_Indicator;digits=2) # round same as tradestation output
        valid = ifelse.(Roofing_Filter_Indicator .== Roofing_Filter_Indicator_benchmark,1,0)
        @test sum(valid) == length(valid)-213 # 213 bar lead in period
    end

    @testset "Modified Stochastic - Equation 7-4" begin
        filename = joinpath(dirname(@__FILE__), "test_7-4_Modified_Stochastic.csv")
        test = CSV.read(filename)
        dat = Float64.(test[:, :x])
        Modified_stochastic_benchmark = Float64.(test[:Modified_stochasitc])
        Modified_stochastic = ModifiedStochastic(dat,n=20)
        Modified_stochastic = round.(Modified_stochastic; digits=2) # round same as tradestation output
        valid = ifelse.(Modified_stochastic .== Modified_stochastic_benchmark,1,0)
        @test sum(valid) == length(valid)-72 # 72 bar lead in period
    end

    @testset "Modified RSI - Equation 7-5" begin
        filename = joinpath(dirname(@__FILE__), "test_7-5_Modified_RSI.csv")
        test = CSV.read(filename)
        dat = Float64.(test[:, :x])
        Modified_RSI_benchmark = Float64.(test[:mod_rsi])
        Modified_RSI = ModifiedRSI(dat,n=10)
        Modified_RSI = round.(Modified_RSI;digits=2) # round same as tradestation output
        valid = ifelse.(Modified_RSI .== Modified_RSI_benchmark,1,0)
        @test sum(valid) == length(valid)-96 # 96 bar lead in period
    end

    @testset "Autocorrelation Indicator - Equation 8-2" begin
        filename = joinpath(dirname(@__FILE__), "test_8-3_Autocorrelation_Indicator.csv")
        test = CSV.read(filename)
        dat = Float64.(test[:, :x])
        Autocor_p48_benchmark = Float64.(test[:cor_48])
        Autocor_p10_benchmark = Float64.(test[:cor_10])
        Autocor_ind_out = AutoCorrelationIndicator(dat)
        autocorr_p48 = Autocor_ind_out[:,48]
        autocorr_p10 = Autocor_ind_out[:,10]
        autocorr_p48 = round.(autocorr_p48;digits=4) # round same as tradestation output
        autocorr_p10 = round.(autocorr_p10;digits=4) # round same as tradestation output
        valid = ifelse.(autocorr_p10 .== Autocor_p10_benchmark,1,0)
        @test sum(valid) == length(valid)-124  # 124 bar lead in period
        valid = ifelse.(autocorr_p48 .== Autocor_p48_benchmark,1,0)
        @test sum(valid) == length(valid)-185  # 185 bar lead in period
    end

    @testset "Single Lag Autocorrelation - Equation 8-3" begin
        filename = joinpath(dirname(@__FILE__), "test_8-3_Autocorrelation_Indicator.csv")
        test = CSV.read(filename)
        dat = Float64.(test[:, :x])
        Autocor_p10_benchmark = Float64.(test[:cor_10])
        Autocor_ind_out = AutoCorrelationIndicator(dat)
        autocorr_p10 = Autocor_ind_out[:,10]
        autocorr_p10 = round.(autocorr_p10;digits=4) # round same as tradestation output
        valid = ifelse.(autocorr_p10 .== Autocor_p10_benchmark,1,0)
        @test sum(valid) == length(valid)-124  # 124 bar lead in period
    end

    # to do - fix MaxPwr calculation for Autocorrelation Peeriodogram 8-3

    @testset "Autocorrelation Reversals - Equation 8-4" begin
        filename = joinpath(dirname(@__FILE__), "test_8-4_Autocorrelation_Reversals.csv")
        test = CSV.read(filename)
        dat = Float64.(test[:, :x])
        Autocor_reversals_benchmark = Float64.(test[:autocor_reversal])
        Autocor_reversals_out = AutoCorrelationReversals(dat)
        valid = ifelse.(Autocor_reversals_out .== Autocor_reversals_benchmark,1,0)
        @test sum(valid) == length(valid)-110  # 110 bar lead in period
    end

    # TO DO- fix MaxPwr calculation for DFT Spectral Estimate

    # TO DO -  Complete Comb Filter Spectral Estimate - Equation 10-1

end
