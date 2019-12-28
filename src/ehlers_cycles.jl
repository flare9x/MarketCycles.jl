##############################
# TO DO
# 5-2. Dominant Cycle Measured by Zero Crossings of the Band-Pass Filter - Validated against TS up to DC portion
# 8-3. Autocorrelation Periodogram - Validated against TS up to Normalization
# Outstanding
# 9-1 onwards
##############################

using Statistics

@doc """
    SuperSmoother(x::Array{Float64}; n::Int64=10)::Array{Float64}

Super Smoother - Equation 3-3
"""
function SuperSmoother(x::Array{Float64}; n::Int64=10)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
    a = exp(-1.414*3.14159 / n)
    b = 2 * a * cosd(1.414 * 180 / n)
    c2 = b
    c3 = -a * a
    c1 = 1 - c2 - c3
    Super = zeros(size(x,1))
    @inbounds for i = 3:length(x)
        Super[i] = c1 * (x[i] + x[i-1]) / 2 + c2 * Super[i-1] + c3 * Super[i-2]
    end
    return Super
end

@doc """
    Decycler(x::Array{Float64}; n::Int64=60)::Array{Float64}

Decycler - Equation 4-1
"""
function Decycler(x::Array{Float64}; n::Int64=60)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
    #Highpass filter cyclic components whose periods are shorter than “cutoff” bars
    alpha1 = ((cosd(360 / n) + sind(360 / n) - 1)) / (cosd(360 / n))
    Decycle = zeros(size(x,1))
    @inbounds for i in 2:length(x)
        Decycle[i] = (alpha1 / 2)*(x[i] + x[i-1]) + (1- alpha1)*Decycle[i-1]
    end
    return Decycle
end

@doc """
    Decycle_OSC(x::Array{Float64}; LPPeriod::Int64=30, HPPeriod::Int64=60)::Array{Float64}

Decycle Oscillator - Equation 4-2
"""
function Decycle_OSC(x::Array{Float64}; LPPeriod::Int64=30, HPPeriod::Int64=60)::Array{Float64}
    @assert HPPeriod<size(x,1) && HPPeriod>0 "Argument HPPeriod out of bounds."
    alpha1 = (cosd(.707*360 / LPPeriod) + sind(.707*360 / LPPeriod) - 1) / cosd(.707*360 / LPPeriod)
    alpha2 = (cosd(.707*360 / HPPeriod) + sind(.707*360 / HPPeriod) - 1) / cosd(.707*360 / HPPeriod)
    HP1 = zeros(size(x,1))
    HP2 = zeros(size(x,1))
    Decycle_OSC = zeros(size(x,1))
    @inbounds for i in 3:length(x)
        HP1[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] + x[i-2]) + 2*(1 - alpha1)*HP1[i-1] - (1 - alpha1)* (1 - alpha1)*HP1[i-2]
        HP2[i] = (1 - alpha2 / 2)*(1 - alpha2 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha2)*HP2[i-1] - (1 - alpha2)*(1 - alpha2)*HP2[i-2]
    end
    Decycle_OSC .= HP2 .- HP1
    return Decycle_OSC
end

@doc """
    BandPassFilter(x::Array{Float64}; n::Int64=30, bandwidth::Float64=.3)::Array{Float64}

Band Pass Filter - Equation 5-1
"""
function BandPassFilter(x::Array{Float64}; n::Int64=30, bandwidth::Float64=.3)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
    alpha2 = (cosd(.25*bandwidth*360 / n) + sind(.25*bandwidth*360 / n) - 1) / cosd(.25*bandwidth*360 /n)
    beta1 = cosd(360 / n);
    gamma1 = 1 / cosd(360*bandwidth / n)
    alpha1 = gamma1 - sqrt(gamma1*gamma1 - 1)
    HP = zeros(size(x,1))
    BP = zeros(size(x,1))
    @inbounds for i in 3:length(x)
    HP[i] = (1 + alpha2 / 2)*(x[i] - x[i-1]) + (1- alpha2)*HP[i-1]
    BP[i] = .5*(1 - alpha1)*(HP[i] - HP[i-2]) + beta1*(1 + alpha1)*BP[i-1] - alpha1*BP[i-2]
    end
    # Signal
    Signal = zeros(size(x,1))
    Peak = zeros(size(x,1))
    @inbounds for i in 2:length(BP)
        Peak[i] = .991*Peak[i-1]
        if abs(BP[i]) > Peak[i]
            Peak[i] = abs(BP[i])
        if Peak[i] != 0
            Signal[i] = BP[i] / Peak[i]
        end
        else
        Signal[i] = BP[i] / Peak[i]
    end
    end

    # Replace Nan to 0
    @inbounds for i in 1:length(Signal)
    if isnan(Signal[i]) == 1
        Signal[i] = 0.0
    else
        Signal[i] = Signal[i]
    end
    end

    # Trigger
    alpha2 = (cosd(1.5*bandwidth*360 / n) + sind(1.5*bandwidth*360 / n) - 1) / cosd(1.5*bandwidth*360 /n)
    BP_Trigger = zeros(size(x,1))
    i=1
    @inbounds for i = 2:length(x)
        BP_Trigger[i] = (1 + alpha2 / 2)*(Signal[i] - Signal[i-1]) +(1 -alpha2)*BP_Trigger[i-1]
    end
    return BP_Trigger
end

#="""
    TO DO - Help Wanted - DC Portion of Calculation

Zero Crossings of the Band-Pass Filter - Equation 5-2
"""
=#

@doc """
    HurstCoefficient(x::Array{Float64}; n::Int64=30, LPPeriod::Int64=20)::Array{Float64}

Hurst Coefficient - Equation 6-1
"""
function HurstCoefficient(x::Array{Float64}; n::Int64=30, LPPeriod::Int64=20)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
    @assert iseven(n) "n must be an even number."
    # Smooth with a Super Smoother Filter from equation 3-3
    half_n = Int64(n/2)
    a1 = exp(-1.414*3.14159 / LPPeriod)
    b1 = 2*a1*cosd(1.414*180 / LPPeriod)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    # Find rolling maximum and minimum
    HH = zeros(size(x,1))
    LL = zeros(size(x,1))
    N3 = zeros(size(x,1))
    @inbounds for i = n:size(x,1)
    HH[i] = maximum(x[i-n+1:i])
    LL[i] = minimum(x[i-n+1:i])
    N3[i] = (HH[i] - LL[i]) / n
    end
    # Rolling min and max half of n
    HH = zeros(size(x,1))
    LL = zeros(size(x,1))
    N1 = zeros(size(x,1))
    @inbounds for i = half_n:size(x,1)
    HH[i] = maximum(x[i-half_n+1:i])
    LL[i] = minimum(x[i-half_n+1:i])
    N1[i] = (HH[i] - LL[i]) / half_n
    end
    # Set trailing close half of n
    HH = [fill(0,half_n); x[1:length(x)-half_n]]
    LL = [fill(0,half_n); x[1:length(x)-half_n]]
    HH_out = zeros(size(x,1))
    LL_out = zeros(size(x,1))
    N2 = zeros(size(x,1))
    @inbounds for i = half_n:size(x,1)
    HH_out[i] = maximum(HH[i-half_n+1:i])
    LL_out[i] = minimum(LL[i-half_n+1:i])
    N2[i] = (HH_out[i] - LL_out[i])/(half_n)
    end
    # Hurst
    Dimen = zeros(size(x,1))
    Hurst = zeros(size(x,1))
    SmoothHurst = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
    if N1[i] > 0 && N2[i] > 0 && N3[i] > 0
    Dimen[i] = .5*((log(N1[i]+ N2[i]) - log(N3[i])) / log(2) + Dimen[i-1])
    Hurst[i] = 2 - Dimen[i]
    SmoothHurst[i] = c1*(Hurst[i] + Hurst[i-1]) / 2 + c2*SmoothHurst[i-1]+ c3*SmoothHurst[i-2];
    end
    end
    return SmoothHurst
end

@doc """
    HPLPRoofingFilter(x::Array{Float64}; HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}

HP LP Roofing Filter - Equation 7-1
"""
function HPLPRoofingFilter(x::Array{Float64}; HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}
    @assert HPPeriod<size(x,1) && HPPeriod>0 "Argument HPPeriod out of bounds."
    # Highpass filter cyclic components whose periods are shorter than 48 bars
    alpha1 = (cosd(360 / HPPeriod) + sind(360 / HPPeriod) - 1) / cosd(360 / HPPeriod)
    HP = zeros(size(x,1))
    @inbounds for i = 2:size(x,1)
    HP[i] = (1 - alpha1 / 2)*(x[i] - x[i-1]) + (1 - alpha1)*HP[i-1]
    end
    # Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / LPPeriod)  # may wish to make this an argument in function
    b1 = 2*a1*cosd(1.414*180 / LPPeriod) # may wish to make this an argument in function
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    LP_HP_Filt = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        LP_HP_Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*LP_HP_Filt[i-1] + c3*LP_HP_Filt[i-2]
    end
    return LP_HP_Filt
end

@doc """
    ZeroMeanRoofingFilterK0(x::Array{Float64}; HPPeriod::Int64=48, Smooth::Int64=10)::Array{Float64}

Zero Mean Roofing Filter - Lag 0 - Equation 7-2
K0 = Lag 0
# Lag 0 Is Most Responsive
# Ehlers describes using Lag 0 and Lag 1 cross overs/unders as a signal trigger for buying / selling
"""
function ZeroMeanRoofingFilterK0(x::Array{Float64}; HPPeriod::Int64=48, Smooth::Int64=10)::Array{Float64}
    @assert HPPeriod<size(x,1) && HPPeriod>0 "Argument HPPeriod out of bounds."
    # Highpass filter cyclic components whose periods are shorter than 48 bars
    alpha1 = (cosd(360 / HPPeriod) + sind(360 / HPPeriod) - 1) /cosd(360 / HPPeriod)
    HP = zeros(size(x,1))
    @inbounds for i = 2:size(x,1)
        HP[i] = (1 - alpha1 / 2)*(x[i] - x[i-1]) +(1 - alpha1)*HP[i-1]
    end
    #Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / Smooth)
    b1 = 2*a1*cosd(1.414*180 / Smooth)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    Zero_Mean_Filt = zeros(size(x,1))
    Zero_Mean_Filt2 = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        Zero_Mean_Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Zero_Mean_Filt[i-1] + c3*Zero_Mean_Filt[i-2]
        Zero_Mean_Filt2[i] = (1 - alpha1 / 2)*(Zero_Mean_Filt[i] - Zero_Mean_Filt[i-1]) + (1 - alpha1)*Zero_Mean_Filt2[i-1]
    end
    return Zero_Mean_Filt
end

@doc """
    ZeroMeanRoofingFilterK1(x::Array{Float64}; HPPeriod::Int64=48, Smooth::Int64=10)::Array{Float64}

Zero Mean Roofing Filter - Lag 1 - Equation 7-2
K1 = Lag 1
"""
function ZeroMeanRoofingFilterK1(x::Array{Float64}; HPPeriod::Int64=48, Smooth::Int64=10)::Array{Float64}
    @assert HPPeriod<size(x,1) && HPPeriod>0 "Argument HPPeriod out of bounds."
    # Highpass filter cyclic components whose periods are shorter than 48 bars
    alpha1 = (cosd(360 / HPPeriod) + sind(360 / HPPeriod) - 1) /cosd(360 / HPPeriod)
    HP = zeros(size(x,1))
    @inbounds for i = 2:size(x,1)
        HP[i] = (1 - alpha1 / 2)*(x[i] - x[i-1]) +(1 - alpha1)*HP[i-1]
    end
    #Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / Smooth)
    b1 = 2*a1*cosd(1.414*180 / Smooth)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    Zero_Mean_Filt = zeros(size(x,1))
    Zero_Mean_Filt2 = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        Zero_Mean_Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Zero_Mean_Filt[i-1] + c3*Zero_Mean_Filt[i-2]
        Zero_Mean_Filt2[i] = (1 - alpha1 / 2)*(Zero_Mean_Filt[i] - Zero_Mean_Filt[i-1]) + (1 - alpha1)*Zero_Mean_Filt2[i-1]
    end
    return Zero_Mean_Filt2
end

@doc """
    RoofingFilterIndicator(x::Array{Float64}; LPPeriod::Int64=40,HPPeriod::Int64=80)::Array{Float64}

Roofing Filter As Indicator - Equation 7-3
"""
function RoofingFilterIndicator(x::Array{Float64}; LPPeriod::Int64=40,HPPeriod::Int64=80)::Array{Float64}
    @assert HPPeriod<size(x,1) && HPPeriod>0 "Argument HPPeriod out of bounds."
    # Highpass filter cyclic components whose periods are shorter than 48 (n) bars
    alpha1 = (cosd(.707*360 / HPPeriod) + sind(.707*360 /HPPeriod) - 1) / cosd(.707*360 / HPPeriod)
    HP = zeros(size(x,1))
    @inbounds for i = 3:length(x)
    HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] + x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
    end
    #Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / LPPeriod)
    b1 = 2*a1*cosd(1.414*180 / LPPeriod)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    Roof_filt_Indicator = zeros(size(x,1))
    @inbounds for i = 3:length(x)
        Roof_filt_Indicator[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Roof_filt_Indicator[i-1] + c3*Roof_filt_Indicator[i-2]
    end
    return Roof_filt_Indicator
end

@doc """
    ModifiedStochastic(x::Array{Float64}; n::Int64=20, HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}

Modified Stochastic - Equation 7-4
"""
function ModifiedStochastic(x::Array{Float64}; n::Int64=20, HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
    #Highpass filter cyclic components whose periods are shorter than 48 bars
    alpha1 = (cosd(.707*360 / HPPeriod) + sind(.707*360 / HPPeriod) - 1) /cosd(.707*360 / HPPeriod)
    HP = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1]+ x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 -alpha1)*HP[i-2]
    end
    # Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / LPPeriod)
    b1 = 2*a1*cosd(1.414*180 / LPPeriod)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    Filt = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
    end
    # Highest and lowest filt over n width
    HighestC = zeros(size(x,1))
    LowestC = zeros(size(x,1))
    Stoc = zeros(size(x,1))
    MyStochastic = zeros(size(x,1))
    @inbounds for i = n:size(x,1)
    HighestC[i] = maximum(Filt[i-n+1:i])
    LowestC[i] = minimum(Filt[i-n+1:i])
    Stoc[i] = (Filt[i] - LowestC[i]) / (HighestC[i] - LowestC[i])
    MyStochastic[i] = c1*(Stoc[i] + Stoc[i-1]) / 2 + c2*MyStochastic[i-1] + c3*MyStochastic[i-2]
    end
    return MyStochastic
end

@doc """
    ModifiedRSI(x::Array{Float64}; n::Int64=10, HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}

Modified RSI - Equation 7-5
"""
function ModifiedRSI(x::Array{Float64}; n::Int64=10, HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
    # Highpass filter cyclic components whose periods areshorter than 48 bars
    alpha1 = (cosd(.707*360 / HPPeriod) + sind(.707*360 / HPPeriod) - 1) /cosd(.707*360 / HPPeriod)
    HP = zeros(size(x,1))
    @inbounds for i =3:size(x,1)
        HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 -alpha1)*HP[i-2]
    end
    # Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / LPPeriod)
    b1 = 2*a1*cosd(1.414*180 / LPPeriod)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    Filt = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
    end
    ClosesUp = zeros(size(x,1))
    ClosesDn = zeros(size(x,1))
    filtdiff = zeros(size(x,1))
    posDiff= zeros(size(x,1))
    negDiff= zeros(size(x,1))
    # pos and neg diffs
    @inbounds for i = 2:size(x,1)
        # difference
        filtdiff[i] = Filt[i] - Filt[i-1]
        if filtdiff[i] > 0
            posDiff[i] = filtdiff[i]
        elseif filtdiff[i] < 0
            negDiff[i] = abs(filtdiff[i])
            end
        end
        # Running Sums of Filt
        posSum = zeros(size(x,1))
        negSum = zeros(size(x,1))
        denom = zeros(size(x,1))
        rsi = zeros(size(x,1))
        @inbounds for i = n:size(x,1)
            posSum[i] = sum(posDiff[i-n+1:i])
            negSum[i] = sum(negDiff[i-n+1:i])
            denom[i] = posSum[i]+negSum[i]
        end
        # RSI
        MyRSI = zeros(size(x,1))
        @inbounds for i = 3:size(x,1)
            if denom != 0 && denom[i-1] != 0
                MyRSI[i] = c1*(posSum[i] /denom[i] + posSum[i-1] / denom[i-1]) / 2 + c2*MyRSI[i-1] +c3*MyRSI[i-2]
            end
        end
        return MyRSI
end

@doc """
    AutoCorrelationIndicator(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48, HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}

Autocorrelation Indicator - Equation 8-2
"""
function AutoCorrelationIndicator(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48, HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}
    @assert max_lag<size(x,1) && max_lag>0 "Argument max_lag out of bounds."
    # Highpass filter cyclic components whose periods areshorter than 48 bars
    alpha1 = (cosd(.707*360 / HPPeriod) + sind(.707*360 / HPPeriod) - 1) / cosd(.707*360 / HPPeriod)
    HP = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
    end
    # Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / LPPeriod)
    b1 = 2*a1*cosd(1.414*180 / LPPeriod)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    Filt = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
    end
    # Pearson correlation for each value of lag
    lags = min_lag:max_lag
    AutoCorrOut = zeros(size(x,1), max_lag)
    @inbounds for j = lags
    # Lag series
        lagged = [fill(0,j); Filt[1:length(Filt)-j]]
    # Roll correlation width of lag and lagged version of itself
    @inbounds for i = max_lag:size(x,1)
        AutoCorrOut[i,j] = cor(lagged[i-j+1:i], Filt[i-j+1:i])
        # Scale each correlation to range between 0 and 1
        AutoCorrOut[i,j]= .5*(AutoCorrOut[i,j] + 1)
        end
    end
    return AutoCorrOut
end

@doc """
    SingleLagAutoCorrelationIndicator(x::Array{Float64}; lag::Int64=10, HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}

Single Lag Autocorrelation Indicator - Equation 8-2
"""
function SingleLagAutoCorrelationIndicator(x::Array{Float64}; lag::Int64=10, HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}
    @assert lag<size(x,1) && lag>0 "Argument n out of bounds."
    # Highpass filter cyclic components whose periods areshorter than 48 bars
    alpha1 = (cosd(.707*360 / HPPeriod) + sind(.707*360 / HPPeriod) - 1) / cosd(.707*360 / HPPeriod)
    HP = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
    end
    # Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / LPPeriod)
    b1 = 2*a1*cosd(1.414*180 / LPPeriod)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    Filt = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
    end
    # Pearson correlation for specified lag
    AutoCorrOut = zeros(size(x,1))
    # Lag series
        lagged = [fill(0,lag); Filt[1:length(Filt)-lag]]
    # Roll correlation width of lag and lagged version of itself
    @inbounds for i = lag:size(x,1)
        AutoCorrOut[i] = cor(lagged[i-lag+1:i], Filt[i-lag+1:i])
        #Scale each correlation to range between 0 and 1
        AutoCorrOut[i]= .5*(AutoCorrOut[i] + 1)
        end
    return AutoCorrOut
end

@doc """
    AutoCorrelationPeriodogram(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48,HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}

Autocorrelation Periodogram- Equation 8-3
"""
function AutoCorrelationPeriodogram(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48,HPPeriod::Int64=48, LPPeriod::Int64=10)::Array{Float64}
        @assert max_lag<size(x,1) && max_lag>0 "Argument max_lag out of bounds."
        alpha1 = (cosd(.707*360 / HPPeriod) + sind(.707*360 / HPPeriod) - 1) / cosd(.707*360 / HPPeriod)
        HP = zeros(size(x,1))
        @inbounds for i = 3:size(x,1)
            HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
        end
        # Smooth with a Super Smoother Filter from equation 3-3
        a1 = exp(-1.414*3.14159 / LPPeriod)
        b1 = 2*a1*cosd(1.414*180 / LPPeriod)
        c2 = b1
        c3 = -a1*a1
        c1 = 1 - c2 - c3
        Filt = zeros(size(x,1))
        @inbounds for i = 3:size(x,1)
            Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
        end
        # Pearson correlation for each value of lag
        # Initialize correlation sums
        lags = min_lag:max_lag
        avglength = 3
        temp = zeros(size(x,1))
        Avg_Corr_Out = zeros(size(x,1), max_lag)
        @inbounds for j = lags
            # Lag series
            lagged = [fill(0,j); Filt[1:length(Filt)-j]]
        # Roll correlation width of lag and lagged version of itself
        @inbounds for i = 96:size(x,1)
            Avg_Corr_Out[i,j] = cor(lagged[i-avglength+1:i], Filt[i-avglength+1:i])
            end
        end
        # Calcualte sine and cosine part
        cosinePart = zeros(size(x,1), max_lag)
        sinePart = zeros(size(x,1), max_lag)
        sqSum = zeros(size(x,1), max_lag)
        @inbounds for j = min_lag:max_lag
            for k = 3:48
                cosinePart[:,j] .= cosinePart[:,j] .+ Avg_Corr_Out[:,k] .* cosd(370 * k / j)
                sinePart[:,j] .= sinePart[:,j] .+ Avg_Corr_Out[:,k] .* sind(370 * k / j)
                sqSum[:,j] .= cosinePart[:,j].^2 .+ sinePart[:,j].^2
            end
        end
        # Iterate over every i in j and smooth R by the .2 and .8 factors
        R = zeros(size(x,1), max_lag)
        @inbounds for j = min_lag:max_lag
            @inbounds for i = 2:size(x,1)
                R[i,j] = (.2*sqSum[i,j]) * (sqSum[i,j]) + (.8 *R[i-1,j])
            end
        end
#### validated against TS above ^^^^^ ###############
## although followed logic for normalization as below could not reproduce same result  = revisit.

    # Find Maximum Power Level for Normalization
    # need to validate this and below!
    MaxPwr = zeros(size(x,1), max_lag)
    #MaxPwr = 0
    @inbounds for j = min_lag:max_lag
        @inbounds for i = 2:size(x,1)
            MaxPwr[i,j] = .995*MaxPwr[i-1,j]
        if R[i,j] > MaxPwr[i,j]
            MaxPwr[i,j]= R[i,j]
            end
        end
    end

    Pwr = zeros(size(x,1), max_lag)
    @inbounds for j = min_lag:max_lag
        @inbounds for i = 1:size(x,1)
            Pwr[i,j] = R[i,j] / MaxPwr[i,j]
        end
    end

    # Replace Nan to 0
    @inbounds for j = 1:max_lag
        @inbounds for i = 1:size(x,1)
    if isnan(Pwr[i,j]) == 1
        Pwr[i,j] = 0.0
    else
        Pwr[i,j] = Pwr[i,j]
    end
    end
    end

    # Compute the dominant cycle using the CG of the spectrum
    Spx = zeros(size(x,1))
    Sp = zeros(size(x,1))
    for j = min_lag:max_lag
        Spx .= ifelse.(Pwr[:,j] .>= 0.5, Spx .+ j .* Pwr[:,j],Spx)
        Sp .= ifelse.(Pwr[:,j] .>= 0.5,Sp .+ Pwr[:,j],Sp)
    end

    DominantCycle = zeros(size(x,1))
    for i = 1:size(x,1)
        if Sp[i] != 0
            DominantCycle[i] = Spx[i] / Sp[i]
        end
    end
    return DominantCycle
end

@doc """
    AutoCorrelationReversals(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48, LPPeriod::Int64=10, HPPeriod::Int64=48, AvgLength::Int64=3)::Array{Float64}

Autocorrelation Reversals - Equation 8-4

The indicated reversals are very sensitive to the smoothing of the price data.
LPLength is made available as an indicator input to decrease or increase the number of indicated reversals as desired.
The AvgLength parameter is also made available as an indicator because this averaging also impacts the number of indicated reversals.
Care should be taken when increasing the value of this input because the lag of the indicator increases in direct proportion to the increase of the value of the AvgLength.
Typical delay of the indicator will be about three bars when the AvgLength parameter is set to a value of 3.
"""
function AutoCorrelationReversals(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48, LPPeriod::Int64=10, HPPeriod::Int64=48, AvgLength::Int64=3)::Array{Float64}
    @assert max_lag<size(x,1) && max_lag>0 "Argument n out of bounds."
    # Highpass filter cyclic components whose periods are shorter than 48 bars
    alpha1 = (cosd(.707*360 / HPPeriod) + sind(.707*360 / HPPeriod) - 1) / cosd(.707*360 / HPPeriod)
    HP = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
    end
    # Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / LPPeriod)
    b1 = 2*a1*cosd(1.414*180 / LPPeriod)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    Filt = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
    end
    # Pearson correlation for each value of lag
    lags = min_lag:max_lag
    Avg_Corr_Rev_Out = zeros(size(x,1), max_lag)
    @inbounds for j = lags
    # Lag series
        lagged = [fill(0,j); Filt[1:length(Filt)-j]]
    # Roll correlation width of lag and lagged version of itself
    @inbounds for i = AvgLength:size(x,1)
        Avg_Corr_Rev_Out[i,j] = cor(lagged[i-AvgLength+1:i], Filt[i-AvgLength+1:i])
        # Scale each correlation to range between 0 and 1
        Avg_Corr_Rev_Out[i,j] = .5*(Avg_Corr_Rev_Out[i,j] + 1)
        end
    end

    # mark all > .5 and <.5 crossings
    SumDeltas = zeros(size(x,1), max_lag)
    @inbounds for j = lags
        @inbounds for i = AvgLength:size(x,1)
            if (Avg_Corr_Rev_Out[i,j] > 0.5) && (Avg_Corr_Rev_Out[i-1,j] < 0.5) || (Avg_Corr_Rev_Out[i,j] < 0.5) && (Avg_Corr_Rev_Out[i-1,j] > 0.5)
                SumDeltas[i,j] = 1.0
            else
                SumDeltas[i,j] = 0.0
                end
            end
        end

        # Sum across the matrix of all correlation 0.5 crossings
        Reversal = zeros(size(x,1))
        test_sum = zeros(size(x,1))
        @inbounds for i = 1:size(x,1)
            test_sum[i] = sum(SumDeltas[i,:])
            if sum(SumDeltas[i,:]) > 24
                Reversal[i] = 1
            else Reversal[i] =  0
            end
        end
        return Reversal
end

@doc """
DFTS(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48, LPLength::Int64=10, HPLength::Int64=48)::Array{Float64}

Discrete Fourier Transform Sprectral Estimate - Equation 9-1
"""
function DFTS(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48, LPLength::Int64=10, HPLength::Int64=48)::Array{Float64}
        @assert HPLength<size(x,1) && HPLength>0 "Argument n out of bounds."
# Highpass filter cyclic components whose periods are shorter than 48 bars
    alpha1 = (cosd(.707*360 / HPLength) + sind(.707*360 / HPLength) - 1) / cosd(.707*360 / HPLength)
    HP = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
    end
    # Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / LPLength)
    b1 = 2*a1*cosd(1.414*180 / LPLength)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    Filt = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
    end

    # Initialize matrix
    CosinePart = zeros(size(x,1), max_lag)
    SinePart = zeros(size(x,1), max_lag)
    Pwr = zeros(size(x,1), max_lag)
    # This is the DFT
    @inbounds for j = min_lag:max_lag
    @inbounds for k = 1:max_lag
        lagged_filt = [fill(0,k); Filt[1:length(Filt)-k]]
        CosinePart[:,j] .= CosinePart[:,j] .+ lagged_filt .* cosd(360 * k / j) / j
        SinePart[:,j] .= SinePart[:,j] .+ lagged_filt .* sind(360 * k / j) / j
        Pwr[:,j] .= CosinePart[:,j] .* CosinePart[:,j] .+ SinePart[:,j] .* SinePart[:,j]
        end
    end
    # Find Maximum Power Level for Normalization
    # Note difers from TS output
    MaxPwr = zeros(size(x,1), max_lag)
    @inbounds for j = min_lag:max_lag
    @inbounds for i = 2:size(x,1)
    MaxPwr[i,j]  = .995*MaxPwr[i-1,j]
        if Pwr[i,j]  > MaxPwr[i,j]
         MaxPwr[i,j] = Pwr[i,j]
        end
    end
    end

#+_+_+_+_+_+_+_+_+_+_+_+_ unable to validate the below against TS
    #Normalize Power Levels and Convert to Decibels
    @inbounds for j = min_lag:max_lag
        @inbounds for i = 1:size(x,1)
            if MaxPwr[i,j] != 0
                Pwr[i,j] = Pwr[i,j] / MaxPwr[i,j]
                end
            end
        end

    # Compute the dominant cycle using the CG of the spectrum
    Spx = zeros(size(x,1))
    Sp = zeros(size(x,1))
    for j = min_lag:max_lag
        Spx .= ifelse.(Pwr[:,j] .>= 0.5, Spx .+ j .* Pwr[:,j],Spx)
        Sp .= ifelse.(Pwr[:,j] .>= 0.5,Sp .+ Pwr[:,j],Sp)
    end

    DominantCycle = zeros(size(x,1))
    for i = 2:size(x,1)
        if Sp[i] != 0
            DominantCycle[i] = Spx[i] / Sp[i]
        else
            DominantCycle[i] = DominantCycle[i-1]  # if its zero carry forwrd previous value
            end
        end
        return DominantCycle
end

#="""
    TO DO

Comb Filter Spectral Estimate - Equation 10-1
This spectral estimate may be used to adust other indicators such as RSI, Stochastic and CCI
For example of how the indicator adjustments are made see Adaptive Indicators in equation 11-*
"""
=#
@doc """
    AdaptiveRSI(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48,LPLength::Int64=10, HPLength::Int64=48, AvgLength::Int64=3)::Array{Float64}

Adaptive RSI - Equation 11-1
Adjust the RSI by a lookback period half the length of the dominant cycle
"""
function AdaptiveRSI(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48,LPLength::Int64=10, HPLength::Int64=48, AvgLength::Int64=3)::Array{Float64}
    @assert max_lag<size(x,1) && max_lag>0 "Argument max_lag is out of bounds."
    #@assert max_lag<size(x,1) && max_lag>0 "Argument n is out of bounds."
    alpha1 = (cosd(.707*360 / HPLength) + sind(.707*360 / HPLength) - 1) / cosd(.707*360 / HPLength)
    HP = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
    end
    # Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / LPLength)
    b1 = 2*a1*cosd(1.414*180 / LPLength)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    Filt = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
    end
    # Pearson correlation for each value of lag
    # Initialize correlation sums
    lags = min_lag:max_lag
    temp = zeros(size(x,1))
    Avg_Corr_Out = zeros(size(x,1), max_lag)
    @inbounds for j = lags
    # Lag series
        lagged = [fill(0,j); Filt[1:length(Filt)-j]]
        # Roll correlation width of lag and lagged version of itself
    @inbounds for i = 96:size(x,1)
        Avg_Corr_Out[i,j] = cor(lagged[i-AvgLength+1:i], Filt[i-AvgLength+1:i])
        end
    end

    # Calcualte sine and cosine part
    cosinePart = zeros(size(x,1), max_lag)
    sinePart = zeros(size(x,1), max_lag)
    sqSum = zeros(size(x,1), max_lag)
    @inbounds for j = min_lag:max_lag
        @inbounds for k = 3:max_lag
            cosinePart[:,j] .= cosinePart[:,j] .+ Avg_Corr_Out[:,k] .* cosd(370 * k / j)
            sinePart[:,j] .= sinePart[:,j] .+ Avg_Corr_Out[:,k] .* sind(370 * k / j)
            sqSum[:,j] .= cosinePart[:,j].^2 .+ sinePart[:,j].^2
        end
    end

    # Iterate over every i in j and smooth R by the .2 and .8 factors
    R = zeros(size(x,1), max_lag)
    @inbounds for j = min_lag:max_lag
        @inbounds for i = 2:size(x,1)
            R[i,j] = (.2*sqSum[i,j]) * (sqSum[i,j]) + (.8 *R[i-1,j])
        end
    end
    #### validated ^^^^^ ###############
    # Find Maximum Power Level for Normalization
    # need to validate this and below!
    MaxPwr = zeros(size(x,1), max_lag)
    #MaxPwr = 0
    @inbounds for j = min_lag:max_lag
        @inbounds for i = 2:size(x,1)
            MaxPwr[i,j] = .995*MaxPwr[i-1,j]
        if R[i,j] > MaxPwr[i,j]
            MaxPwr[i,j]= R[i,j]
            end
        end
    end
    Pwr = zeros(size(x,1), max_lag)
    @inbounds for j = min_lag:max_lag
        @inbounds for i = 1:size(x,1)
            Pwr[i,j] = R[i,j] / MaxPwr[i,j]
        end
    end

    # Replace Nan to 0
    @inbounds for j = 1:max_lag
        @inbounds for i = 1:size(x,1)
    if isnan(Pwr[i,j]) == 1
        Pwr[i,j] = 0.0
    else
        Pwr[i,j] = Pwr[i,j]
    end
    end
    end

    # Compute the dominant cycle using the CG of the spectrum
    Spx = zeros(size(x,1))
    Sp = zeros(size(x,1))
    for j = min_lag:max_lag
        Spx .= ifelse.(Pwr[:,j] .>= 0.5, Spx .+ j .* Pwr[:,j],Spx)
        Sp .= ifelse.(Pwr[:,j] .>= 0.5,Sp .+ Pwr[:,j],Sp)
    end

    dominant_cycle = zeros(size(x,1))
    for i = 1:size(x,1)
        if Sp[i] != 0
            dominant_cycle[i] = Spx[i] / Sp[i]
        end
    if dominant_cycle[i] < 5
        dominant_cycle[i] = 5
    end
    if dominant_cycle[i] > max_lag
        dominant_cycle[i] = max_lag
    end
    end

    # Adaptive RSI starts here, using half the measured dominant
    ClosesUp = zeros(size(x,1))
    ClosesDn = zeros(size(x,1))
    filtdiff = zeros(size(x,1))
    posDiff= zeros(size(x,1))
    negDiff= zeros(size(x,1))
    # pos and neg diffs
    @inbounds for i = 2:size(x,1)
        # difference
        filtdiff[i] = Filt[i] - Filt[i-1]
        if filtdiff[i] > 0
            posDiff[i] = filtdiff[i]
        elseif filtdiff[i] < 0
            negDiff[i] = abs(filtdiff[i])
            end
        end

        # Running Sums of Filt
        posSum = zeros(size(x,1))
        negSum = zeros(size(x,1))
        denom = zeros(size(x,1))
        rsi= zeros(size(x,1))
        # Set width of look back 50% of the dominant cycle
        n = Int64.(round.(dominant_cycle ./ 2;digits = 0))
        for i = 1:size(n,1)
        if isnan(n[i]) == 1
            n[i] = 2.0
        else
            n[i] == n[i]
        end
    end
        @inbounds for i = n[1]:size(x,1)
            k = n[i]
            posSum[i] = sum(posDiff[i-k+1:i])
            negSum[i] = sum(negDiff[i-k+1:i])
            denom[i] = posSum[i]+negSum[i]
        end

    # RSI
    MyRSI = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
    if denom != 0 && denom[i-1] != 0
        MyRSI[i] = c1*(posSum[i] /denom[i] + posSum[i-1] / denom[i-1]) / 2 + c2*MyRSI[i-1] +c3*MyRSI[i-2]
        end
    end
    return MyRSI
end

@doc """
    AdaptiveStochastic(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48,LPLength::Int64=10, HPLength::Int64=48, AvgLength::Int64=3)::Array{Float64}

Adaptive Stochastic - Equation 11-2
Adjust the stochastic lookback period by the same value as the dominant cycle
"""
function AdaptiveStochastic(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48,LPLength::Int64=10, HPLength::Int64=48, AvgLength::Int64=3)::Array{Float64}
    @assert max_lag<size(x,1) && max_lag>0 "Argument n out of bounds."
    alpha1 = (cosd(.707*360 / HPLength) + sind(.707*360 / HPLength) - 1) / cosd(.707*360 / HPLength)
    HP = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
    end
    # Smooth with a Super Smoother Filter from equation 3-3
    a1 = exp(-1.414*3.14159 / LPLength)
    b1 = 2*a1*cosd(1.414*180 / LPLength)
    c2 = b1
    c3 = -a1*a1
    c1 = 1 - c2 - c3
    Filt = zeros(size(x,1))
    @inbounds for i = 3:size(x,1)
        Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
    end
    # Pearson correlation for each value of lag
    # Initialize correlation sums
    lags = min_lag:max_lag
    temp = zeros(size(x,1))
    Avg_Corr_Out = zeros(size(x,1), max_lag)
    @inbounds for j = lags
    # Lag series
        lagged = [fill(0,j); Filt[1:length(Filt)-j]]
        # Roll correlation width of lag and lagged version of itself
    @inbounds for i = 96:size(x,1)
        Avg_Corr_Out[i,j] = cor(lagged[i-AvgLength+1:i], Filt[i-AvgLength+1:i])
        end
    end

    # Calcualte sine and cosine part
    cosinePart = zeros(size(x,1), max_lag)
    sinePart = zeros(size(x,1), max_lag)
    sqSum = zeros(size(x,1), max_lag)
    @inbounds for j = min_lag:max_lag
        @inbounds for k = 3:max_lag
            cosinePart[:,j] .= cosinePart[:,j] .+ Avg_Corr_Out[:,k] .* cosd(370 * k / j)
            sinePart[:,j] .= sinePart[:,j] .+ Avg_Corr_Out[:,k] .* sind(370 * k / j)
            sqSum[:,j] .= cosinePart[:,j].^2 .+ sinePart[:,j].^2
        end
    end

    # Iterate over every i in j and smooth R by the .2 and .8 factors
    R = zeros(size(x,1), max_lag)
    @inbounds for j = min_lag:max_lag
        @inbounds for i = 2:size(x,1)
            R[i,j] = (.2*sqSum[i,j]) * (sqSum[i,j]) + (.8 *R[i-1,j])
        end
    end
    #### validated ^^^^^ ###############
    # Find Maximum Power Level for Normalization
    # need to validate this and below!
    MaxPwr = zeros(size(x,1), max_lag)
    #MaxPwr = 0
    @inbounds for j = min_lag:max_lag
        @inbounds for i = 2:size(x,1)
            MaxPwr[i,j] = .995*MaxPwr[i-1,j]
        if R[i,j] > MaxPwr[i,j]
            MaxPwr[i,j]= R[i,j]
            end
        end
    end
    Pwr = zeros(size(x,1), max_lag)
    @inbounds for j = min_lag:max_lag
        @inbounds for i = 1:size(x,1)
            Pwr[i,j] = R[i,j] / MaxPwr[i,j]
        end
    end

    # Replace Nan to 0
    @inbounds for j = 1:max_lag
        @inbounds for i = 1:size(x,1)
    if isnan(Pwr[i,j]) == 1
        Pwr[i,j] = 0.0
    else
        Pwr[i,j] = Pwr[i,j]
            end
        end
    end

    # Compute the dominant cycle using the CG of the spectrum
    Spx = zeros(size(x,1))
    Sp = zeros(size(x,1))
    for j = min_lag:max_lag
        Spx .= ifelse.(Pwr[:,j] .>= 0.5, Spx .+ j .* Pwr[:,j],Spx)
        Sp .= ifelse.(Pwr[:,j] .>= 0.5,Sp .+ Pwr[:,j],Sp)
    end

    dominant_cycle = zeros(size(x,1))
    for i = 1:size(x,1)
        if Sp[i] != 0
            dominant_cycle[i] = Spx[i] / Sp[i]
        end
    if dominant_cycle[i] < 10
        dominant_cycle[i] = 10
    end
    if dominant_cycle[i] > max_lag
        dominant_cycle[i] = max_lag
    end
    end
    # Stochastic Computation starts here
    # Highest and lowest filt over same period as dominant cycle
    HighestC = zeros(size(x,1))
    LowestC = zeros(size(x,1))
    Stoc = zeros(size(x,1))
    adaptive_stochastic = zeros(size(x,1))
    n = Int64.(round.(dominant_cycle; digits=0))
    @inbounds for i = n[1]:size(x,1)
        k = n[i]
        HighestC[i] = maximum(Filt[i-k+1:i])
        LowestC[i] = minimum(Filt[i-k+1:i])
        Stoc[i] = (Filt[i] - LowestC[i]) / (HighestC[i] - LowestC[i])
        adaptive_stochastic[i] = c1*(Stoc[i] + Stoc[i-1]) / 2 + c2*adaptive_stochastic[i-1] + c3*adaptive_stochastic[i-2]
    end
    return adaptive_stochastic
end

@doc """
    PFish(h::Array{Float64}, l::Array{Float64}; n::Int64=10)::Array{Float64}

- Fisher Transform - Equation 15-2 (Variation of)
- Price Normalization
- n = length of normalization period
"""
function PFish(h::Array{Float64}, l::Array{Float64}; n::Int64=10)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
    # Fisher Transform Calculation
    price_input = zeros(size(h,1))
    price_input .= (h .+ l) ./ 2
    MaxH = zeros(size(price_input,1))
    MinL = zeros(size(price_input,1))
    @inbounds for i =n:size(price_input,1)
    MaxH[i] =  maximum(price_input[i-n+1:i])
    MinL[i] = minimum(price_input[i-n+1:i])
    end
    Value1 = zeros(size(price_input,1))
    Fish_out = zeros(size(price_input,1))
    @inbounds for i = n:size(price_input,1)
        Value1[i] = .33*2*((price_input[i] - MinL[i])/(MaxH[i] - MinL[i]) - .5) + .67*Value1[i-1]
        if Value1[i] > .99
            Value1[i] = .999
        else
            Value1[i] = Value1[i]
        end
        if Value1[i] < -.99
            Value1[i] = -.999
        else
            Value1[i] = Value1[i]
        end
        Fish_out[i] = .5*log((1 + Value1[i])/(1 - Value1[i])) + .5*Fish_out[i-1]
    end
        return Fish_out
end
