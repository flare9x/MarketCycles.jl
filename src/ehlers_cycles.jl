##############################
# TO DO
# 5-2. Dominant Cycle Measured by Zero Crossings of the Band-Pass Filter - Validated against TS up to DC portion
# 8-3. Autocorrelation Periodogram - Validated against TS up to Normalization
# Outstanding
# 9-1 onwards
##############################

@doc """
Super Smoother - Equation 3-3
`SuperSmoother(x::Array{Float64}; n::Int64=10)::Array{Float64}`
""" ->
function SuperSmoother(x::Array{Float64}; n::Int64=10)::Array{Float64}
a = exp(-1.414*3.14159 / n)
b = 2 * a * cosd(1.414 * 180 / n)
c2 = b
c3 = -a * a
c1 = 1 - c2 - c3
@assert n<size(x,1) && n>0 "Argument n out of bounds."
Super = zeros(x)
 @inbounds for i = 3:length(x)
Super[i] = c1 * (x[i] + x[i-1]) / 2 + c2 * Super[i-1] + c3 * Super[i-2]
end
return Super
end

@doc """
Decycler - Equation 4-1
`Decycler(x::Array{Float64}; n::Int64=60)::Array{Float64}`
""" ->
function Decycler(x::Array{Float64}; n::Int64=60)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
#Highpass filter cyclic components whose periods are shorter than “cutoff” bars
alpha1 = ((cosd(360 / n) + sind(360 / n) - 1)) / (cosd(360 / n))
Decycle = zeros(x)
 @inbounds for i in 2:length(x)
     Decycle[i] = (alpha1 / 2)*(x[i] + x[i-1]) + (1- alpha1)*Decycle[i-1]
end
return Decycle
end

@doc """
Decycle Oscillator - Equation 4-2
`Decycle_OSC(x::Array{Float64}; n1::Int64=30, n2::Int64=60)::Array{Float64}`
""" ->
function Decycle_OSC(x::Array{Float64}; n1::Int64=30, n2::Int64=60)::Array{Float64}
        @assert n2<size(x,1) && n2>0 "Argument n out of bounds."
alpha1 = (cosd(.707*360 / n1) + sind(.707*360 / n1) - 1) / cosd(.707*360 / n1)
alpha2 = (cosd(.707*360 / n2) + sind(.707*360 / n2) - 1) / cosd(.707*360 / n2)
HP1 = zeros(x)
HP2 = zeros(x)
Decycle_OSC = zeros(x)
@inbounds for i in 3:length(x)
HP1[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] + x[i-2]) + 2*(1 - alpha1)*HP1[i-1] - (1 - alpha1)* (1 - alpha1)*HP1[i-2]
HP2[i] = (1 - alpha2 / 2)*(1 - alpha2 / 2)*(x[i] - 2*x[i-1] + x[i-2]) + 2*(1 - alpha2)*HP2[i-1] - (1 - alpha2)*(1 - alpha2)*HP2[i-2]
end
Decycle_OSC .= HP2 .- HP1
return Decycle_OSC
end

@doc """
Band Pass Filter - Equation 5-1
`BandPassFilter(x::Array{Float64}; n::Int64=30, bandwidth::Float64=.3)::Array{Float64}`
""" ->
function BandPassFilter(x::Array{Float64}; n::Int64=30, bandwidth::Float64=.3)::Array{Float64}
            @assert n<size(x,1) && n>0 "Argument n out of bounds."
alpha2 = (cosd(.25*bandwidth*360 / n) + sind(.25*bandwidth*360 / n) - 1) / cosd(.25*bandwidth*360 /n)
beta1 = cosd(360 / n);
gamma1 = 1 / cosd(360*bandwidth / n)
alpha1 = gamma1 - sqrt(gamma1*gamma1 - 1)
HP = zeros(x)
BP = zeros(x)
@inbounds for i in 3:length(x)
HP[i] = (1 + alpha2 / 2)*(x[i] - x[i-1]) + (1- alpha2)*HP[i-1]
BP[i] = .5*(1 - alpha1)*(HP[i] - HP[i-2]) + beta1*(1 + alpha1)*BP[i-1] - alpha1*BP[i-2]
end
# Signal
Signal = zeros(x)
Peak = zeros(x)
@inbounds for i in 2:length(BP)
Peak[i] = .991*Peak[i-1]
if abs(BP[i]) > Peak[i]
Peak[i] = abs(BP[i])
    Signal[i] = BP[i] / Peak[i]
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
BP_Trigger = zeros(x)
i=1
@inbounds for i = 2:length(x)
BP_Trigger[i] = (1 + alpha2 / 2)*(Signal[i] - Signal[i-1]) +(1 -alpha2)*BP_Trigger[i-1]
end
return BP_Trigger
end


@doc """
Hurst Coefficient - Equation 6-1
`HurstCoefficient(x::Array{Float64}; n::Int64=30,)::Array{Float64}`
""" ->
function HurstCoefficient(x::Array{Float64}; n::Int64=30,)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
        @assert iseven(n) "n must be an even number."
half_n = Int64(n/2)
a1 = exp(-1.414*3.14159 / 20)  # smoothes by 20 period same as equation 3-3- may wish to make this a argument in the function?
b1 = 2*a1*cosd(1.414*180 / 20) # smoothes by 20 period same as equation 3-3- may wish to make this a argument in the function?
c2 = b1
c3 = -a1*a1
c1 = 1 - c2 - c3
# Find rolling maximum and minimum
HH = zeros(x)
LL = zeros(x)
N3 = zeros(x)
@inbounds for i = n:size(x,1)
                HH[i] = maximum(x[i-n+1:i])
                LL[i] = minimum(x[i-n+1:i])
                N3[i] = (HH[i] - LL[i]) / n
            end
# Rolling min and max half of n
HH = zeros(x)
LL = zeros(x)
N1 = zeros(x)
@inbounds for i = half_n:size(x,1)
                HH[i] = maximum(x[i-half_n+1:i])
                LL[i] = minimum(x[i-half_n+1:i])
                N1[i] = (HH[i] - LL[i]) / half_n
            end

# Set trailing close half of n
HH = [fill(0,half_n); x[1:length(x)-half_n]]
LL = [fill(0,half_n); x[1:length(x)-half_n]]
HH_out = zeros(x)
LL_out = zeros(x)
N2 = zeros(x)
@inbounds for i = half_n:size(x,1)
    HH_out[i] = maximum(HH[i-half_n+1:i])
    LL_out[i] = minimum(LL[i-half_n+1:i])
    N2[i] = (HH_out[i] - LL_out[i])/(half_n)
end

# Hurst
Dimen = zeros(x)
Hurst = zeros(x)
SmoothHurst = zeros(x)
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
HP LP Roofing Filter - Equation 7-1
`HPLPRoofingFilter(x::Array{Float64}; n::Int64=30,)::Array{Float64}`
""" ->
function HPLPRoofingFilter(x::Array{Float64})::Array{Float64}
        @assert n<size(x,1) && n>0 "Argument n out of bounds."
# Highpass filter cyclic components whose periods are shorter than 48 bars
alpha1 = (cosd(360 / 48) + sind(360 / 48) - 1) / cosd(360 / 48)
HP = zeros(x)
@inbounds for i = 2:size(x,1)
    HP[i] = (1 - alpha1 / 2)*(x[i] - x[i-1]) + (1 - alpha1)*HP[i-1]
end
# Smooth with a Super Smoother Filter from equation 3-3
a1 = exp(-1.414*3.14159 / 10)  # may wish to make this an argument in function
b1 = 2*a1*cosd(1.414*180 / 10) # may wish to make this an argument in function
c2 = b1
c3 = -a1*a1
c1 = 1 - c2 - c3
LP_HP_Filt = zeros(x)
@inbounds for i = 3:size(x,1)
LP_HP_Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*LP_HP_Filt[i-1] + c3*LP_HP_Filt[i-2]
end
return LP_HP_Filt
end

@doc """
Zero Mean Roofing Filter - Equation 7-2
P0 = Lag 0
`ZeroMeanRoofingFilterP0(x::Array{Float64})::Array{Float64}`
""" ->
function ZeroMeanRoofingFilterP0(x::Array{Float64})::Array{Float64}
            @assert n<size(x,1) && n>0 "Argument n out of bounds."
# Highpass filter cyclic components whose periods are shorter than 48 bars
alpha1 = (cosd(360 / 48) + sind(360 / 48) - 1) /cosd(360 / 48)
HP = zeros(x)
@inbounds for i = 2:size(data1_c,1)
HP[i] = (1 - alpha1 / 2)*(x[i] - x[i-1]) +(1 - alpha1)*HP[i-1]
end
#Smooth with a Super Smoother Filter from equation 3-3
a1 = exp(-1.414*3.14159 / 10)
b1 = 2*a1*cosd(1.414*180 / 10)
c2 = b1
c3 = -a1*a1
c1 = 1 - c2 - c3
Zero_Mean_Filt = zeros(x)
Zero_Mean_Filt2 = zeros(x)
@inbounds for i = 3:size(x,1)
Zero_Mean_Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Zero_Mean_Filt[i-1] + c3*Zero_Mean_Filt[i-2]
Zero_Mean_Filt2[i] = (1 - alpha1 / 2)*(Filt[i] - Filt[i-1]) + (1 - alpha1)*Filt2[i-1]
end
return Zero_Mean_Filt
end

@doc """
Zero Mean Roofing Filter - Equation 7-2
P1 = Lag 1
`ZeroMeanRoofingFilterP0(x::Array{Float64})::Array{Float64}`
""" ->
function ZeroMeanRoofingFilterP1(x::Array{Float64})::Array{Float64}
            @assert n<size(x,1) && n>0 "Argument n out of bounds."
# Highpass filter cyclic components whose periods are shorter than 48 bars
alpha1 = (cosd(360 / 48) + sind(360 / 48) - 1) /cosd(360 / 48)
HP = zeros(x)
@inbounds for i = 2:size(data1_c,1)
HP[i] = (1 - alpha1 / 2)*(x[i] - x[i-1]) +(1 - alpha1)*HP[i-1]
end
#Smooth with a Super Smoother Filter from equation 3-3
a1 = exp(-1.414*3.14159 / 10)
b1 = 2*a1*cosd(1.414*180 / 10)
c2 = b1
c3 = -a1*a1
c1 = 1 - c2 - c3
Zero_Mean_Filt = zeros(x)
Zero_Mean_Filt2 = zeros(x)
@inbounds for i = 3:size(x,1)
Zero_Mean_Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Zero_Mean_Filt[i-1] + c3*Zero_Mean_Filt[i-2]
Zero_Mean_Filt2[i] = (1 - alpha1 / 2)*(Zero_Mean_Filt[i] - Zero_Mean_Filt[i-1]) + (1 - alpha1)*Zero_Mean_Filt2[i-1]
end
return Zero_Mean_Filt2
end

@doc """
Roofing Filter As Indicator - Equation 7-3
`RoofingFilterIndicator(x::Array{Float64}; n1::Int64=40,n2::Int64=80)::Array{Float64}`
""" ->
function RoofingFilterIndicator(x::Array{Float64}; n1::Int64=40,n2::Int64=80)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
    # Highpass filter cyclic components whose periods are shorter than 48 bars
LPPeriod = n1
HPPeriod = n2
alpha1 = (cosd(.707*360 / HPPeriod) + sind(.707*360 /HPPeriod) - 1) / cosd(.707*360 / HPPeriod)
HP = zeros(x)
@inbounds for i = 3:length(x)
    HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] + x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
end
#Smooth with a Super Smoother Filter from equation 3-3
a1 = exp(-1.414*3.14159 / LPPeriod)
b1 = 2*a1*cosd(1.414*180 / LPPeriod)
c2 = b1
c3 = -a1*a1
c1 = 1 - c2 - c3
Roof_filt_Indicator = zeros(x)
@inbounds for i = 3:length(x)
Roof_filt_Indicator[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Roof_filt_Indicator[i-1] + c3*Roof_filt_Indicator[i-2]
end
return Roof_filt_Indicator
end

@doc """
Modified Stochastic - Equation 7-4
`ModifiedStochastic(x::Array{Float64}; n::Int64=20)::Array{Float64}`
""" ->
function ModifiedStochastic(x::Array{Float64}; n::Int64=20)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
#Highpass filter cyclic components whose periods are shorter than 48 bars
alpha1 = (cosd(.707*360 / 48) + sind(.707*360 / 48) - 1) /cosd(.707*360 / 48)
HP = zeros(x)
@inbounds for i = 3:size(x,1)
HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1]+ x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 -alpha1)*HP[i-2]
end
# Smooth with a Super Smoother Filter from equation 3-3
a1 = exp(-1.414*3.14159 / 10)
b1 = 2*a1*cosd(1.414*180 / 10)
c2 = b1
c3 = -a1*a1
c1 = 1 - c2 - c3
Filt = zeros(x)
@inbounds for i = 3:size(x,1)
Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
end
# Highest and lowest filt over n width
HighestC = zeros(x)
LowestC = zeros(x)
Stoc = zeros(x)
MyStochastic = zeros(x)
@inbounds for i = n:size(x,1)
    HighestC[i] = maximum(Filt[i-n+1:i])
    LowestC[i] = minimum(Filt[i-n+1:i])
    Stoc[i] = (Filt[i] - LowestC[i]) / (HighestC[i] - LowestC[i])
    MyStochastic[i] = c1*(Stoc[i] + Stoc[i-1]) / 2 + c2*MyStochastic[i-1] + c3*MyStochastic[i-2]
end
return MyStochastic
end

@doc """
Modified RSI - Equation 7-5
`ModifiedRSI(x::Array{Float64}; n::Int64=20)::Array{Float64}`
""" ->
function ModifiedRSI(x::Array{Float64}; n::Int64=10)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
#Highpass filter cyclic components whose periods areshorter than 48 bars
alpha1 = (cosd(.707*360 / 48) + sind(.707*360 / 48) - 1) /cosd(.707*360 / 48)
HP = zeros(x)
@inbounds for i =3:size(x,1)
HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 -alpha1)*HP[i-2]
end
#Smooth with a Super Smoother Filter from equation 3-3
a1 = exp(-1.414*3.14159 / 10)
b1 = 2*a1*cosd(1.414*180 / 10)
c2 = b1
c3 = -a1*a1
c1 = 1 - c2 - c3
Filt = zeros(x)
@inbounds for i = 3:size(x,1)
Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
end
ClosesUp = zeros(x)
ClosesDn = zeros(x)
filtdiff = zeros(x)
posDiff= zeros(x)
negDiff= zeros(x)
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
posSum = zeros(x)
negSum = zeros(x)
denom = zeros(x)
rsi= zeros(x)
@inbounds for i = n:size(x,1)
    posSum[i] = sum(posDiff[i-n+1:i])
    negSum[i] = sum(negDiff[i-n+1:i])
     denom[i] = posSum[i]+negSum[i]
end

# RSI
MyRSI = zeros(x)
@inbounds for i = 3:size(x,1)
if denom != 0 && denom[i-1] != 0
    MyRSI[i] = c1*(posSum[i] /denom[i] + posSum[i-1] / denom[i-1]) / 2 + c2*MyRSI[i-1] +c3*MyRSI[i-2]
end
end
return MyRSI
end

@doc """
Autocorrelation - Equation 8-2
Computes Matrix for all min_lag:max_lag
`AutoCorrelationIndicator(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48)::Array{Float64}`
""" ->
function AutoCorrelationIndicator(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
#Highpass filter cyclic components whose periods areshorter than 48 bars
alpha1 = (cosd(.707*360 / 48) + sind(.707*360 / 48) - 1) / cosd(.707*360 / 48)
HP = zeros(x)
@inbounds for i = 3:size(x,1)
HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
end
#Smooth with a Super Smoother Filter from equation 3-3
a1 = exp(-1.414*3.14159 / 10)
b1 = 2*a1*cosd(1.414*180 / 10)
c2 = b1
c3 = -a1*a1
c1 = 1 - c2 - c3
Filt = zeros(x)
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
    @inbounds for i = 96:size(x,1)
        AutoCorrOut[i,j] = cor(lagged[i-j+1:i], Filt[i-j+1:i])
        #//Scale each correlation to range between 0 and 1
        AutoCorrOut[i,j]= .5*(AutoCorrOut[i,j] + 1)
        end
    end
    return AutoCorrOut
end

@doc """
Single Lag Autocorrelation - Equation 8-3
Extract A Single Lag Autocorrelation To Use as Standalone Indicator
`SingleLagAutoCorrelationIndicator(x::Array{Float64}; lag::Int64=10)::Array{Float64}`
""" ->
function SingleLagAutoCorrelationIndicator(x::Array{Float64}; lag::Int64=10)::Array{Float64}
    @assert n<size(x,1) && n>0 "Argument n out of bounds."
#Highpass filter cyclic components whose periods areshorter than 48 bars
alpha1 = (cosd(.707*360 / 48) + sind(.707*360 / 48) - 1) / cosd(.707*360 / 48)
HP = zeros(x)
@inbounds for i = 3:size(x,1)
HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
end
#Smooth with a Super Smoother Filter from equation 3-3
a1 = exp(-1.414*3.14159 / 10)
b1 = 2*a1*cosd(1.414*180 / 10)
c2 = b1
c3 = -a1*a1
c1 = 1 - c2 - c3
Filt = zeros(x)
@inbounds for i = 3:size(x,1)
Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
end
#Pearson correlation for specified lag
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
Autocorrelation Periodogram- Equation 8-3
Computes Matrix for all min_lag:max_lag
    # //////////// note from normalization down did not reproduce TS results - reuires revisit - Help wanted
`AutoCorrelationPeriodogram(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48)::Array{Float64}`
""" ->
function AutoCorrelationPeriodogram(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48)::Array{Float64}
        @assert max_lag<size(x,1) && max_lag>0 "Argument n out of bounds."
alpha1 = (cosd(.707*360 / 48) + sind(.707*360 / 48) - 1) / cosd(.707*360 / 48)
HP = zeros(x)
@inbounds for i = 3:size(x,1)
HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
end
#Smooth with a Super Smoother Filter from equation 3-3
a1 = exp(-1.414*3.14159 / 10)
b1 = 2*a1*cosd(1.414*180 / 10)
c2 = b1
c3 = -a1*a1
c1 = 1 - c2 - c3
Filt = zeros(x)
@inbounds for i = 3:size(x,1)
Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
end
#//Pearson correlation for each value of lag
#//Initialize correlation sums
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


cosinePart = Array{Float64}(length(x),max_lag)
sinePart = Array{Float64}(length(x),max_lag)
sqSum = Array{Float64}(length(x),max_lag)

# Calcualte sine and cosine part
cosinePart = Array{Float64}(length(x),max_lag)
sinePart = Array{Float64}(length(x),max_lag)
sqSum = Array{Float64}(length(x),max_lag)
@inbounds for j = min_lag:max_lag
    for k = 3:48
    cosinePart[:,j] .= cosinePart[:,j] .+ Avg_Corr_Out[:,k] .* cosd(370 * k / j)
    sinePart[:,j] .= sinePart[:,j] .+ Avg_Corr_Out[:,k] .* sind(370 * k / j)
    sqSum[:,j] .= cosinePart[:,j].^2 .+ sinePart[:,j].^2
end
end

# Iterate over every i in j and smooth R by the .2 and .8 factors
R = Array{Float64}(length(x),max_lag)
@inbounds for j = min_lag:max_lag
    @inbounds for i = 2:size(x,1)
    R[i,j] = (.2*sqSum[i,j]) * (sqSum[i,j]) + (.8 *R[i-1,j])
end
end
#### validated against TS above ^^^^^ ###############
 ## although followed logic for normalization could not reproduce same result  = revisit.

# Find Maximum Power Level for Normalization
# need to validate this and below!
MaxPwr = Array{Float64}(length(x),max_lag)
#MaxPwr = 0
@inbounds for j = min_lag:max_lag
    @inbounds for i = 2:size(x,1)
            MaxPwr[i,j] = .995*MaxPwr[i-1,j]
if R[i,j] > MaxPwr[i,j]
    MaxPwr[i,j]= R[i,j]
        end
    end
end

Pwr = Array{Float64}(length(x),max_lag)
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
for j = 10:48
Spx .= ifelse.(Pwr[:,j] .>= 0.5, Spx .+ j .* Pwr[:,j],Spx)
Sp .= ifelse.(Pwr[:,j] .>= 0.5,Sp .+ Pwr[:,j],Sp)
end

DominantCycle = zeros(x)
for i = 1:size(x,1)
    if Sp[i] != 0
        DominantCycle[i] = Spx[i] / Sp[i]
    end
end
return DominantCycle
end

@doc """
Autocorrelation Reversals - Equation 8-3
# The indicated reversals are very sensitive to the smoothing of the price data.
# LPLength is made available as an indicator input to decrease or increase the number of indicated reversals as desired.
# The AvgLength parameter is also made available as an indicator because this averaging also impacts the number of indicated reversals.
# Care should be taken when increasing the value of this input because the lag of the indicator increases in direct proportion to the increase of the value of the AvgLength.
# Typical delay of the indicator will be about three bars when the AvgLength parameter is set to a value of 3.
`AutoCorrelationReversals(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48, LPLength::Int64=10, HPLength::Int64=48, AvgLength::Int64=3)::Array{Float64}`
""" ->
function AutoCorrelationReversals(x::Array{Float64}; min_lag::Int64=1, max_lag::Int64=48, LPLength::Int64=10, HPLength::Int64=48, AvgLength::Int64=3)::Array{Float64}
        @assert max_lag<size(x,1) && max_lag>0 "Argument n out of bounds."
# Highpass filter cyclic components whose periods are shorter than 48 bars
alpha1 = (cosd(.707*360 / HPLength) + sind(.707*360 / HPLength) - 1) / cosd(.707*360 / HPLength)
HP = zeros(x)
@inbounds for i = 3:size(x,1)
HP[i] = (1 - alpha1 / 2)*(1 - alpha1 / 2)*(x[i] - 2*x[i-1] +x[i-2]) + 2*(1 - alpha1)*HP[i-1] - (1 - alpha1)*(1 - alpha1)*HP[i-2]
end
#Smooth with a Super Smoother Filter from equation 3-3
a1 = exp(-1.414*3.14159 / LPLength)
b1 = 2*a1*cosd(1.414*180 / LPLength)
c2 = b1
c3 = -a1*a1
c1 = 1 - c2 - c3
Filt = zeros(x)
@inbounds for i = 3:size(x,1)
Filt[i] = c1*(HP[i] + HP[i-1]) / 2 + c2*Filt[i-1] + c3*Filt[i-2]
end
#Pearson correlation for each value of lag
lags = min_lag:max_lag
Avg_Corr_Rev_Out = zeros(size(x,1), max_lag)
@inbounds for j = lags
    # Lag series
        lagged = [fill(0,j); Filt[1:length(Filt)-j]]
        # Roll correlation width of lag and lagged version of itself
    @inbounds for i = 96:size(x,1)
        Avg_Corr_Rev_Out[i,j] = cor(lagged[i-AvgLength+1:i], Filt[i-AvgLength+1:i])
        # Scale each correlation to range between 0 and 1
                Avg_Corr_Rev_Out[i,j] = .5*(Avg_Corr_Rev_Out[i,j] + 1)
        end
    end

    # mark all > .5 and <.5 crossings
    SumDeltas = zeros(size(x,1), max_lag)
@inbounds for j = lags
    @inbounds for i = 96:size(x,1)
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

# Help Wanted
##### DC portion needs to be completed correctly ###
## Data validation vs TS passes with exception of DC
# 5-2. Dominant Cycle Measured by Zero Crossings of the Band-Pass Filter
#function DominantCycle(x::Array{Float64}; n::Int64=20, bandwidth::Float64=.7)::Array{Float64}
#alpha2 = (cosd(.25*bandwidth*360 / n) + sind(.25*bandwidth*360 / n) - 1) / cosd(.25*bandwidth*360 /n)
#HP = zeros(x)
#@inbounds for i = 2:length(x)
#HP[i] = (1 + alpha2 / 2)*(x[i] - x[i-1]) + (1- alpha2)*HP[i-1]
#end
#beta1 = cosd(360 / n)
#gamma1 = 1 / cosd(360*bandwidth / n)
#alpha1 = gamma1 - sqrt(gamma1*gamma1 - 1)
#BP = zeros(x)
#@inbounds for i = 3:length(x)
#    if i == 3
#        BP[i-1] = 0
#        BP[i-2] = 0
#    else
#BP[i] = .5*(1 - alpha1)*(HP[i] - HP[i-2]) + beta1*(1 + alpha1)*BP[i-1] -alpha1*BP[i-2]
#end
#end

#peak = 0.0
#peak_save = zero(x)
#DC = zeros(x)
#Real = zeros(x)
#@inbounds for i = 1:length(BP)
#    peak = .991*peak
#    peak_save[i] = .991*peak
#    if abs(BP[i]) > peak
#        peak = abs(BP[i])
#        peak_save[i] = peak
#    elseif peak != 0.0
#        Real[i] = BP[i] / peak
#                peak_save[i] = peak
#    end
#end

#DC = fill(6.0,length(x))
        # mark all 0 crossings
#        co_cu = zeros(BP)
#@inbounds for i = 2:length(BP)
#           if (Real[i] > 0.0) && (Real[i-1] < 0.0) || (Real[i] < 0.0) && (Real[i-1] > 0.0)# cross over
#               co_cu[i] = 1.0
#           else
#               co_cu[i] = 0.0
#           end
#       end

# Count until each cross over or cross under
#function co_cu_count(x::Array{Float64})::Array{Float64}
#    count_save = zero(x)
#    count = 0
#    for i = 2:length(x)
#    if x[i] == 1.0 && x[i-1] == 0.0
#        count = 1.0 + count
#        count_save[i] = count
#    end
#    if x[i] == 0.0 && x[i-1] == 1.0
#        count = 1.0
#        #count_save[i] = count
#    end
#    if x[i] == 0.0 && x[i-1] == 1.0
#        count = 1.0
#    elseif x[i] == 0.0
#        count = count + 1.0
#        #count_save[i] = count
#    end
#end
#    return count_save
#end

#count = co_cu_count(x)

#function locf(x::Array{Float64})
#    dx = zeros(x)
#    for i in 2:length(x)
#        if i == 2
#            # fill index 1
#            dx[i-1] = x[i-1]
#        end
#    if (x[i] == 0.0 && x[i-1] != 0.0)
#        dx[i] = x[i-1]
#    else
#        dx[i] = x[i]
#    end
#    if (x[i] == 0.0 && x[i-1] == 0.0)
#                dx[i] = dx[i-1]
#    end
#end
#    return dx
#end

# run na locf function
#count = locf(count)
#i=67400 - 13
#count = co_cu_count(co_cu)
#DC .= 2 .* count
#DC = fill(6.0,length(x))
#@inbounds for i = 2:length(count)
#    DC[i] = 2*count[i]
#    if (2*count[i]) > (1.25*DC[i-1])
#        DC[i] = 1.25*DC[i-1]
#    elseif (2*count[i]) < (.8*DC[i-1])
#         DC[i] = .8*DC[i-1]
#     end
 #end
#        return DC
#    end
