using Random
using Statistics
using StatsBase
using Plots
using DSP

import Dates
import Base.+
import Base.-

# include("Point_Struct.jl")
# include("fun.jl")
# include("constants.jl")
# include("integrate.jl")
# include("Corr_Coef.jl")


per=Periodograms.periodogram(Coo, fs=1/dt, window=hamming)
pow = power(per)
freqq = freq(per)

plot(freqq .* 60, pow/(0.1*findmax(pow)[1]), xlims=(0,5), ylims=(0,0.02), label="Raster - dir = 0 deg; W = 3 m/s", xlabel="Freq[Hz] / fs")
