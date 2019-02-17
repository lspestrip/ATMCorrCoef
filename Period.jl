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

plot(freqq, pow/findmax(pow)[1], xlims=(0,0.8), yscale=:log10, label="dir = 90 deg; W = 4")
