using Random
using Statistics
using StatsBase
using Plots
import Dates
import Base.+
import Base.-

include("Point_Struct.jl")
include("fun.jl")
include("constants.jl")
include("integrate.jl")

#################### Integration #####################

calls = Int64(1E3)
it    = 5
dimensions = 6.

Δₐ = 0.4
ss = 0.04 #0.02 rad  / sec

fₛ = ss / Δₐ


Φₛ = π/2         # start point
vₛ = 2*π / 60.   # speed rad/sec
dt = 0.3         # Time interval

Start_Time = 0   #parse(Float64, ARGS[1])
Stop_Time  = 35 #parse(Float64, ARGS[2])
file_name  = "ciao.txt" #ARGS[3]

Step_start = Start_Time / dt
Step_stop  = Stop_Time / dt


time_s = Array{Float64, 1}(undef, 0)
Coo  = Array{Float64, 1}(undef, 0)
std_err = Array{Float64,1}(undef, 0)
el_arr = Array{Float64,1}(undef, 0)
az_arr = Array{Float64,1}(undef, 0)

el = π/(4.5)

# Raster scan
# xu = Point{Float64}(4000.0,  π/(4.5) + θᵦ, Φₛ + θᵦ, undef, undef, undef)
# xl = Point{Float64}(0.,      π/(4.5) - θᵦ,  Φₛ - θᵦ, undef, undef, undef )

xu = Point{Float64}(2000.0,  el + θᵦ, Φₛ + θᵦ, undef, undef, undef)
xl = Point{Float64}(00.,     el - θᵦ, Φₛ - θᵦ, undef, undef, undef)

start = Dates.value(Dates.now())

#hit=Array{Float64}(undef,0)


for t = Step_start:Step_stop
    tempo = Float64(t*dt)
    # xu = Point{Float64}(2000.0,  el + θᵦ, Φₛ + θᵦ, undef, undef, undef)
    # xl = Point{Float64}(00.,     el - θᵦ, Φₛ - θᵦ, undef, undef, undef)
    # Rivedere gli estremi di integrazione.
    # xu_p = Point{Float64}(5035.,   2.8902652413026098,  (Φₛ + vₛ*t*dt) + θᵦ, undef, undef, undef)
    # xl_p = Point{Float64}(0.,     0.12566370614359174,  (Φₛ + vₛ*t*dt) - θᵦ , undef, undef ,undef)

    # Raster scan
    # xu_p = Point{Float64}(4000.,   π/(4.5) + θᵦ,  (Φₛ + (Δₐ/2)*sin(2*π*fₛ*t*dt)) + θᵦ, undef, undef, undef)
    # xl_p = Point{Float64}(0.,      π/(4.5) - θᵦ,  (Φₛ + (Δₐ/2)*sin(2*π*fₛ*t*dt)) - θᵦ , undef, undef ,undef)

    xu_p = Point{Float64}(2000.,   el + θᵦ,  (Φₛ + (Δₐ/2)*sin(2*π*fₛ*t*dt)) + θᵦ, undef, undef, undef)
    xl_p = Point{Float64}(0.,      el - θᵦ,  (Φₛ + (Δₐ/2)*sin(2*π*fₛ*t*dt)) - θᵦ, undef, undef ,undef)

    result, chi_square = Integrate(corr, xl, xu, xl_p, xu_p, calls, it, dimensions, tempo )
    error = (chi_square/(it+1))^0.5
    println(t*dt," ", result, " ± ", error)

    append!(time_s, t*dt)
    append!(Coo, result)
    append!(std_err,error)

    if t*dt % 10 == 0
        global el -= 2*θᵦ
    end
    #println(el, " ", Φₛ + (Δₐ/2)*sin(2*π*fₛ*t*dt))
    append!(el_arr, el)
    append!(az_arr, Φₛ + (Δₐ/2)*sin(2*π*fₛ*t*dt))
end

ord = sortperm(time_s)
time_s = time_s[ord]
Coo  = Coo[ord]


plot(time_s, Coo/findmax(Coo)[1], seriestype=:scatter, ylims=(0,1.1), title="Test", xlabel = "Time [seconds]", ylabel="Autocorrelation Level [Normalized Units]", m=(1, :circle, 1), bg=RGB(1,1,1))#, yerr=(std_err ./ findmax(Coo)[1]))
