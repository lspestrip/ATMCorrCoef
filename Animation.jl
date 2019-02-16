using Random
using Statistics
using StatsBase
using Plots
import Dates
import Base.+
import Base.-

# include("Point_Struct.jl")
# include("fun.jl")
# include("constants.jl")
# include("integrate.jl")
# include("Corr_Coef.jl")


x = Array{Float64}(undef, 0)
y = Array{Float64}(undef, 0)

x_w = Array{Float64}(undef, 0)
y_w = Array{Float64}(undef, 0)

time_arr = Array{Float64}(undef, 0)
Coo_arr  = Array{Float64}(undef, 0)


phi_cir = ((0:1000) ./ 1000) .* π
punti_cir = [Point{Float64}(2500, π/4., i, undef, undef, undef)() for i in phi_cir]

x_c = [i.x/2000 for i in punti_cir ]
y_c = [i.y/2000 for i in punti_cir ]


@gif for i = 1:length(az_arr)

    p1 = plot([az_arr[i],az_arr[i]], [0,1], proj=:polar, linewidth = 4, arrow=arrow(), yaxis=nothing, label="Telescope line of sight"  )
    p1 = plot!([ (π/2) - Wind_direction + π, (π/2) + Wind_direction + π ], [1, 0], proj=:polar, linewidth = 4, arrow=arrow(), label="Wind direction" )

    append!(time_arr, time_s[i])
    append!(Coo_arr, Coo[i]/findmax(Coo)[1])

    p2 = plot(time_arr, Coo_arr, xlims=(0,35), ylims=(0,1.2), ylabel="Autocorrelation [normalized unit]", xlabel="time [s]")
    p2 = plot!(time_arr, Coo_arr, seriestype=:scatter)


    plot(p1, p2, layout =grid(1,2, width=[0.1,0.9]), size=(1366,768))

end every 1
