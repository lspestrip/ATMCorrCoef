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
include("Corr_Coef.jl")


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

    pun = Point{Float64}(2500., π/4., az_arr[i], undef, undef, undef )()
    append!(x, pun.x/2000)
    append!(y, pun.y/2000)

    append!(x_w, -1-(45/2000)*i*cos((π/2)+(π/7)) )
    append!(y_w, 0+(45/2000)*i*sin((π/2)+(π/7)) )


    p1=plot(x, y, xlims=(-1000/2000,1000/2000), ylims=(-1,1), arrow=arrow(), legend=false, ylabel="Y", xlabel="X")
    p1=plot!(x_c,y_c, xlims=(-1000/2000,1000/2000), ylims=(-1,1))


    for j = -5:5
        for k = 0:8
            p1=plot!(x_w .+ k*(500/2000) , y_w .- j*(500/2000), arrow=arrow(), legend = false)
        end
    end
    append!(time_arr, time_s[i])
    append!(Coo_arr, Coo[i]/findmax(Coo)[1])

    p2 = plot(time_arr, Coo_arr, xlims=(0,35), ylims=(0,1.2), arrow=arrow(),  ylabel="Autocorrelation C_00^0t [normalized unit]", xlabel="time [s]")

    plot(p1,p2, layout =(1,2))

end every 1
