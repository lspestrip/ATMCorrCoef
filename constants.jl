const λ  = 0.002 #in mm => 150GHz
const θᵦ = (π / 180) * 0.15#   10 primi
const ω₀ = λ / ( π * θᵦ ) # θᵦ is the beam opening angle in radians
const z_atm = 3E4 #depends on the observation site
const Tᵧ = 300 # ground temperature in kelvin
const χ₁_0 = 1 # in m^{-2}
const χ₂_0 = 1 # in m^{-2}
const z_0  = 1E4 # the quote where water vapor vanish
const L₀   = 450 # in m, the turbulance correlation length
const Wind_direction = 1*((π/2))
const Wind_Intensity = 45#45  #m/s



#Fig 4  ss = 0.02 delta 0= 0.6
