from waves import Medium, Sine
f_0     = 1.8e9 # [Hz]
E_0     = 10.0  # [V/m]

#medium = [Medium(ε_r=1, μ_r=1, σ=0),Medium(ε_r=2, μ_r=1, σ=.81),Medium(ε_r=2, μ_r=1, σ=.81)]
medium = [Medium(ε_r=1, μ_r=1, σ=0),Medium(ε_r=2, μ_r=1, σ=.81)]
wave = Sine(mediums=medium,f=f_0, A=E_0)

# wave = Gaussian(rms=1.3)
# wave = Rect(width=4)

# wave.add_mediums(medium1=free_space, medium2=medium2)
wave.print_data()
wave.show()