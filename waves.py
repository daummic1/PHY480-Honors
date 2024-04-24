import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import time

ε_0 = 8.854187817620389e-12
μ_0 = 4*np.pi*1e-7

### untouched, medium base is fine
class Medium:
    """docstring for Medium"""
    def __init__(self, ε_r, μ_r, σ):
        self._ε_r = ε_r
        self._μ_r = μ_r
        self._σ = σ

    def ε_eq(self, wave):
        """
        dielectric constant [F/m]
        """
        return ε_0 * self._ε_r * (1 + self._σ/(1j*wave._ω*ε_0 * self._ε_r)) 

    def μ_eq(self):
        """
        magnetic constant [H/m]
        """
        return μ_0 * self._μ_r

    def ζ_eq(self, wave):
        """
        characteristic impedance [Ω]
        epsilon_eq:   
        """
        return (self.μ_eq()/self.ε_eq(wave))**0.5

    def type(self, wave):
        U = self._σ / (wave._ω * self._ε_r)
        _type = 'Good conductor' if (U >= 1e2)              else \
                'Dielectric'     if (U < 1e2 and U >= 1e-2) else \
                'Insulator'
        return (U, _type)
    
    def __repr__(self):
        return f"<ε_r={self._ε_r}, μ_r={self._μ_r}, σ={self._σ}>"


### convert mediums and related objects to be list based
class Wave:
    """docstring for Wave object"""

    # pending
    def __init__(self, mediums=None, f=1.8e9, A=10):
        self._f = f             # [Hz]
        self._ω = 2*np.pi*f     # [rad/s]
        self._A = A             # [V/m]
        self._mediums = []
        if mediums == None:
            self._mediums.append(Medium(ε_r=1, μ_r=1, σ=0)) # vacuum is the default
            self._mediums.append(Medium(ε_r=1, μ_r=1, σ=0))
        elif type(mediums) == type(Medium):
            self._mediums.append(Medium(ε_r=1, μ_r=1, σ=0))
            self._mediums.append(mediums)
        elif type(mediums) == list:
            if len(mediums) == 0:
                raise RuntimeError("List must contain medium class member(s)")
            elif len(mediums) == 1:
                self._mediums.append(Medium(ε_r=1, μ_r=1, σ=0))
                self._mediums.append(mediums[0])
            else:
                for medium in mediums:
                    self._mediums.append(medium)
        else:
            raise RuntimeError("Mediums must be passed as a single medium class member or list containing medium class member(s)")
        self._num_mediums = len(self._mediums)

    # pending
    def add_mediums(self, mediums):
        if type(mediums) == type(Medium):
            self._mediums.append(mediums)
            self._mediums += 1
        elif type(mediums) == list:
            for i in mediums:
                if type(i) != type(Medium):
                    raise RuntimeError("All members of list must be of the medium class")
                else:
                    self._mediums.append(i)
            self._num_mediums += len(mediums)
        else:
            raise RuntimeError("Mediums argument must be of the medium class or a list of medium class members")

    # pending
    def get_mediums(self):
        return self._mediums

    # pending
    def remove_mediums(self, medium_indices):
        if type(medium_indices) == int:
            del self._mediums[medium_indices]
            self._num_mediums -= 1
        elif type(medium_indices) == list:
            for i in medium_indices:
                if type(i) != int:
                    raise RuntimeError("All elements of list must be an int index")
                else:
                    del self._mediums[i]
            self._num_mediums -= len(medium_indices)
        else:
            raise RuntimeError("Medium indices must be either an int or list of int")
        
    def k(self, medium):
        """
        wavenumber [1/m]
        """
        return self._ω * np.sqrt(medium.μ_eq() * medium.ε_eq(self))


    # def Γ(self, medium1, medium2):
    #     """
    #     reflection coefficient
    #     """
        
    #     return (medium2.ζ_eq(self)-medium1.ζ_eq(self))/(medium2.ζ_eq(self)+medium1.ζ_eq(self))

    # pending
    def Γ(self, mediums):
        """
        reflection coefficient
        """
        Γc = []
        for i in range(len(mediums)-1):
            Γc.append((mediums[i+1].ζ_eq(self)-mediums[i].ζ_eq(self))/(mediums[i+1].ζ_eq(self)+mediums[i].ζ_eq(self)))
        return Γc

    # def τ(self, medium1, medium2):
    #     """
    #     transmission coefficient
    #     """
    #     return 2*medium2.ζ_eq(self)/(medium2.ζ_eq(self)+medium1.ζ_eq(self))

    # pending
    def τ(self, mediums):
        """
        transmission coefficient
        """
        τc = []
        for i in range(len(mediums)-1):
            τc.append(2*mediums[i+1].ζ_eq(self)/(mediums[i+1].ζ_eq(self)+mediums[i].ζ_eq(self)))
        return τc

    # def δ(self, medium2):
    #     """
    #     skin depth
    #     """
    #     # return np.sqrt(2/(self._ω*medium2.μ_eq()*medium2._σ))
        
    #     α = self.k(medium2).imag
    #     return -1/α if α != 0 else np.inf

    # pending
    def δ(self, mediums):
        """
        skin depth
        """
        δc = []
        for i in range(len(mediums)-1):
            α = self.k(mediums[i+1]).imag
            δc.append(-1/α if α != 0 else np.inf)
        return δc

    def v(self, medium):
        """
        signal velocity [m/s]
        """
        return (1/np.sqrt(medium.ε_eq(self) * medium.μ_eq())).real

    def λ(self, medium):
        """
        wavelength [m]
        """
        return self.v(medium)/self._f

    def power_density_inc(self, medium):
        return 0.5 * 1/abs(medium.ζ_eq(self)) * abs(self._A)**2

    # pending
    def power_density_trans(self, mediums):
        pdt = []
        Γc = self.Γ(mediums)
        for i in range(len(mediums)-1):
            pdt.append(self.power_density_inc(mediums[i]) * (1-(abs(Γc[i]))**2))
        return pdt

    def __repr__(self):
        return f"Wave: f = {self._f} Hz; A = {self._A}"

    # pending
    def print_data(self):
        for i in range(len(self._mediums)):
            t = self._mediums[i].type(self)
            print(f"U_{i+1} := σ_{i+1}/(ω*ε_0*ε_r_{i+1}) = {t[0]:.4g}  ==> medium {i+1} is a(n) \033[92m{t[1]}\x1b[0m".format(i+1))

        for i in range(len(self._mediums)):
            μ = self._mediums[i].μ_eq()
            print(f"μ_eq_{i+1} = {μ:.4g}")

        for i in range(len(self._mediums)):
            ε = self._mediums[i].ε_eq(self)
            print(f"ε_eq_{i+1} = {ε:.4g}")

        for i in range(len(self._mediums)):
            ζ = self._mediums[i].ζ_eq(self)
            print(f"ζ_eq_{i+1} = {ζ:.4g}")

        for i in range(len(self._mediums)):
            k = self.k(self._mediums[i])
            print(f"k_{i+1} = {k:.4g}")

        Γc = self.Γ(self._mediums)
        for i in range(len(self._mediums)-1):
            print(f"Γ_e for mediums {i+1} and {i+2} = {Γc[i]:.4g} = {abs(Γc[i]):.4g} ∠ {np.angle(Γc[i]):.4g}")

        τc = self.τ(self._mediums)
        for i in range(len(self._mediums)-1):
            print(f"τ_e for mediums {i+1} and {i+2} = {τc[i]:.4g} = {abs(τc[i]):.4g} ∠ {np.angle(τc[i]):.4g}")

        δc = self.δ(self._mediums)
        for i in range(len(self._mediums)-1):
            print(f"δ for medium {i+2} = {δc[i]:.4g}")

        print(f"S_1 = {self.power_density_inc(self._mediums[0]):.4g}")

        pdt = self.power_density_trans(self._mediums)
        for i in range(len(self._mediums)-1):        
            print(f"S_{i+2} = {pdt[i]:.4g} = {100*pdt[i]/self.power_density_inc(self._mediums[i]):.4g}% S_{i+1}")
        

    # pending
    def show(self, t, E1_i, ylim):
        fig, ax = plt.subplots(figsize=(10,8))
        fig.set_dpi(100)

        λc = []
        for i in range(len(self._mediums)):
            λc.append(self.λ(self._mediums[i]))
        n = len(self._mediums)
        spacing = 6*λc[0]/n
        L = n*spacing
        z_list = []
        
        # for i in range(len(self._mediums)):
        #     z_list.append(np.linspace(i*spacing - L/2,(i+1)*spacing - L/2,300))

        z_list.append(np.linspace(-spacing,0,300))
        z_list.append(np.linspace(0,spacing,300))
        z_list.append(np.linspace(spacing,2*spacing,300))

        k_list = []
        for i in range(len(self._mediums)):
            k_list.append(self.k(self._mediums[i]))

        Γ_list = self.Γ(self._mediums)
        τ_list = self.τ(self._mediums)

        ei = []
        er = []
        et = []
        
        for l in range(len(self._mediums)-1):
            if l == 0:
                ei.append(lambda z, t, k=k_list[l]: (E1_i(k,-z,t)).real)
                er.append(lambda z, t, k=k_list[l]: (Γ_list[l] * E1_i(k,z,t)).real)
                et.append(lambda k,z,t: (τ_list[l] * E1_i(k,-z,t)).real)
            else:
                ei.append(lambda k,z,t: (τ_list[l-1] * E1_i(k,-z,t)).real)
                er.append(lambda k,z,t: (Γ_list[l] * τ_list[l] * E1_i(k,z-z.max(),t)).real)
                et.append(lambda k,z,t: (τ_list[l] * E1_i(k,-z,t)).real)
                # et.append(lambda k,z,t: τ_list[l] * E1_i(k,-z,t).real/(ei[l](k_list[l],z_list[1][-1], 0)*τ_list[l-1]))

        lines1 = []
        lines2 = []
        lines3 = []
        
        plot_points = []
        for i in range(len(self._mediums)):
            plot_points.append(np.linspace(i*spacing - L/2,(i+1)*spacing - L/2,300))

        for i in range(len(self._mediums)-1):
            if i == 0:
                line1, = ax.plot(plot_points[0],ei[0](z_list[0], t[0]), "--", color='tab:blue', label='incident', linewidth=1)
                line2, = ax.plot(plot_points[0],er[0](z_list[0], t[0]), "-.", color='tab:orange', label='reflected', linewidth=1)
                line3, = ax.plot(plot_points[1],et[0](k_list[1],z_list[1], t[0]), "-", color='tab:purple', label='transmitted', linewidth=1.5)
                
            else:
                # line1, = ax.plot(plot_points[i],ei[i](k_list[i], z_list[1], t[0]), "--", color='tab:blue', linewidth=1)
                line2, = ax.plot(plot_points[i],er[i](k_list[i], z_list[i], t[0]), "-.", color='tab:orange', linewidth=1)
                line3, = ax.plot(plot_points[i+1],et[i](k_list[i+1],z_list[i+1], t[0]), "-", color='tab:purple', linewidth=1.5)

            lines1.append(line1)
            lines2.append(line2)
            lines3.append(line3) 

        media = []
        for i in range(len(self._mediums)):
            media.append([self._mediums[i]._ε_r,self._mediums[i]._μ_r,self._mediums[i]._σ])
        
        string = "Traveling wave"
        for i in range(len(media)):
            string = string + "\nmedium {}: $\epsilon_r$=".format(i+1) + str(media[i][0]) + ", $\mu_r$=" + str(media[i][1]) + ", $\sigma$ =" + str(media[i][2])

        for m in range(len(plot_points)-1):
            if m == 0:
                ax.axvline(plot_points[m+1][0],color='k',alpha=.5,linestyle='--',label='boundary') 
            else:
                ax.axvline(plot_points[m+1][0],color='k',alpha=.5,linestyle='--')

        plt.title(string)
        
        plt.xlabel('space (z)')
        plt.ylabel('E [V/m]')
        plt.legend(loc='upper right')
        plt.ylim(ylim)
        plt.xlim([-L/2, L/2])
        plt.grid(True)

        def animate(i):

            for j in range(len(self._mediums)-1):
                line1 = lines1[j]
                line2 = lines2[j]
                line3 = lines3[j]
                if j == 0:
                    line1.set_data(plot_points[0], ei[0](z_list[0], t[i]))
                    line2.set_data(plot_points[0], er[0](z_list[0], t[i]))
                    line3.set_data(plot_points[1], et[0](k_list[1], z_list[1], t[i]))
                else:
                    # line1.set_data(plot_points[j], ei[j](k_list[j], z_list[1], t[i]))
                    line2.set_data(plot_points[j], er[j](k_list[j], z_list[j], t[i]))
                    line3.set_data(plot_points[j+1], et[j](k_list[j+1], z_list[j+1], t[i]))
                lines1[j] = line1
                lines2[j] = line2
                lines3[j] = line3

            return *lines1, *lines2, *lines3
        anim = animation.FuncAnimation(fig, animate, frames=len(t), interval=40, blit=True)
        plt.show()
        return anim

    # here 
    def save(self, t, E1_i, ylim):
        raise Exception("Unimplemented")
        anim = self.show(t, E1_i, ylim)
        ### Save animation
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=dict(artist='fu'), bitrate=1800)
        anim.save('wave{time.time()}.mp4', writer=writer, dpi=200)

# changed
class Sine(Wave):
    """docstring for Sine"""
    def function(self, k, z, t, A=None):
        if type(A) == type(None):
            return self._A * np.exp(1j*k*z) * np.exp(1j*self._ω*t)
        else:
            return A * np.exp(1j*k*z) * np.exp(1j*self._ω*t)

    def show(self, mediums=[Medium(ε_r=1, μ_r=1, σ=0), Medium(ε_r=2, μ_r=1, σ=.81)]):
        super().show(
            t=np.linspace(0, super().λ(mediums[0])/super().v(mediums[0]), 45),
            E1_i=self.function,
            ylim=[-2*self._A, 2*self._A]
        )

# class Gaussian(Wave):
#     """docstring for Gaussian"""
#     def __init__(self, f=1.8e9, A=10, rms=2.20):
#         super().__init__(f, A)
#         self._rms = rms

#     def function(self, k, z, t):
#         return self._A * 1/np.sqrt(2*np.pi*self._rms**2) * np.exp(-((self._ω*t + k.real*z)**2)/(2 * self._rms**2)) * np.exp(-k.imag*z)

#     def show(self, medium1=Medium(ε_r=1, μ_r=1, σ=0), medium2=Medium(ε_r=1.5, μ_r=1, σ=.21)):
#         peak=self._A/(self._rms*(2*np.pi)**0.5)
#         super().show(
#             t=np.linspace(-.8e-9, 1e-9, 160),
#             E1_i=self.function,
#             ylim=[-peak,1.2*peak]
#         )

# class Rect(Wave):
#     """docstring for Rectangle"""
#     def __init__(self, f=1.8e9, A=10, width=6.5):
#         super().__init__(f, A)
#         self._width = width

#     def function(self, k, z, t):
#         return self._A * (np.heaviside(self._ω*t + k.real*z+self._width,1e-6) - np.heaviside(self._ω*t + k.real*z-self._width, 1e-6)) * np.exp(-k.imag*z)

#     def show(self, medium1=Medium(ε_r=1, μ_r=1, σ=0), medium2=Medium(ε_r=2, μ_r=1, σ=.81)):
#         peak=self._A
#         super().show(
#             t=np.linspace(-.8e-9, 3e-9, 160),
#             E1_i=self.function,
#             ylim=[-peak,1.2*peak]
#         )
