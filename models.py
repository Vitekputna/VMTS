from scipy.optimize import brentq
import numpy as np

class properties_model:

    def __init__(self, molar_mass : float) -> None:
        self.R = 8.314 #J/K.mol
        self.gas_constant = self.R/molar_mass*1000
        self.Mm = molar_mass

    def pressure(self,density : float, temperature : float) -> float:
        pass

    def density(self,pressure : float, temperature : float) -> float:
        pass

class ideal_gas(properties_model):

    def __init__(self, molar_mass : float) -> None:
        super().__init__(molar_mass)

    def pressure(self, density : float, temperature: float) -> float:
        return density*self.gas_constant*temperature

    def density(self, pressure : float, temperature : float) -> float:
        return pressure/(self.gas_constant*temperature)
            
class cubic_EOS(properties_model):
# Initialization
    @staticmethod
    def get_m(acentric_factor : float) -> float:

        if acentric_factor <= 0.49:
            return 0.37464 + 1.54226*acentric_factor -0.26992*acentric_factor**2
        else:
            return 0.3796+1.485*acentric_factor-0.1644*acentric_factor**2+0.01667*acentric_factor**3

    def init_pure(self,molar_mass,critical_temperature,critical_pressure,acentric_factor) -> None:

        self.Mm = molar_mass
        self.Tc = critical_temperature
        self.Pc = critical_pressure
        self.acentric_factor = acentric_factor

        self.a_c = (self.Omega_a*(self.R*critical_temperature)**2)/critical_pressure
        self.b = self.Omega_b*self.R*critical_temperature/critical_pressure
        
        self.m = self.get_m(acentric_factor)

    def init_mixture(self,molar_mass,critical_temperature,critical_pressure,acentric_factor,molar_fractions) -> None:
        
        size = len(molar_fractions)

        self.Mm = np.sum(np.multiply(molar_fractions,molar_mass))
        
        self.X = molar_fractions
        self.Tc = critical_temperature
        self.Pc = critical_pressure
        self.acentric_factor = acentric_factor
        self.N = len(self.X)

        self.a_c = [(self.Omega_a*(self.R*critical_temperature[i])**2)/critical_pressure[i] for i in range(size)]
        self.b = [self.Omega_b*self.R*critical_temperature[i]/critical_pressure[i] for i in range(size)]
        self.m = [self.get_m(acentric_factor[i]) for i in range(size)]
            
    def __init__(self, molar_mass, critical_temperature, critical_pressure, acentric_factor, molar_fractions : list = [], c = 1) -> None:
        self.R = 8.314

        self.Omega_a = 0.45724
        self.Omega_b = 0.07780

        self.delta_1 = 0.5*(c+1+np.sqrt((c+1)**2+4*c))
        self.delta_2 = 0.5*(c+1-np.sqrt((c+1)**2+4*c))

        self.is_mixture = False
        
        if len(molar_fractions) > 1:
            self.is_mixture = True
            self.init_mixture(molar_mass,critical_temperature,critical_pressure,acentric_factor,molar_fractions)
        else:
            self.init_pure(molar_mass,critical_temperature,critical_pressure,acentric_factor)

# Pure and mixture a,b coefficients
    def get_a_pure(self, temperature : float) -> float:
        return self.a_c*(1+self.m*(1-np.sqrt(temperature/self.Tc)))**2
    
    def get_a_mixture(self, temperature : float) -> float:
        a_pure = [self.a_c[i]*(1+self.m[i]*(1-np.sqrt(temperature/self.Tc[i])))**2 for i in range(self.N)]

        a_mix = 0

        for i in range(self.N):
            for j in range(self.N):
                a_mix += self.X[i]*self.X[j]*np.sqrt(a_pure[i]*a_pure[j])

        return a_mix

    def get_a(self, temperature : float):
        
        if self.is_mixture:
            return self.get_a_mixture(temperature)
        else:
            return self.get_a_pure(temperature)

    def get_b_pure(self) -> float:
        return self.b
    
    def get_b_mixture(self) -> float:
        
        b_mix = 0

        for i in range(self.N):
            for j in range(self.N):
                b_mix += self.X[i]*self.X[j]*0.5*(self.b[i]+self.b[j])

        return b_mix

    def get_b(self) -> float:
        
        if self.is_mixture:
            return self.get_b_mixture()
        else:
            return self.get_b_pure()

# Basic properties
    def pressure(self, density: float, temperature: float) -> float:
        Vm = self.Mm*1e-3/density
        a = self.get_a(temperature)
        b = self.get_b()

        return self.R*temperature/(Vm-b)-a/((Vm+self.delta_1*b)*(Vm+self.delta_2*b))

    def compressibility_factor(self, density : float, temperature : float) -> float:
        a = self.get_a(temperature)
        b = self.get_b()

        Vm = self.Mm*1e-3/density
        return 1/(1-b/Vm) - (a*b/Vm)/(self.R*temperature*b*(1+self.delta_1*b/Vm)*(1+self.delta_2*b/Vm))
    
    def density(self, pressure: float, temperature: float) -> list:
        
        a = self.get_a(temperature)
        b = self.get_b()
        
        R = self.R
        d1 = self.delta_1
        d2 = self.delta_2
        T = temperature
        P = pressure

        A = P
        B = P*b*(d1+d2-1)-R*T
        C = P*(b**2)*(d1*d2-d1-d2) -R*T*b*(d1+d2) + a
        D = -P*d1*d2*b**3 -R*T*d1*d2*b**2 -a*b

        a = B/A
        b = C/A
        c = D/A

        Q = (a**2 -3*b)/9
        R = (2*a**3-9*a*b+27*c)/54

        if R**2 < Q**3:
            #3 roots
            O = np.arccos(R/np.sqrt(Q**3))

            x1 = -2*np.sqrt(Q)*np.cos(O/3) -a/3
            x2 = -2*np.sqrt(Q)*np.cos((O+2*np.pi)/3) -a/3
            x3 = -2*np.sqrt(Q)*np.cos((O-2*np.pi)/3) -a/3

            res = [self.Mm/x1/1000,self.Mm/x2/1000,self.Mm/x3/1000]

            return [np.min(res),np.max(res)]
        else:
            A = -np.sign(R)*np.power(abs(R)+np.sqrt(R**2-Q**3),1/3)
            
            if A == 0:
                B = 0
            else:
                B = Q/A

            x1 = A+B - a/3
            return [self.Mm/x1/1000]

# Helmholz function and derivatives
    def reduced_residual_helmholz(self, density : float, temperature : float) -> float:
        a = self.get_a(temperature)
        b = self.get_b()

        Vm = self.Mm*1e-3/density
        return -np.log(1-b/Vm)-(a/(self.R*temperature*b*(self.delta_1-self.delta_2)))*np.log((1+self.delta_1*b/Vm)/(1+self.delta_2*b/Vm))
    # Wrong
    def fugacity_coefficient(self, density : float, temperature : float) -> float:
        Z = self.compressibility_factor(density,temperature)
        RAr = self.reduced_residual_helmholz(density,temperature)

        return np.exp(RAr + Z -1 -np.log(Z))
        
    # Wrong
    def saturated_pressure(self, temperature : float, N_divs : int = 100) -> float:
        # Find P at which is phi_l = phi_v for given T

        Tc = np.max(self.Tc)

        if temperature - Tc > 0:
            return None
        elif Tc-temperature < 0.1*Tc:
            N_divs = 50000

        def delta(x,args):
            P = x
            self = args[0]
            T = args[1]
            density = self.density(P,T)
            phi_l = self.fugacity_coefficient(density[1],T)
            phi_v = self.fugacity_coefficient(density[0],T)
            return phi_l-phi_v


        Pc = np.max(self.Pc)
        P = np.linspace(Pc,1,N_divs)

        d_phi = 0

        for i in range(0,len(P)):
            density = self.density(P[i],temperature)

            if len(density) > 1:
                vapor_density = density[0]
                liquid_density = density[1]

                phi_l = self.fugacity_coefficient(liquid_density,temperature)
                phi_v = self.fugacity_coefficient(vapor_density,temperature)

                d_phi = phi_l-phi_v

                # print(P[i],d_phi)

                if d_phi*last_d_phi < 0:
                    args = [self,temperature]
                    return brentq(delta,P[i-1],P[i],args)

                    

            last_d_phi = d_phi

