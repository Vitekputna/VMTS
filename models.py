from scipy.optimize import fsolve
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
    
class van_der_waals(properties_model):

    def __init__(self, molar_mass : float, critical_temperature : float, critical_pressure : float) -> None:
        super().__init__(molar_mass)
        self.critical_temperature = critical_temperature
        self.critical_pressure = critical_pressure
        self.constant_a = ((27/64)*(self.R*critical_temperature)**2)/critical_pressure
        self.constant_b = self.R*critical_temperature/critical_pressure/8

    def pressure(self, density: float, temperature: float) -> float:
        return density*(self.R*temperature/(self.Mm*1e-3-density*self.constant_b) - density*self.constant_a/((self.Mm*1e-3)**2))

    def density(self, pressure : float, temperature : float) -> float:
        func = lambda x : self.pressure(x,temperature)-pressure

        roots = []
        x0 = 0

        #vapor root
        root = fsolve(func,x0=x0)
        roots.append(root[0])
        x0=10*roots[0]

        #middle root
        root = fsolve(func,x0=x0)
        x0=20*roots[0]

        #liquid root
        root = fsolve(func,x0=x0)
        roots.append(root[0])

        return roots

class Peng_Robinson(van_der_waals):

    def __init__(self, molar_mass: float, critical_temperature : float, critical_pressure : float, acentric_factor : float) -> None:
        #source https://wiki.whitson.com/eos/eos_models/pr_eos/
        self.R = 8.314 #J/K.mol
        self.gas_constant = self.R/molar_mass*1000
        self.Mm = molar_mass
    
        self.critical_temperature = critical_temperature
        self.critical_pressure = critical_pressure

        C_a = 0.45724
        C_b = 0.07780

        self.constant_a = (C_a*(self.R*critical_temperature)**2)/critical_pressure
        self.constant_b = C_b*self.R*critical_temperature/critical_pressure
        
        self.acentric_factor = acentric_factor  
        if acentric_factor <= 0.49:
            self.m_PR = 0.37464 + 1.54226*acentric_factor -0.26992*acentric_factor**2
        else:
            self.m_PR = 0.3796+1.485*acentric_factor-0.1644*acentric_factor**2+0.01667*acentric_factor**3

    def get_a(self, temperature : float) -> float:

        return self.constant_a*(1+self.m_PR*(1-np.sqrt(temperature/self.critical_temperature)))**2

    def pressure(self, density: float, temperature: float) -> float:
        Vm = self.Mm*1e-3/density
        a = self.get_a(temperature)

        return self.R*temperature/(Vm-self.constant_b)-a/(Vm*(Vm+self.constant_b) + self.constant_b*(Vm - self.constant_b))
    
    def density(self, pressure : float, temperature : float) -> float:
        # func = lambda x : self.pressure(x,temperature)-pressure
        
        pass        
    
    def compressibility_factor(self, density : float, temperature : float) -> float:
        d1 = 1+np.sqrt(2)
        d2 = 1-np.sqrt(2)

        a = self.get_a(temperature)
        b = self.constant_b

        Vm = self.Mm*1e-3/density
        return 1/(1-b/Vm) - (a*b/Vm)/(self.R*temperature*b*(1+d1*b/Vm)*(1+d2*b/Vm))

    def reduced_residual_helmholz(self, density : float, temperature : float) -> float:
        d1 = 1+np.sqrt(2)
        d2 = 1-np.sqrt(2)

        a = self.get_a(temperature)
        b = self.constant_b

        Vm = self.Mm*1e-3/density
        return -np.log(1-b+Vm)-(a/(self.R*temperature*b*(d1-d2)))*np.log((1+d1*b/Vm)/(1+d2*b/Vm))
    
    def fugacity_coefficient(self, density : float, temperature : float) -> float:
        Z = self.compressibility_factor(density,temperature)
        RAr = self.reduced_residual_helmholz(density,temperature)

        return np.exp(RAr + Z -1 -np.log(Z))
        
class cubic_EOS(properties_model):

    @staticmethod
    def get_m(acentric_factor : float) -> float:

        if acentric_factor <= 0.49:
            return 0.37464 + 1.54226*acentric_factor -0.26992*acentric_factor**2
        else:
            return 0.3796+1.485*acentric_factor-0.1644*acentric_factor**2+0.01667*acentric_factor**3

    def init_pure(self,molar_mass,critical_temperature,critical_pressure,acentric_factor):

        self.Mm = molar_mass
        self.Tc = critical_temperature
        self.Pc = critical_pressure
        self.acentric_factor = acentric_factor

        self.a_c = (self.Omega_a*(self.R*critical_temperature)**2)/critical_pressure
        self.b = self.Omega_b*self.R*critical_temperature/critical_pressure
        
        self.m = self.get_m(acentric_factor)

    def init_mixture(self,molar_mass,critical_temperature,critical_pressure,acentric_factor,molar_fractions):
        
        size = len(molar_fractions)

        self.Mm = np.sum(np.multiply(molar_fractions,molar_mass))
        
        self.Tc = critical_temperature
        self.Pc = critical_pressure
        self.acentric_factor = acentric_factor

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

    def get_a_pure(self, temperature : float) -> float:
        return self.a_c*(1+self.m*(1-np.sqrt(temperature/self.Tc)))**2
    
    def get_a_mixture(self, temperature : float) -> float:
        pass

    def get_a(self, temperature : float):
        
        if self.is_mixture:
            return self.get_a_mixture(temperature)
        else:
            return self.get_a_pure(temperature)

    def get_b_pure(self) -> float:
        return self.b
    
    def get_b_mixture(self) -> float:
        pass

    def get_b(self):
        
        if self.is_mixture:
            return self.get_b_mixture()
        else:
            return self.get_b_pure()
        
    def pressure(self, density: float, temperature: float) -> float:
        Vm = self.Mm*1e-3/density
        a = self.get_a(temperature)
        b = self.get_b()

        return self.R*temperature/(Vm-b)-a/((Vm+self.delta_1*b)*(Vm+self.delta_2*b))



        
