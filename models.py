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
            
class cubic_EOS(properties_model):
# Initialization
    @staticmethod
    def get_m(acentric_factor : float) -> float:

        if acentric_factor <= 0.49:
            return 0.37464 + 1.54226*acentric_factor -0.26992*acentric_factor**2
        else:
            return 0.3796+1.485*acentric_factor-0.1644*acentric_factor**2+0.01667*acentric_factor**3

    def init_mixture(self,molar_mass,critical_temperature,critical_pressure,acentric_factor,number_moles) -> None:
        
        size = len(number_moles)

        self.Mm = np.sum(np.multiply(number_moles,molar_mass))/self.N_moles
        
        self.X = np.multiply(number_moles,1/self.N_moles)
        self.Tc = critical_temperature
        self.Pc = critical_pressure
        self.acentric_factor = acentric_factor
        self.N = len(self.X)

        self.a_c = [(self.Omega_a*(self.R*critical_temperature[i])**2)/critical_pressure[i] for i in range(size)]
        self.b = [self.Omega_b*self.R*critical_temperature[i]/critical_pressure[i] for i in range(size)]
        self.m = [self.get_m(acentric_factor[i]) for i in range(size)]
            
    def __init__(self, molar_mass, critical_temperature, critical_pressure, acentric_factor, number_moles : list = [1], c = 1) -> None:
        self.R = 8.314

        self.Omega_a = 0.45724
        self.Omega_b = 0.07780

        self.delta_1 = 0.5*(c+1+np.sqrt((c+1)**2+4*c))
        self.delta_2 = 0.5*(c+1-np.sqrt((c+1)**2+4*c))

        self.N_moles = sum(number_moles)
        self.N = len(number_moles)

        self.init_mixture(molar_mass,critical_temperature,critical_pressure,acentric_factor,number_moles)

# Pure and mixture a,b coefficients
    def get_a(self, temperature : float) -> float:
        a_pure = [self.a_c[i]*(1+self.m[i]*(1-np.sqrt(temperature/self.Tc[i])))**2 for i in range(self.N)]

        a_mix = 0

        for i in range(self.N):
            for j in range(self.N):
                a_mix += self.X[i]*self.X[j]*np.sqrt(a_pure[i]*a_pure[j])

        return a_mix
    
    def get_b(self) -> float:
        
        b_mix = 0

        for i in range(self.N):
            for j in range(self.N):
                b_mix += self.X[i]*self.X[j]*0.5*(self.b[i]+self.b[j])

        return b_mix

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

# Helmholz function, derivatives and helper functions
    def __get_D(self, temperature : float) -> float:
        return self.get_a(temperature)*self.N_moles**2
    
    def __get_B(self) -> float:
        return self.N_moles*self.get_b()

    def __get_g(self, volume : float) -> float:
        B = self.__get_B()
        return np.log(1-B/volume)

    def __get_f(self, volume : float) -> float:
        B = self.__get_B()
        V = volume
        R = self.R
        d1 = self.delta_1
        d2 = self.delta_2

        return np.log((1+d1*B/V)/(1+d2*B/V))/(R*B*(d1-d2))

    def reduced_residual_helmholz(self, density : float, temperature : float) -> float:
        volume = self.N_moles*self.Mm*1e-3/density
        T = temperature

        g = self.__get_g(volume)
        D = self.__get_D(T)
        f = self.__get_f(volume)

        return -self.N_moles*g - D*f/T

    def __get_Fn(self, volume : float) -> float:
        return -self.__get_g(volume)
    
    def __get_gB(self, volume : float) -> float:
        return -1/(volume-self.__get_B())

    def __get_fV(self, volume : float) -> float:
        V = volume
        d1 = self.delta_1
        d2 = self.delta_2
        B = self.__get_B()

        return -1/(self.R*(V+d1*B)*(V+d2*B))

    def __get_fB(self, volume : float) -> float:
        f = self.__get_f(volume)
        V = volume
        fV = self.__get_fV(volume)
        B = self.__get_B()

        return -(f+V*fV)/B

    def __get_FB(self, temperature : float, volume : float) -> float:
        n = self.N_moles
        gB = self.__get_gB(volume)
        D = self.__get_D(temperature)
        fB = self.__get_fB(volume)

        return -n*gB-D*fB/temperature

    def __get_FD(self, temperature : float, volume : float) -> float:
        f = self.__get_f(volume)

        return -f/temperature

    def __get_specie_D(self, specie_idx : int, temperature : float) -> float:
        a_pure = [self.a_c[i]*(1+self.m[i]*(1-np.sqrt(temperature/self.Tc[i])))**2 for i in range(self.N)]

        Di = 0
        for i in range(self.N):
            Di += 2*self.N_moles*self.X[i]*np.sqrt(a_pure[specie_idx]*a_pure[i])

        return Di
    
    def __get_specie_B(self, specie_idx : int) -> float:
        B = self.__get_B()
        
        sum = 0
        for i in range(self.N):
            sum += 2*self.N_moles*self.X[i]*0.5*(self.b[i]+self.b[specie_idx])

        return(sum-B)/self.N_moles
        
    def helmholz_dFdn(self, specie_idx : int, temperature : float, volume : float) -> float:
        Fn = self.__get_Fn(volume)
        FB = self.__get_FB(temperature,volume)
        Bi = self.__get_specie_B(specie_idx)
        FD = self.__get_FD(temperature,volume)
        Di = self.__get_specie_D(specie_idx,temperature)

        return Fn+FB*Bi+FD*Di
        
    def specie_fugacity_coefficient(self, specie_idx : int, density : float, temperature : float) -> float:
        Z = self.compressibility_factor(density,temperature)
        V = self.N_moles*self.Mm*1e-3/density
        helmholz_dni = self.helmholz_dFdn(specie_idx,temperature,V)

        return np.exp(helmholz_dni - np.log(Z))
        
    def saturated_pressure(self, temperature : float, P_step : float = 1e5) -> float:
        # Find P at which is phi_l = phi_v for given T

        if self.N > 1:
            print("This function is only for single specie saturated pressure")
            return None

        Tc = np.max(self.Tc)

        if temperature - Tc > 0:
            return None
        elif Tc-temperature < 0.1*Tc:
            P_step = 1e3

        P_new = np.max(self.Pc)
        P_old = P_new
        d_phi = 0

        counter = 0
        while counter < 1e6:
            density = self.density(P_new,temperature)

            if len(density) > 1:
                vapor_density = density[0]
                liquid_density = density[1]

                phi_l = self.specie_fugacity_coefficient(0,liquid_density,temperature)
                phi_v = self.specie_fugacity_coefficient(0,vapor_density,temperature)

                d_phi = phi_l-phi_v

                # print(P_new,d_phi)

                if d_phi*last_d_phi < 0:
                    return P_old - last_d_phi*(P_new-P_old)/(d_phi-last_d_phi)
                
            P_old = P_new
            P_new = P_new-P_step

            if P_new < 0:
                P_new = 1e3

            last_d_phi = d_phi
            counter += 1

    def PT_flash(self, pressure : float, temperature : float) -> float:

        def residual(X,args) -> float:
            z = args[0]
            K = args[1]

            res = 0
            for i in range(len(K)):
                res += z[i]*(K[i]-1)/(1-X+X*K[i])
    
            return res

        # Compute vapor and liquid densities
        density = self.density(pressure,temperature)
        if len(density) > 1:
            rho_v = density[0]
            rho_l = density[1]
        else: 
            return None

        # Compute equilibrium factors
        K = []
        for i in range(self.N):
            phi_v = self.specie_fugacity_coefficient(i,rho_v,temperature)
            phi_l = self.specie_fugacity_coefficient(i,rho_l,temperature)
            K.append(phi_l/phi_v)

        # Check validity
        if residual(0,[self.X,K]) > 0 and residual(1,[self.X,K]) > 0: # Superheated vapor
            return 1
        elif residual(0,[self.X,K]) < 0 and residual(1,[self.X,K]) < 0: # Subcooled liquid
            return 0
        else:
            args = [self.X,K]
            return fsolve(residual,0.5,args)[0]

