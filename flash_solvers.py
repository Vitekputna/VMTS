from models import *
from prettytable import PrettyTable

class PT_flash:

    def __init__(self, Composition, molar_mass, critical_temperature, critical_pressure, acentric_factor, names = []) -> None:
        self.Z = Composition
        self.Mm = molar_mass
        self.Tc = critical_temperature
        self.Pc = critical_pressure
        self.omega = acentric_factor
        self.names = names

        self.__check_validity()

        self.N_species = len(self.Z)

    def check_stability(self, pressure : float, temperature : float) -> bool:
        return True
    
    def output_table(self,beta,x,y):
        T = PrettyTable()

        T.field_names = ["Specie","Liquid mole fraction","Vapor mole fraction"]
        T.align["Liquid mole fraction"] = 'l'
        T.align["Vapor mole fraction"] = 'l'

        if len(self.names) == self.N_species:
            T.align["Specie"] = 'l'
            for i in range(self.N_species):
                T.add_row([self.names[i],x[i],y[i]])       
        else:
            for i in range(self.N_species):
                T.add_row([str(i),x[i],y[i]])
        
        return T

    def solve(self, pressure : float, temperature : float, beta_init : float = 0.5, precision = 1e-6, max_iterations = 1e3, log = False):
        
        # Compute initial estimate of equilibrium factors
        K = self.__initial_equilibrium_factors(pressure,temperature)

        if log:
            print("/-------------------------------/")
            print("PT FLASH CALCULATION:")
            print("pressure =\t" + str(pressure))
            print("temperature =\t" + str(temperature))

        beta = beta_init
        eps = 1
        iteration = 0

        x = np.ones(self.N_species)
        y = np.ones(self.N_species)

        liquid = cubic_EOS(self.Mm,self.Tc,self.Pc,self.omega,x)
        vapor  = cubic_EOS(self.Mm,self.Tc,self.Pc,self.omega,y) 

        while eps > precision and iteration < max_iterations:
            print(K)
            # Compute 
            x = self.__liquid_composition(self.Z,K,beta)
            y = self.__vapor_composition(self.Z,K,beta)

            # Test
            liquid.set_composition(x)
            vapor.set_composition(y)

            # Update equilibrium factors
            liquid_density = max(liquid.density(pressure,temperature))
            vapor_density = min(vapor.density(pressure,temperature))

            for j in range(self.N_species):
                phi_l = liquid.specie_fugacity_coefficient(j,liquid_density,temperature)
                phi_v = vapor.specie_fugacity_coefficient(j,vapor_density,temperature)

                K[j] = phi_l/phi_v

            # Compute beta
            args = [self.Z,K]
            beta_old = beta
            beta = fsolve(self.__residual,beta,args,full_output=True)[0][0]

            eps = np.abs(beta_old-beta)
            iteration += 1
        

        check = self.__check_bounds(self.Z,K)
        if not check[0]:
            
            beta = check[1]

            x = self.__liquid_composition(self.Z,K,beta)
            y = self.__vapor_composition(self.Z,K,beta)

            if log:
                print()
                print("Solution out of bounds")
                print("final error =\t" + str(eps))
                print("number of iterations =\t" + str(iteration))
                print()
                print("Vapor fraction = " + str(beta))
                print(self.output_table(beta,x,y))
                print()

            return beta, x, y

        else:
            if log:
                print()
                print("Solution converged")
                print("final error =\t" + str(eps))
                print("number of iterations =\t" + str(iteration))
                print()
                print("Vapor fraction = " + str(beta))
                print(self.output_table(beta,x,y))
                print()

            return beta, x, y
    
    def __check_validity(self) -> bool:
        if len(self.Z) == len(self.Mm) and len(self.Mm) == len(self.Tc) and len(self.Tc) == len(self.Pc) and len(self.omega):
            return True
        else:
            print("Size not valid")
            exit(ValueError)

    def __initial_equilibrium_factors(self, pressure : float, temperature : float) -> np.array:

        P = pressure
        T = temperature

        K = np.zeros(self.N_species)
        for i in range(self.N_species):
            K[i] = np.exp( np.log(self.Pc[i]/P) + 5.373*(1+self.omega[i])*(1-self.Tc[i]/T) )

        return K

    @staticmethod
    def __check_bounds(Z,K) -> tuple[bool, int]:

        g0 = -1
        for i in range(len(Z)):
            g0 += Z[i]*K[i]

        if g0 < 0:
            return False, 0
        
        g1 = 1
        for i in range(len(Z)):
            g1 += -Z[i]/K[i]

        if g1 > 0:
            return False, 1
        else:
            return True, 0.5
        
    @staticmethod
    def __residual(X,args)->float:
        z = args[0]
        K = args[1]

        res = 0
        for i in range(len(K)):
            res += z[i]*(K[i]-1)/(1-X+X*K[i])

        return res
    
    @staticmethod
    def __vapor_composition(Z, K, beta):
        
        N = len(Z)
        X = np.zeros(N)

        for i in range(N):
            X[i] = K[i]*Z[i]/(1-beta+beta*K[i])

        return X

    @staticmethod
    def __liquid_composition(Z, K, beta):
        
        N = len(Z)
        X = np.zeros(N)

        for i in range(N):
            X[i] = Z[i]/(1-beta+beta*K[i])

        return X