import matplotlib.pyplot as plt
import numpy as np
from CRSprop import CRSprop
from models import*
Mm = 44
Tc = 309.56
Pc = 72.3e5
acf = 0.162

Mm = [40,30]
Tc = [309.56,350]
Pc = [72.3e5,82e5]
acf = [0.162,0.15]

pg = cubic_EOS(Mm,Tc,Pc,acf,[0.5,0.5])
# print(pg.pressure(200,300))


# idg = ideal_gas(Mm)
# vdW = van_der_waals(Mm,309.56,72.3e5)
# pr = Peng_Robinson(Mm,309.56,72.3e5,0.162)

# print(vdW.density(9e6,300))
# print(pr.density(9e6,300))
# print(pr.compressibility_factor(200,300))
# print(pr.fugacity_coefficient(200,300))

# print(pr.pressure(763,300)-9e6)

# # # props = CRSprop(["N2O"])

# rhos = np.linspace(1,1000,100)
# T = 309.56

# for rho in rhos:
#     plt.plot(rho,idg.pressure(rho,T),'rx')
#     plt.plot(rho,vdW.pressure(rho,T),'k.')
#     plt.plot(rho,pr.pressure(rho,T),'go')

# # print(idg.pressure(rho,T))
# # print(vdW.pressure(rho,T))
# plt.ylim([0,50e6])
# plt.show()

