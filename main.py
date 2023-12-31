import matplotlib.pyplot as plt
import numpy as np
from models import*
# Mm = 58.1222
# Tc = 425.125
# Pc = 3.796e6
# acf = 0.201

Mm = [40,50]
Tc = [309.56,280]
Pc = [72.3e5,50e5]
acf = [0.162,0.15]

pg = cubic_EOS(Mm,Tc,Pc,acf,number_moles=[0.5,0.5])
print(pg.specie_fugacity_coefficient(1,50,300))
# print(pg.saturated_pressure(200,N_divs=500))

# Mm = 44
# Tc = 309.56
# Pc = 72.3e5
# acf = 0.162

# n2o = cubic_EOS(Mm,Tc,Pc,acf)

# T = np.linspace(200,275,50)
# for t in T:
#     plt.plot(t,pg.saturated_pressure(t),'rx')
#     plt.plot(t,n2o.saturated_pressure(t),'bx')

# plt.show()

