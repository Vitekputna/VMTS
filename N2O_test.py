import matplotlib.pyplot as plt
import numpy as np
from models import*

Mm = [44]
Tc = [309.56]
Pc = [72.3e5]
acf = [0.162]

pg = cubic_EOS(Mm,Tc,Pc,acf)
print(pg.saturated_pressure(150,P_step=50e3))

T_vec = np.linspace(150,300)

for T in T_vec:
    plt.plot(T,pg.saturated_pressure(T),'rx')

plt.show()