import matplotlib.pyplot as plt
import numpy as np
from models import cubic_EOS

Mm = [44.0956,58.1222]
Tc = [369.89,425.125]
Pc = [4.2512e6,3.796e6]
acf = [0.1521,0.201]

gas = cubic_EOS(Mm,Tc,Pc,acf,[0.75,0.25]) # Winter mix

P = 3e6
T = 372
# print(gas.PT_flash(P,T),gas.density(P,T))
# print(gas.density(P,T))

T_vec = np.linspace(200,350,num=500)
P_vec = [1e5,2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,1e6,1.1e6,1.2e6,1.3e6,1.4e6,1.5e6]

for P in P_vec:
    Vals = []
    X_vec = []
    for T in T_vec:
        Value = gas.PT_flash(P,T)
        if Value != None:
            Vals.append(Value)
            X_vec.append(T)
    
    plt.plot(X_vec,Vals)

plt.show()