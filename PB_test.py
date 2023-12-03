import matplotlib.pyplot as plt
import numpy as np
from models import cubic_EOS

Mm = [44.0956,58.1222]
Tc = [369.89,425.125]
Pc = [4.2512e6,3.796e6]
acf = [0.1521,0.201]

gas = cubic_EOS(Mm,Tc,Pc,acf,[0.75,0.25]) # Winter mix

P = 2.07e6
T = 372
# print(gas.PT_flash(P,T),gas.density(P,T))
# print(gas.density(P,T))

T_vec = np.linspace(320,390,num=500)
P_vec = np.linspace(2e6,4e6)

# for P in P_vec:
#     Vals = []
#     X_vec = []
#     for T in T_vec:
#         Value = gas.PT_flash(P,T)
#         if Value != None:
#             Vals.append(Value)
#             X_vec.append(T)
    
#     plt.plot(X_vec,Vals)

# plt.show()


mix = np.linspace(0,1)
last_beta = 0
for x in mix:
    Vals = []
    X_vec = []
    gas = cubic_EOS(Mm,Tc,Pc,acf,[x,1-x]) # Winter mix

    for T in T_vec:
        Value = gas.PT_flash(P,T)
        if Value != None:

            if last_beta == 0 and Value > 0:

                Vals.append(T)
                X_vec.append(x)
            
            elif last_beta < 1 and Value == 1:

                Vals.append(T)
                X_vec.append(x)

            last_beta = Value
            
    
    # plt.plot(x,)
    plt.plot(X_vec,Vals,'k.')

plt.xlabel("Propane x[-]")
plt.ylabel("T[K]")
plt.show()