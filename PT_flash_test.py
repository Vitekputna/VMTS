import matplotlib.pyplot as plt
import numpy as np
from flash_solvers import PT_flash

Mm = [44.0956,58.1222]
Tc = [369.89,425.125]
Pc = [4.2512e6,3.796e6]
acf = [0.1521,0.201]
names = ["Propane","Butane"]

P = 2.07e6
T = 350

Comp = [0.5,0.5]

flash = PT_flash(Comp,Mm,Tc,Pc,acf,names=names)
flash.solve(P,T,log=True)