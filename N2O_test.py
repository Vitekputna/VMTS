import matplotlib.pyplot as plt
import numpy as np
from models import*

Mm = [44]
Tc = [309.56]
Pc = [72.3e5]
acf = [0.162]

pg = cubic_EOS(Mm,Tc,Pc,acf)
print(pg.pressure(1039.26,250))
print(pg.density(1622715,250))
print(pg.saturated_pressure(250))