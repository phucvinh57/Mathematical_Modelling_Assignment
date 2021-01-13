import math
from math import sqrt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def msg(a):
    n = 0
    for idx in range(0, 2013):
        n += math.pow((float(df.at[idx, a]) - float(de.at[idx, 'CO2Air_Real'])),2)
    return n/2013

data = pd.read_excel("Output_4.xlsx")
df = pd.DataFrame(data)
da = pd.read_excel("data.xlsx")
de = pd.DataFrame(da)
time = df[['time']]
euler_air = df[['CO2Air_E']]
real_co2 = de[['CO2Air_Real']]
rk4_air = df[['CO2Air_RK4']]

print(msg('CO2Air_E'))
print(msg('CO2Air_RK4'))

plt.plot(time, euler_air, label='Euler')
plt.plot(time, real_co2 , label='Thực Tế')
plt.xlabel('Time')
plt.ylabel('CO2Air')
plt.title("Biểu đồ CO2Air đo bằng phương pháp Euler và CO2Air thực tế")
plt.legend(loc='best')
plt.show()



plt.plot(time, rk4_air, 'b' , label='RK4')
plt.plot(time, real_co2, 'r' ,label='Thực Tế')
plt.xlabel('Time')
plt.ylabel('CO2Air')
plt.title("Biểu đồ CO2Air đo bằng phương pháp RK4 và CO2Air thực tế")
plt.legend(loc='best')
plt.show()

euler_top = df[['CO2Top_E']]
plt.plot(time, euler_top , label='Euler')
plt.xlabel('Time')
plt.ylabel('CO2Top')
plt.title("Biểu đồ CO2Top đo bằng phương pháp Euler")
plt.legend(loc='best')
plt.show()

rk4_top = df[['CO2Top_RK4']]
plt.plot(time, rk4_top , label='RK4')
plt.xlabel('Time')
plt.ylabel('CO2Top')
plt.title("Biểu đồ CO2Top đo bằng phương pháp RK4")
plt.legend(loc='best')
plt.show()