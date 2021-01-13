import math
from math import sqrt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def mse(a):
    n = 0
    for idx in range(0, nRow):
        n += math.pow((float(df.at[idx, a]) - float(de.at[idx, 'CO2Air_real'])),2)
    return n/nRow

data = pd.read_excel("data_out\Output_4.xlsx")
df = pd.DataFrame(data)
da = pd.read_excel("data_in\data_4.xlsx")
de = pd.DataFrame(da)
time = df[['time']]
euler_air = df[['CO2Air_E']]
real_co2 = de[['CO2Air_real']]
rk4_air = df[['CO2Air_RK4']]
nRow = len(df.index)
print(nRow)
print(mse('CO2Air_E'))
print(mse('CO2Air_RK4'))

plt.plot(time, euler_air, label='Euler')
plt.plot(time, real_co2[:nRow] , label='Thực Tế')
plt.xlabel('Time')
plt.ylabel('CO2Air')
plt.title("Biểu đồ CO2Air đo bằng phương pháp Euler và CO2Air thực tế")
plt.legend(loc='best')
plt.show()



plt.plot(time, rk4_air, 'b' , label='RK4')
plt.plot(time, real_co2[:nRow], 'r' ,label='Thực Tế')
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