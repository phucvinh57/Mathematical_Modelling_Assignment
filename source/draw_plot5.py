import math
from math import sqrt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def mse(a):
    n = 0
    for idx in range(0, nRow):
        n += math.pow((float(df.at[idx, a]) - float(de.at[idx, 'VPAir_Real'])),2)
    return n/nRow

data = pd.read_excel("../data_out/Output_5.xlsx")
df = pd.DataFrame(data)
da = pd.read_excel("../data_in/data_5.xlsx")
de = pd.DataFrame(da)
time = df[['time']]
euler_vpair = df[['VPAir_euler']]
real_vpair = de[['VPAir_Real']]
rk4_vpair = df[['VPAir_rk4']]
nRow = len(df.index)
print(mse('VPAir_euler'))
print(mse('VPAir_rk4'))

plt.plot(time, euler_vpair, label='Euler')
plt.plot(time, real_vpair[:nRow] , label='Thực Tế')
plt.xlabel('Time')
plt.ylabel('VPAir')
plt.title("Biểu đồ VPAir đo bằng phương pháp Euler và VPAir thực tế")
plt.legend(loc='best')
plt.show()

plt.plot(time, rk4_vpair, 'b' , label='RK4')
plt.plot(time, real_vpair[:nRow], 'r' ,label='Thực Tế')
plt.xlabel('Time')
plt.ylabel('VPAir')
plt.title("Biểu đồ VPAir đo bằng phương pháp RK4 và VPAir thực tế")
plt.legend(loc='best')
plt.show()

euler_top = df[['VPTop_euler']]
plt.plot(time, euler_top , label='Euler')
plt.xlabel('Time')
plt.ylabel('VPTop')
plt.title("Biểu đồ VPTop đo bằng phương pháp Euler")
plt.legend(loc='best')
plt.show()

rk4_top = df[['VPTop_rk4']]
plt.plot(time, rk4_top , label='RK4')
plt.xlabel('Time')
plt.ylabel('VPTop')
plt.title("Biểu đồ VPTop đo bằng phương pháp RK4")
plt.legend(loc='best')
plt.show()