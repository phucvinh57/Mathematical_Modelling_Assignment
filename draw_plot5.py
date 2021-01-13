import math
from math import sqrt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def msg(a):
    n = 0
    for idx in range(0, 2013):
        n += math.pow((float(df.at[idx, a]) - float(de.at[idx, 'VPAir_Real'])),2)
    return n/2013

data = pd.read_excel("Output_5.xlsx")
df = pd.DataFrame(data)
da = pd.read_excel("data_5.xlsx")
de = pd.DataFrame(da)
time = df[['time']]
euler_vpair = df[['VPAir_euler']]
real_vpair = de[['VPAir_Real']]
rk4_vpair = df[['VPAir_rk4']]

print(msg('VPAir_euler'))
print(msg('VPAir_rk4'))

plt.plot(time, euler_vpair, label='Euler')
plt.plot(time, real_vpair , label='Thực Tế')
plt.xlabel('Time')
plt.ylabel('VPAir')
plt.title("Biểu đồ VPAir đo bằng phương pháp Euler và VPAir thực tế")
plt.legend(loc='best')
plt.show()

plt.plot(time, rk4_vpair, 'b' , label='RK4')
plt.plot(time, real_vpair, 'r' ,label='Thực Tế')
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