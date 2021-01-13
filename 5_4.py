import math
from math import sqrt
import pandas as pd
from xlsxwriter import Workbook


# f''VentRoof
def cal_ppfVentRoof(Cd, URoof, ARoof, AFlr, g, hVent, TAir, TOut, Cw, vWind):
    TMeanAir = (TAir + TOut) / 2
    part1 = Cd * URoof * ARoof / (2 * AFlr)
    part2 = g * hVent * (TAir - TOut) / 2 / TMeanAir + Cw * pow(vWind, 2)
    return part1 * sqrt(part2)


# formula 1
def cal_MVCanAir(VECCanAir, VPCan, VPAir):
    return VECCanAir * (VPCan - VPAir)


# formula 2
def cal_VECCanAir(pAir, c_p_Air, LAI, delta_H, y, rb, rs):
    # y is gamma
    return 2 * pAir * c_p_Air * LAI / (delta_H * y * (rb + rs))


def cal_rs(VPCan, VPAir, R_Can, CO2Air):
    rs_min = 82
    c_evap1 = 4.3
    c_evap2 = 0.54
    rf_RCan = (R_Can+c_evap1)/(R_Can+c_evap2)
    R_Can_SP = 5
    S_rs = 1/(1+math.exp(-1*(R_Can-R_Can_SP)))
    c_night_evap3 = 1.1*pow(10, -11)
    c_night_evap4 = 5.2*pow(10, -6)
    c_evap3 = c_night_evap3*(1-S_rs)+c_night_evap3*S_rs
    c_evap4 = c_night_evap4*(1-S_rs)+c_night_evap4*S_rs
    n_mg_ppm = 0.554
    rf_CO2Air = 1+c_evap3*math.pow(n_mg_ppm*CO2Air-200, 2)
    rf_VPCan_VPAir = 1+c_evap4*pow(VPCan - VPAir, 2)
    return rs_min*rf_RCan*rf_CO2Air*rf_VPCan_VPAir


# formula 3
def cal_MVPadAir(pAir, fPad, nPad, xPad, xOut):
    return pAir * fPad * (nPad * (xPad - xOut) + xOut)


# use for formula 3 & formula 10
def cal_fPad(UPad, phiPad, AFlr):
    return UPad * phiPad / AFlr


# formula 4
def cal_MVFogAir(UFog, phiFog, AFlr):
    return UFog * phiFog / AFlr


# formula 5
def cal_MVBlowAir(nHeatVap, UBlow, PBlow, AFlr):
    return nHeatVap * UBlow * PBlow / AFlr


# formula 6
def cal_MVAirThScr(HECAirThScr, VPAir, VPThScr):
    if VPAir <= VPThScr:
        return 0
    else:
        return 6.4 * pow(10, -9) * HECAirThScr * (VPAir - VPThScr)


# use for formula 7
def cal_HECAirThScr(UThScr, TAir, TThScr):
    return 1.7 * UThScr * pow(abs(TAir - TThScr), 0.33)


# formula 7
def cal_MVAirTop(MWater, R, fThScr, VPAir, VPTop, TAir, TTop):
    return (MWater / R) * fThScr * (VPAir / TAir - VPTop / TTop)


# use for formula 6
def cal_fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop):
    a = UThScr * KThScr * pow(abs(TAir - TTop), 2/3)
    PMean_Air = (pAir + pTop) / 2
    b = (1 - UThScr) * pow(g * (1 - UThScr) * abs(pAir - pTop) / (2 * PMean_Air), 1/2)
    return a + b


# formula 8
def cal_MVAirOut(MWater, R, fVentSide, fVentForced, VPAir, VPOut, TAir, TOut):
    return (MWater / R) * (fVentSide + fVentForced) * (VPAir / TAir - VPOut / TOut)


def cal_fVentSide(nInsScr, ppfVentSide, fleakage, UThScr, ppfVentRoofSide, nSide, nSide_Thr):  # use for formula 8
    if nSide >= nSide_Thr:
        return nInsScr * ppfVentSide + 0.5 * fleakage
    else:
        return nInsScr * (UThScr * ppfVentSide + (1 - UThScr) * ppfVentRoofSide * nSide) + 0.5 * fleakage


def cal_ppfVentSide(Cd, USide, ASide, vWind, AFlr, Cw):  # use for fVentSide
    return Cd * USide * ASide * vWind * sqrt(Cw) / (2 * AFlr)

# use for fVentSide
def cal_ppfVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind):
    a = Cd / AFlr
    b = pow(URoof * USide * ARoof * ASide, 2) / (pow(URoof * ARoof, 2) + pow(USide * ASide, 2))
    TMean_Air = (TAir + TOut) / 2
    c = 2 * g * hSideRoof * (TAir - TOut) / TMean_Air
    _d = (URoof * ARoof + USide * ASide) / 2
    d = pow(_d, 2) * Cw * pow(vWind, 2)
    return a * sqrt(b * c + d)


def cal_fVentForced(nInsScr, UVentForced, phiVentForced, AFlr):  # use for formula 8
    return nInsScr * UVentForced * phiVentForced / AFlr


def cal_nInsScr(sInsScr):
    return sInsScr * (2 - sInsScr)


def cal_fleakage(cleakage, vWind):
    if vWind < 0.25:
        return 0.25 * cleakage
    else:
        return vWind * cleakage


# formula 9
def cal_MVTopOut(MWater, R, fVentRoof, VPTop, VPOut, TTop, TOut):
    return (MWater / R) * fVentRoof * (VPTop / TTop - VPOut / TOut)


# use for formula 9
def cal_fVentRoof(nInsScr, fleakage, UThScr, ppfVentRoofSide, nRoof, nSide, nRoof_Thr, ppfVentRoof):
    # nRoof_Thr la nguong Stack
    if nRoof >= nRoof_Thr:
        return nInsScr * ppfVentRoof + 0.5 * fleakage
    else:
        return nInsScr * (UThScr * ppfVentRoof + (1 - UThScr) * ppfVentRoofSide * nSide) + 0.5 * fleakage


# formula 10
def cal_MVAirOut_Pad(fPad, MWater, R, VPAir, TAir):
    return fPad * MWater / R * VPAir / TAir


# formula 11
def cal_MVAirMech(HECAirMech, VPAir, VPMech):
    if VPAir <= VPMech:
        return 0
    else:
        return 6.4 * pow(10, -9) * HECAirMech * (VPAir - VPMech)


# use for formula 11
def cal_HECAirMech(UMechCool, COPMechCool, PMechCool, AFlr, TAir, TMechCool, delta_H, VPAir, VPMechCool):
    A1 = UMechCool * COPMechCool * PMechCool / AFlr
    A2 = -(TAir - TMechCool + 6.4 * pow(10, -9) * delta_H * (VPAir - VPMechCool))
    return A1 / A2


# formula 12
def cal_MVTopCov_in(HECTopCov_in, VPTop, VPCov_in):
    if VPTop <= VPCov_in:
        return 0
    else:
        return 6.4 * pow(10, -9) * HECTopCov_in * (VPTop - VPCov_in)


# use for formula 12
def cal_HECTopCov_in(cHECin, TTop, TCov_in, ACov, AFlr):
    return cHECin * pow(TTop - TCov_in, 0.33) * ACov / AFlr


###### Calculate dx ######
def dxVPAir(VPAir, VPTop):
    rs = cal_rs(VPCan, VPAir, R_Can, CO2Air)
    VECCanAir = cal_VECCanAir(pAir, c_p_Air, LAI, delta_H, y, rb, rs)
    MVCanAir = cal_MVCanAir(VECCanAir, VPCan, VPAir)
    MVPadAir = cal_MVPadAir(pAir, fPad, nPad, xPad, xOut)
    MVFogAir = cal_MVFogAir(UFog, phiFog, AFlr)
    MVAirTop = cal_MVAirTop(MWater, R, fThScr, VPAir, VPTop, TAir, TTop)
    MVAirThScr = cal_MVAirThScr(HECAirThScr, VPAir, VPThScr)
    MVBlowAir = cal_MVBlowAir(nHeatVap, UBlow, PBlow, AFlr)
    MVAirOut = cal_MVAirOut(MWater, R, fVentSide, fVentForced, VPAir, VPOut, TAir, TOut)
    MVAirOut_Pad = cal_MVAirOut_Pad(fPad, MWater, R, VPAir, TAir)
    MVAirMech = cal_MVAirMech(HECAirMech, VPAir, VPMech)
    return (MVCanAir+MVPadAir+MVFogAir+MVBlowAir-MVAirThScr-MVAirTop-MVAirOut-MVAirOut_Pad-MVAirMech)/capVPAir


def dxVPTop(VPAir, VPTop):
    ppfVentRoof = cal_ppfVentRoof(Cd, URoof, ARoof, AFlr, g, hVent, TAir, TOut, Cw, vWind)
    fVentRoof = cal_fVentRoof(nInsScr, fleakage, UThScr, ppfVentRoofSide, nRoof, nSide, nRoof_Thr, ppfVentRoof)
    MVTopOut = cal_MVTopOut(MWater, R, fVentRoof, VPTop, VPOut, TTop, TOut)
    MVAirTop = cal_MVAirTop(MWater, R, fThScr, VPAir, VPTop, TAir, TTop)
    MVTopCov_in = cal_HECTopCov_in(cHECin, TTop, TCov_in, ACov, AFlr)
    return (MVAirTop-MVTopCov_in-MVTopOut)/capVPTop


# Explicit Euler
def euler(VPAir_0, VPTop_0, h, time):
    n = int(time/h)
    VPAir = VPAir_0
    VPTop = VPTop_0

    for idx in range(1, n + 1):
        k = h * dxVPAir(VPAir, VPTop)
        t = h * dxVPTop(VPAir, VPTop)

        VPAir += k
        VPTop += t

        if idx % 300 == 0:
            row = idx // 300
            worksheet.write(row+1, 0, idx)
            worksheet.write(row + 1, 1, VPAir)
            worksheet.write(row + 1, 2, VPTop)

    return VPAir, VPTop


# Explicit Runge-Kutta 4th order
def rk4(VPAir_0, VPTop_0, h, time):
    n = int(time/h)
    VPAir = VPAir_0
    VPTop = VPTop_0

    for idx in range(1, n + 1):
        k1 = h * dxVPAir(VPAir, VPTop)
        t1 = h * dxVPTop(VPAir, VPTop)
        k2 = h * dxVPAir(VPAir + 0.5 * k1, VPTop + 0.5 * k1)
        t2 = h * dxVPTop(VPAir + 0.5 * t1, VPTop + 0.5 * t1)
        k3 = h * dxVPAir(VPAir + 0.5 * k2, VPTop + 0.5 * k2)
        t3 = h * dxVPTop(VPAir + 0.5 * t2, VPTop + 0.5 * t2)
        k4 = h * dxVPAir(VPAir + k3, VPTop + k3)
        t4 = h * dxVPTop(VPAir + t3, VPTop + t3)

        VPAir += (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        VPTop += (1.0 / 6.0) * (t1 + 2 * t2 + 2 * t3 + t4)

        if idx % 300 == 0:
            row = idx // 300
            worksheet.write(row + 1, 4, VPAir)
            worksheet.write(row + 1, 5, VPTop)

    return VPAir, VPTop


##########   main   ##########
def main():
    print('Place: ', df.at[0, 'Place'])
    print('VPAir_0: ', VPAir)
    print('VPTop_0: ', VPTop)

    step = float(input('Input step: '))
    time = float(input('Input time: '))

    worksheet.write(0, 0, 'time')
    worksheet.write(0, 1, 'VPAir_euler')
    worksheet.write(0, 2, 'VPTop_euler')
    worksheet.write(0, 4, 'VPAir_rk4')
    worksheet.write(0, 5, 'VPTop_rk4')
    worksheet.write(1, 0, 0)
    worksheet.write(1, 1, VPAir)
    worksheet.write(1, 2, VPTop)
    worksheet.write(1, 4, VPAir)
    worksheet.write(1, 5, VPTop)

    air, top = euler(VPAir, VPTop, step, time)
    print('\nExplicit Euler')
    print('The VPAir: ', round(air, 10))
    print('The VPTop: ', round(top, 10))

    air, top = rk4(VPAir, VPTop, step, time)
    print('\nExplicit Runge-Kutta 4th order')
    print('The VPAir: ', round(air, 10))
    print('The VPTop: ', round(top, 10))


##########   START READING DATA   ##########
i = int(input())

data = pd.read_excel("data_VP.xlsx")
df = pd.DataFrame(data)

rb = 275
CO2Air = float(df.at[i, 'CO2Air'])
pAir = float(df.at[i, 'pAir'])
LAI = float(df.at[i, 'LAI'])
VPCan = float(df.at[i, 'VPCan'])
R_Can = float(df.at[i, 'RCan'])
delta_H = 2450000
y = 65.8
c_p_Air = 1000

UPad = float(df.at[i, 'UPad'])
phiPad = float(df.at[i, 'phiPad'])
AFlr = float(df.at[i, 'AFlr'])
fPad = UPad*phiPad/AFlr
nPad = float(df.at[i, 'nPad'])
xPad = float(df.at[i, 'xPad'])
xOut = float(df.at[i, 'xOut'])

UFog = float(df.at[i, 'UFog'])
phiFog = 1.39

nHeatVap = 4.43*pow(10, -8)
UBlow = float(df.at[i, 'UBlow'])
PBlow = float(df.at[i, 'PBlow'])

UThScr = float(df.at[i, 'UThScr'])
TAir = float(df.at[i, 'TAir'])
TThScr = float(df.at[i, 'TThScr'])
HECAirThScr = 1.7*UThScr*pow(abs(TAir - TThScr), 0.33)
VPThScr = float(df.at[i, 'VPThScr'])

MWater = float(df.at[i, 'MWater'])
R = float(df.at[i, 'R'])
KThScr = float(df.at[i, 'KThScr'])
TOut = float(df.at[i, 'TOut'])
p_Mean_Air = float(df.at[i, 'p_Mean_Air'])
pOut = float(df.at[i, 'pOut'])
g = float(df.at[i, 'g'])
fThScr = UThScr*KThScr*pow(abs(TAir - TOut), 0.66)+(1-UThScr)/p_Mean_Air*pow((0.5*p_Mean_Air*(1-UThScr)*g*abs(pAir-pOut)), 0.5)
TTop = float(df.at[i, 'TTop'])
sInsScr = float(df.at[i, 'sInsScr'])
nInsScr = cal_nInsScr(sInsScr)
Cd = float(df.at[i, 'Cd'])
USide = float(df.at[i, 'USide'])
ASide = float(df.at[i, 'ASide'])
vWind = float(df.at[i, 'vWind'])
Cw = float(df.at[i, 'Cw'])
URoof = float(df.at[i, 'URoof'])
ARoof = float(df.at[i, 'ARoof'])
hSideRoof = float(df.at[i, 'hSideRoof'])
ppfVentSide = cal_ppfVentSide(Cd, USide, ASide, vWind, AFlr, Cw)
cleakage = float(df.at[i, 'cleakage'])
fleakage = cal_fleakage(cleakage, vWind)
ppfVentRoofSide = cal_ppfVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind)
nSide = float(df.at[i, 'nSide'])
nSide_Thr = float(df.at[i, 'nSide_Thr'])
fVentSide = cal_fVentSide(nInsScr, ppfVentSide, fleakage, UThScr, ppfVentRoofSide, nSide, nSide_Thr)

UVentForced = float(df.at[i, 'UVentForced'])
phiVentForced = float(df.at[i, 'phiVentForced'])
fVentForced = cal_fVentForced(nInsScr, UVentForced, phiVentForced, AFlr)
UMechCool = float(df.at[i, 'UMechCool'])
COPMechCool = float(df.at[i, 'COPMechCool'])
PMechCool = float(df.at[i, 'PMechCool'])
TMechCool = float(df.at[i, 'TMechCool'])
VPAir = float(df.at[i, 'VPAir'])
VPTop = float(df.at[i, 'VPTop'])
VPOut = float(df.at[i, 'VPOut'])
VPMechCool = float(df.at[i, 'VPMechCool'])
HECAirMech = cal_HECAirMech(UMechCool, COPMechCool, PMechCool, AFlr, TAir, TMechCool, delta_H, VPAir, VPMechCool)
VPMech = float(df.at[i, 'VPMech'])
capVPAir = float(df.at[i, 'capVPAir'])
cHECin = float(df.at[i, 'cHECin'])
TCov_in = float(df.at[i, 'TCov_in'])
ACov = float(df.at[i, 'ACov'])
nRoof = float(df.at[i, 'nRoof'])
nRoof_Thr = float(df.at[i, 'nRoof_Thr'])
hVent = float(df.at[i, 'hVent'])

capVPTop = float(df.at[i, 'capVPTop'])
##########   END READING DATA   ##########


# Calculate VPAir and VPTop, then write result in file Output_5.xlsx
wb = Workbook('Output_5.xlsx')
worksheet = wb.add_worksheet()
main()
wb.close()
