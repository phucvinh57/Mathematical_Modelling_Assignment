import math
from math import sqrt
import pandas as pd
from xlsxwriter import Workbook

# formula 3
def cal_MCBlowAir(nHeatCO2, UBlow, PBlow, AFlr):
    return nHeatCO2 * UBlow * PBlow / AFlr


# formula 4
def cal_MCExtAir(UExtCO2, phiExtCO2, AFlr):
    return UExtCO2 * phiExtCO2 / AFlr


# formula 5
def cal_MCPadAir_1(fPad, CO2Out, CO2Air):
    return fPad * (CO2Out - CO2Air)


def cal_MCPadAir_2(UPad, phiPad, AFlr, CO2Out, CO2Air):
    fPad = UPad * phiPad / AFlr
    return fPad * (CO2Out - CO2Air)


# formula 6
def cal_MCAirTop(fThScr, CO2Air, CO2Top):
    return fThScr * (CO2Air - CO2Top)


# formula 7
def cal_fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop):
    a = UThScr * KThScr * pow(abs(TAir - TTop), 2 / 3)
    PMean_Air = (pAir + pTop) / 2
    b = (1 - UThScr) * pow(g * (1 - UThScr) * abs(pAir - pTop) / (2 * PMean_Air), 1 / 2)
    return a + b


# formula 9
def cal_MCAirOut(fVentSide, fVentForce, CO2Air, CO2Out):
    return (fVentSide + fVentForce) * (CO2Air - CO2Out)


# formula 10
def cal_ppfVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind):
    a = Cd / AFlr
    b = pow(URoof * USide * ARoof * ASide, 2) / (pow(URoof * ARoof, 2) + pow(USide * ASide, 2))
    TMean_Air = (TAir + TOut) / 2
    c = 2 * g * hSideRoof * (TAir - TOut) / TMean_Air
    _d = (URoof * ARoof + USide * ASide) / 2
    d = pow(_d, 2) * Cw * pow(vWind, 2)
    return a * sqrt(b * c + d)


# formula 11
def cal_nInsScr(sInsScr):
    return sInsScr * (2 - sInsScr)


# formula 12
def cal_fleakage(cleakage, vWind):
    if vWind < 0.25:
        return 0.25 * cleakage
    else:
        return vWind * cleakage


# formula (**) calculate ppfVentSide
def cal_ppfVentSide(Cd, USide, ASide, vWind, AFlr, Cw):
    return Cd * USide * ASide * vWind * sqrt(Cw) / (2 * AFlr)


# formula 13, use formula 10 and formula (**)
def cal_fVentSide(nInsScr, ppfVentSide, fleakage, UThScr, ppfVentRoofSide, nSide, nSide_Thr):
    # nSide_Thr la nguong Stack
    # pp_fVentSide la f"VentSide tinh bang ppfVentRoofSide tai ARoof = 0
    if nSide >= nSide_Thr:
        return nInsScr * ppfVentSide + 0.5 * fleakage
    else:
        return nInsScr * (UThScr * ppfVentSide + (1 - UThScr) * ppfVentRoofSide * nSide) + 0.5 * fleakage


# formula 14
def cal_fVentForced(nInsScr, UVentForced, phiVentForced, AFlr):
    return nInsScr * UVentForced * phiVentForced / AFlr


# formula 15
def cal_MCTopOut(fVentRoof, CO2Top, CO2Out):
    return fVentRoof * (CO2Top - CO2Out)


# formula 16
def cal_fVentRoof(nInsScr, fleakage, UThScr, ppfVentRoofSide, nRoof, nSide, nRoof_Thr, ppfVentRoof):
    # nRoof_Thr la nguong Stack
    if nRoof >= nRoof_Thr:
        return nInsScr * ppfVentRoof + 0.5 * fleakage
    else:
        return nInsScr * (UThScr * ppfVentRoof + (1 - UThScr) * ppfVentRoofSide * nSide) + 0.5 * fleakage


# formula 17
def cal_ppfVentRoof(Cd, URoof, ARoof, AFlr, g, hVent, TAir, TOut, Cw, vWind):
    TMeanAir = (TAir + TOut) / 2
    part1 = Cd * URoof * ARoof / (2 * AFlr)
    part2 = g * hVent * (TAir - TOut) / 2 / TMeanAir + Cw * pow(vWind, 2)
    return part1 * sqrt(part2)


# formular 18 include 19
# tinh P
def cal_P(CO2Air, LAI):
    T_Can_K = 20 + 273
    J_POT = LAI * 210 * math.exp(37000 * (T_Can_K - 298.15) / (8.314 * T_Can_K * 298.15)) * (
                1 + math.exp((710 * 298.15 - 220000) / (8.314 * 298.15))) / (
                        1 + math.exp((710 * T_Can_K - 220000) / (8.314 * T_Can_K)))
    J = (J_POT + 38.5 - math.sqrt(math.pow(J_POT + 38.5, 2) - 2.8 * J_POT * 38.5)) / 1.4
    CO2Stom = 0.67 * CO2Air
    P = (J * (CO2Stom - 498.1)) / (4 * (CO2Stom + 2 * 498.1))
    return P


# tinh R
def cal_R(CO2Air, P):
    CO2Stom = 0.67 * CO2Air
    R = P * 498.1 / CO2Stom
    return R


def cal_MCAirCan(P, R, CBuf, CMaxBuf):
    MCH2O = 0.03
    hCBuf = 1
    if CBuf > CMaxBuf:
        hCBuf = 0
    return MCH2O * hCBuf * (P - R)

print("Input line of data: ")
i = int(input())
# Read data from excel file
data = pd.read_excel("data_CO2.xlsx")
df = pd.DataFrame(data)

nHeatCO2 = float(df.at[i, "nHeatCO2"])
UBlow = float(df.at[i, "UBlow"])
PBlow = float(df.at[i, 'PBlow'])
AFlr = float(df.at[i, 'AFlr'])
UExtCO2 = float(df.at[0, 'UExtCO2'])
phiExtCO2 = float(df.at[0, 'phiExtCO2'])
UPad = float(df.at[0, 'UPad'])
phiPad = float(df.at[0, 'phiPad'])
CO2Out = float(df.at[0, 'CO2Out'])
LAI = float(df.at[0, 'LAI'])
CBuf = float(df.at[0, 'CBuf'])
CMax_Buf = float(df.at[0, 'CMax_Buf'])
UThScr = float(df.at[0, 'UThScr'])
KThScr = float(df.at[0, 'KThScr'])
TAir = float(df.at[0, 'TAir'])
TTop = float(df.at[0, 'TTop'])
g = float(df.at[0, 'g'])
pAir = float(df.at[0, 'pAir'])
pTop = float(df.at[0, 'pTop'])
cleakage = float(df.at[0, 'cleakage'])
vWind = float(df.at[0, 'vWind'])
Cd = float(df.at[0, 'Cd'])
URoof = float(df.at[0, 'URoof'])
USide = float(df.at[0, 'USide'])
ARoof = float(df.at[0, 'ARoof'])
ASide = float(df.at[0, 'ASide'])
hSideRoof = float(df.at[0, 'hSideRoof'])
TOut = float(df.at[0, 'TOut'])
Cw = float(df.at[0, 'Cw'])
nSide = float(df.at[0, 'nSide'])
nSide_Thr = float(df.at[0, 'nSide_Thr'])
sInsScr = float(df.at[0, 'sInsScr'])
UVentForced = float(df.at[0, 'UVentForced'])
phiVentForced = float(df.at[0, 'phiVentForced'])
capCO2Air = float(df.at[0, 'capCO2Air'])
hVent = float(df.at[0, 'hVent'])
nRoof = float(df.at[0, 'nRoof'])
nRoof_Thr = float(df.at[0, 'nRoof_Thr'])
capCO2Top = float(df.at[0, 'capCO2Top'])
CO2Air = float(df.at[0, 'CO2Air'])
CO2Top = float(df.at[0, 'CO2Top'])

# formular 1
def dxCO2Air(CO2Air, CO2Top):
    # TODO

    ######## Calculate MCBlowAir ########
    MCBlowAir = cal_MCBlowAir(nHeatCO2, UBlow, PBlow, AFlr)
    # print("MCBlowAir = ", MCBlowAir)
    ######## Calculate MCExtAir ########
    MCExtAir = cal_MCExtAir(UExtCO2, phiExtCO2, AFlr)
    # print("MCExtAir = ", MCExtAir)
    ######## Calculate MCPadAir ########
    MCPadAir = cal_MCPadAir_2(UPad, phiPad, AFlr, CO2Out, CO2Air)
    # print("MCPadAir = ", MCPadAir)
    ######## Calculate MCAirCan ########
    P = cal_P(CO2Air, LAI)
    R = 0
    MCAirCan = cal_MCAirCan(P, R, CBuf, CMax_Buf)
    # print("MCAirCan = ", MCAirCan)
    ######## Calculate MCAirTop ########
    fThScr = cal_fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop)
    MCAirTop = cal_MCAirTop(fThScr, CO2Air, CO2Top)
    # print("MCAirTop = ", MCAirTop)
    ######## Calculte MCAirOut ########
    # Calculate fleakage
    fleakage = cal_fleakage(cleakage, vWind)

    # Calculate ppfVentRoofSide
    ppfVentRoofSide = cal_ppfVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind)

    # Calculate ppfVentSide
    ppfVentSide = cal_ppfVentSide(Cd, USide, ASide, vWind, AFlr, Cw)

    # Calculate fVentSide
    nInsScr = cal_nInsScr(sInsScr)
    fVentSide = cal_fVentSide(nInsScr, ppfVentSide, fleakage, UThScr, ppfVentRoofSide, nSide, nSide_Thr)

    # Calculate fVentForce
    fVentForced = float(cal_fVentForced(nInsScr, UVentForced, phiVentForced, AFlr))

    MCAirOut = cal_MCAirOut(fVentSide, fVentForced, CO2Air, CO2Out)
    # print("MCAirOut = ", MCAirOut)
    return (MCBlowAir + MCExtAir + MCPadAir - MCAirCan - MCAirTop - MCAirOut) / capCO2Air

# formula 2
def dxCO2Top(CO2Air, CO2Top):
    # TODO

    ######## Calculate MCAirTop ########
    fThScr = cal_fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop)
    MCAirTop = cal_MCAirTop(fThScr, CO2Air, CO2Top)
    # print("MCAirTop = ", MCAirTop)
    ######## Calculate MCTopOut ########
    # Calculate ppfVentRoofSide
    ppfVentRoofSide = cal_ppfVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind)

    # Calculate ppfVentRoof
    ppfVentRoof = cal_ppfVentRoof(Cd, URoof, ARoof, AFlr, g, hVent, TAir, TOut, Cw, vWind)

    # Calculate fleakage
    fleakage = cal_fleakage(cleakage, vWind)

    # Calculate fVentRoof
    nInsScr = cal_nInsScr(sInsScr)
    fVentRoof = cal_fVentRoof(nInsScr, fleakage, UThScr, ppfVentRoofSide, nRoof, nSide, nRoof_Thr, ppfVentRoof)
    MCTopOut = cal_MCTopOut(fVentRoof, CO2Top, CO2Out)
    # print("MCTopOut = ", MCTopOut)

    return (MCAirTop - MCTopOut) / capCO2Top

# print("dxCO2Air = ", dxCO2Air(CO2Air, CO2Top))
# print("dxCO2Top = ", dxCO2Top(CO2Air, CO2Top))
def euler(CO2Air, CO2Top, h, time):
    n = int(time/h)
    CO2Air_0 = CO2Air
    CO2Top_0 = CO2Top

    for idx in range(1, n + 1):
        k = h * dxCO2Air(CO2Air_0, CO2Top_0)
        t = h * dxCO2Top(CO2Air_0, CO2Top_0)

        CO2Air_0 += k
        CO2Top_0 += t

        if idx % 300 == 0:
            row = idx // 300
            worksheet.write(row, 0, idx)
            worksheet.write(row, 1, CO2Air_0)
            worksheet.write(row, 2, CO2Top_0)
            global TAir, TTop, TOut, vWind
            TAir = float(df.at[row, 'TAir'])
            TTop = float(df.at[row, 'TTop'])
            TOut = float(df.at[row, 'TOut'])
            vWind = float(df.at[row, 'vWind'])

    return CO2Air_0, CO2Top_0


def rk4(CO2Air, CO2Top, h, time):
    n = int(time/h)
    CO2Air_0 = CO2Air
    CO2Top_0 = CO2Top

    for idx in range(1, n + 1):
        k1 = h * dxCO2Air(CO2Air_0, CO2Top_0)
        t1 = h * dxCO2Top(CO2Air_0, CO2Top_0)
        k2 = h * dxCO2Air(CO2Air_0+0.5*k1, CO2Top_0+0.5*k1)
        t2 = h * dxCO2Top(CO2Air_0+0.5*t1, CO2Top_0+0.5*t1)
        k3 = h * dxCO2Air(CO2Air_0+0.5*k2, CO2Top_0+0.5*k2)
        t3 = h * dxCO2Top(CO2Air_0+0.5*t2, CO2Top_0+0.5*t2)
        k4 = h * dxCO2Air(CO2Air_0+k3, CO2Top_0+k3)
        t4 = h * dxCO2Top(CO2Air_0+t3, CO2Top_0+t3)

        CO2Air_0 += (1.0/6.0) * (k1+2*k2+2*k3+k4)
        CO2Top_0 += (1.0/6.0) * (t1+2*t2+2*t3+t4)

        if idx % 300 == 0:
            row = idx // 300
            worksheet.write(row, 4, CO2Air_0)
            worksheet.write(row, 5, CO2Top_0)
            global TAir, TTop, TOut, vWind
            TAir = float(df.at[row, 'TAir'])
            TTop = float(df.at[row, 'TTop'])
            TOut = float(df.at[row, 'TOut'])
            vWind = float(df.at[row, 'vWind'])

    return CO2Air_0, CO2Top_0


############## main ##############
def main():
    print(df.at[0, 'Place'])
    print('CO2Air_0: ', CO2Air)
    print('CO2Top_0: ', CO2Top)

    step = float(input('Input step: '))
    time = float(input('Input time: '))

    air, top = euler(CO2Air, CO2Top, step, time)
    print('\nExplicit Euler')
    print('The CO2Air: ', round(air, 10))
    print('The CO2Top: ', round(top, 10))

    air, top = rk4(CO2Air, CO2Top, step, time)
    print('\nExplicit Runge-Kutta 4th order')
    print('The CO2Air: ', round(air, 10))
    print('The CO2Top: ', round(top, 10))

    print(CO2Air)
    print(CO2Top)


wb = Workbook('Output.xlsx')
worksheet = wb.add_worksheet()
worksheet.write(0, 1, 'CO2Air')
worksheet.write(0, 2, 'CO2Top')
worksheet.write(0, 4, 'CO2Air')
worksheet.write(0, 5, 'CO2Top')
main()
wb.close()
