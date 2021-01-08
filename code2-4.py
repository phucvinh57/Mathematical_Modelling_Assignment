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


# # formula 19
# def cal_hCBuf(CBuf, CBufMax):
#     return (int)(CBufMax >= CBuf)

# # formula 22
# def get_abc(Res,CO2Air,CO2_05,PMax):
#     return Res,-(CO2Air+CO2_05+Res*PMax),CO2Air*PMax
# def cal_P(get_abc):
#     a,b,c = get_abc
#     return (-b-math.sqrt(b*b-4*a*c))/(2*a)

# # formula 23
# def cal_k(T, T0, kT0, Ha, R):
#     return kT0 * math.exp(-Ha / R * (1 / T - 1 / T0))

# # formula 24
# def cal_f(T, T0, Hd, R, S):
#     return (1 + math.exp(-Hd / R * (1 / T0 - S / Hd)))/(1 + math.exp(-Hd / R * (1 / T - S / Hd)))

# # formula 25
# def cal_PMax_T(k, f):
#     return k * f

# # formula 27
# def cal_L(L0, K, LAI, m):
#     return L0 * (1 - (K * math.exp(-K * LAI)) / (1 - m))

# # formula 28
# def cal_k_expand(LAI, k):
#     return LAI*k

# # formula 29
# # PMax_T tinh theo k_expand
# def cal_PMax_LT(P_MLT, PMax_T, L, L05):
#     return (P_MLT * PMax_T * L) / (L + L05)


# formular 1
def dxCO2Air(CO2Air, CO2Top):
    # TODO

    ######## Calculate MCBlowAir ########
    nHeatCO2 = nHeatCO2_0
    UBlow = UBlow_0
    PBlow = PBlow_0
    AFlr = AFlr_0
    MCBlowAir = cal_MCBlowAir(nHeatCO2, UBlow, PBlow, AFlr)

    ######## Calculate MCExtAir ########
    UExtCO2 = UExtCO2_0
    phiExtCO2 = phiExtCO2_0
    MCExtAir = cal_MCExtAir(UExtCO2, phiExtCO2, AFlr)

    ######## Calculate MCPadAir ########
    UPad = UPad_0
    phiPad = phiPad_0
    CO2Out = CO2Out_0
    MCPadAir = cal_MCPadAir_2(UPad, phiPad, AFlr, CO2Out, CO2Air)

    ######## Calculate MCAirCan ########
    LAI = LAI_0
    P = cal_P(CO2Air, LAI)
    R = 0
    CBuf = CBuf_0
    CMax_Buf = CMax_Buf_0
    MCAirCan = cal_MCAirCan(P, R, CBuf, CMax_Buf)

    ######## Calculate MCAirTop ########
    UThScr = UThScr_0
    KThScr = KThScr_0
    TAir = TAir_0
    TTop = TTop_0
    g = g_0
    pAir = pAir_0
    pTop = pTop_0
    fThScr = cal_fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop)
    MCAirTop = cal_MCAirTop(fThScr, CO2Air, CO2Top)

    ######## Calculte MCAirOut ########
    # Calculate fleakage
    cleakage = cleakage_0
    vWind = vWind_0
    fleakage = cal_fleakage(cleakage, vWind)

    # Calculate ppfVentRoofSide
    Cd = Cd_0
    URoof = URoof_0
    USide = USide_0
    ARoof = ARoof_0
    ASide = ASide_0
    hSideRoof = hSideRoof_0
    TOut = TOut_0
    Cw = Cw_0
    ppfVentRoofSide = cal_ppfVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind)

    # Calculate ppfVentSide
    ppfVentSide = cal_ppfVentSide(Cd, USide, ASide, vWind, AFlr, Cw)

    # Calculate fVentSide
    nSide = nSide_0
    nSide_Thr = nSide_Thr_0
    sInsScr = sInsScr_0
    nInsScr = cal_nInsScr(sInsScr)
    fVentSide = cal_fVentSide(nInsScr, ppfVentSide, fleakage, UThScr, ppfVentRoofSide, nSide, nSide_Thr)

    # Calculate fVentForce
    UVentForced = UVentForced_0
    phiVentForced = phiVentForced_0
    fVentForced = float(cal_fVentForced(nInsScr, UVentForced, phiVentForced, AFlr))

    MCAirOut = cal_MCAirOut(fVentSide, fVentForced, CO2Air, CO2Out)

    capCO2Air = capCO2Air_0
    return (MCBlowAir + MCExtAir + MCPadAir - MCAirCan - MCAirTop - MCAirOut) / capCO2Air


# formula 2
def dxCO2Top(CO2Air, CO2Top):
    # TODO

    ######## Calculate MCAirTop ########
    UThScr = UThScr_0
    KThScr = KThScr_0
    TAir = TAir_0
    TTop = TTop_0
    g = g_0
    pAir = pAir_0
    pTop = pTop_0
    fThScr = cal_fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop)
    MCAirTop = cal_MCAirTop(fThScr, CO2Air, CO2Top)

    ######## Calculate MCTopOut ########
    # Calculate ppfVentRoofSide
    AFlr = AFlr_0
    Cd = Cd_0
    URoof = URoof_0
    USide = USide_0
    ARoof = ARoof_0
    ASide = ASide_0
    hSideRoof = hSideRoof_0
    TOut = TOut_0
    Cw = Cw_0
    vWind = vWind_0
    ppfVentRoofSide = cal_ppfVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind)

    # Calculate ppfVentRoof
    hVent = hVent_0
    ppfVentRoof = cal_ppfVentRoof(Cd, URoof, ARoof, AFlr, g, hVent, TAir, TOut, Cw, vWind)

    # Calculate fleakage
    cleakage = cleakage_0
    fleakage = cal_fleakage(cleakage, vWind)

    # Calculate fVentRoof
    sInsScr = sInsScr_0
    nInsScr = cal_nInsScr(sInsScr)
    nRoof = nRoof_0
    nSide = nSide_0
    nRoof_Thr = nRoof_Thr_0
    fVentRoof = cal_fVentRoof(nInsScr, fleakage, UThScr, ppfVentRoofSide, nRoof, nSide, nRoof_Thr, ppfVentRoof)

    CO2Out = CO2Out_0
    MCTopOut = cal_MCTopOut(fVentRoof, CO2Top, CO2Out)

    capCO2Top = capCO2Top_0
    return (MCAirTop - MCTopOut) / capCO2Top


def euler(CO2Air_0, CO2Top_0, h, time):
    n = int(time/h)
    CO2Air = CO2Air_0
    CO2Top = CO2Top_0

    for idx in range(1, n + 1):
        k = h * dxCO2Air(CO2Air, CO2Top)
        t = h * dxCO2Top(CO2Air, CO2Top)

        CO2Air += k
        CO2Top += t

        if idx % 300 == 0:
            row = idx // 300
            worksheet.write(row, 0, idx)
            worksheet.write(row, 1, CO2Air)
            worksheet.write(row, 2, CO2Top)
            global TAir_0, TTop_0, TOut_0, vWind_0
            TAir_0 = float(df.at[row, 'TAir'])
            TTop_0 = float(df.at[row, 'TTop'])
            TOut_0 = float(df.at[row, 'TOut'])
            vWind_0 = float(df.at[row, 'vWind'])

    return CO2Air, CO2Top


def rk4(CO2Air_0, CO2Top_0, h, time):
    n = int(time/h)
    CO2Air = CO2Air_0
    CO2Top = CO2Top_0

    for idx in range(1, n + 1):
        k1 = h * dxCO2Air(CO2Air, CO2Top)
        t1 = h * dxCO2Top(CO2Air, CO2Top)
        k2 = h * dxCO2Air(CO2Air+0.5*k1, CO2Top+0.5*k1)
        t2 = h * dxCO2Top(CO2Air+0.5*t1, CO2Top+0.5*t1)
        k3 = h * dxCO2Air(CO2Air+0.5*k2, CO2Top+0.5*k2)
        t3 = h * dxCO2Top(CO2Air+0.5*t2, CO2Top+0.5*t2)
        k4 = h * dxCO2Air(CO2Air+k3, CO2Top+k3)
        t4 = h * dxCO2Top(CO2Air+t3, CO2Top+t3)

        CO2Air += (1.0/6.0) * (k1+2*k2+2*k3+k4)
        CO2Top += (1.0/6.0) * (t1+2*t2+2*t3+t4)

        if idx % 300 == 0:
            row = idx // 300
            worksheet.write(row, 4, CO2Air)
            worksheet.write(row, 5, CO2Top)
            global TAir_0, TTop_0, TOut_0, vWind_0
            TAir_0 = float(df.at[row, 'TAir'])
            TTop_0 = float(df.at[row, 'TTop'])
            TOut_0 = float(df.at[row, 'TOut'])
            vWind_0 = float(df.at[row, 'vWind'])

    return CO2Air, CO2Top


############## main ##############
def main():
    print(df.at[0, 'Place'])
    print('CO2Air_0: ', CO2Air_0)
    print('CO2Top_0: ', CO2Top_0)

    step = float(input('Input step: '))
    time = float(input('Input time: '))

    air, top = euler(CO2Air_0, CO2Top_0, step, time)
    print('\nExplicit Euler')
    print('The CO2Air: ', round(air, 10))
    print('The CO2Top: ', round(top, 10))

    global TAir_0, TTop_0, TOut_0, vWind_0
    TAir_0 = float(df.at[0, 'TAir'])
    TTop_0 = float(df.at[0, 'TTop'])
    TOut_0 = float(df.at[0, 'TOut'])
    vWind_0 = float(df.at[0, 'vWind'])

    air, top = rk4(CO2Air_0, CO2Top_0, step, time)
    print('\nExplicit Runge-Kutta 4th order')
    print('The CO2Air: ', round(air, 10))
    print('The CO2Top: ', round(top, 10))


# Read data from excel file
data = pd.read_excel("data.xlsx")
df = pd.DataFrame(data)
nHeatCO2_0 = float(df.at[0, "nHeatCO2"])
UBlow_0 = float(df.at[0, "UBlow"])
PBlow_0 = float(df.at[0, 'PBlow'])
AFlr_0 = float(df.at[0, 'AFlr'])
UExtCO2_0 = float(df.at[0, 'UExtCO2'])
phiExtCO2_0 = float(df.at[0, 'phiExtCO2'])
UPad_0 = float(df.at[0, 'UPad'])
phiPad_0 = float(df.at[0, 'phiPad'])
CO2Out_0 = float(df.at[0, 'CO2Out'])
LAI_0 = float(df.at[0, 'LAI'])
CBuf_0 = float(df.at[0, 'CBuf'])
CMax_Buf_0 = float(df.at[0, 'CMax_Buf'])
UThScr_0 = float(df.at[0, 'UThScr'])
KThScr_0 = float(df.at[0, 'KThScr'])
TAir_0 = float(df.at[0, 'TAir'])
TTop_0 = float(df.at[0, 'TTop'])
g_0 = float(df.at[0, 'g'])
pAir_0 = float(df.at[0, 'pAir'])
pTop_0 = float(df.at[0, 'pTop'])
cleakage_0 = float(df.at[0, 'cleakage'])
vWind_0 = float(df.at[0, 'vWind'])
Cd_0 = float(df.at[0, 'Cd'])
URoof_0 = float(df.at[0, 'URoof'])
USide_0 = float(df.at[0, 'USide'])
ARoof_0 = float(df.at[0, 'ARoof'])
ASide_0 = float(df.at[0, 'ASide'])
hSideRoof_0 = float(df.at[0, 'hSideRoof'])
TOut_0 = float(df.at[0, 'TOut'])
Cw_0 = float(df.at[0, 'Cw'])
nSide_0 = float(df.at[0, 'nSide'])
nSide_Thr_0 = float(df.at[0, 'nSide_Thr'])
sInsScr_0 = float(df.at[0, 'sInsScr'])
UVentForced_0 = float(df.at[0, 'UVentForced'])
phiVentForced_0 = float(df.at[0, 'phiVentForced'])
capCO2Air_0 = float(df.at[0, 'capCO2Air'])
hVent_0 = float(df.at[0, 'hVent'])
nRoof_0 = float(df.at[0, 'nRoof'])
nRoof_Thr_0 = float(df.at[0, 'nRoof_Thr'])
capCO2Top_0 = float(df.at[0, 'capCO2Top'])
CO2Air_0 = float(df.at[0, 'CO2Air'])
CO2Top_0 = float(df.at[0, 'CO2Top'])

wb = Workbook('Output.xlsx')
worksheet = wb.add_worksheet()
worksheet.write(0, 1, 'CO2Air')
worksheet.write(0, 2, 'CO2Top')
worksheet.write(0, 4, 'CO2Air')
worksheet.write(0, 5, 'CO2Top')
main()
wb.close()
