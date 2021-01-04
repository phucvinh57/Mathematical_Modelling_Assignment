import math
from math import sqrt
import pandas as pd
 
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
    a = UThScr * KThScr * pow(abs(TAir - TTop), 2/3)
    PMean_Air = (pAir + pTop) / 2
    b = (1 - UThScr) * pow(g * (1 - UThScr) * abs(pAir - pTop) / (2 * PMean_Air), 1/2)
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
    part2 = g * hVent * (TAir - TOut)  / 2 / TMeanAir + Cw * pow(vWind, 2)
    return part1 * sqrt(part2)

# formular 18 include 19
def cal_MCAirCan(P, R, CBuf, CMaxBuf):
    MCH2O = 30
    hCBuf = 1
    if CBuf > CMaxBuf:
        hCBuf = 0
    return MCH2O * hCBuf * (P - R)

# formula 19
def cal_hCBuf(CBuf, CBufMax):
    return (int)(CBufMax >= CBuf)

# formula 22
def get_abc(Res,CO2Air,CO2_05,PMax):
    return Res,-(CO2Air+CO2_05+Res*PMax),CO2Air*PMax
def cal_P(get_abc):
    a,b,c = get_abc
    return (-b-math.sqrt(b*b-4*a*c))/(2*a)

# formula 23
def cal_k(T, T0, kT0, Ha, R):
    return kT0 * math.exp(-Ha / R * (1 / T - 1 / T0))

# formula 24
def cal_f(T, T0, Hd, R, S):
    return (1 + math.exp(-Hd / R * (1 / T0 - S / Hd)))/(1 + math.exp(-Hd / R * (1 / T - S / Hd)))

# formula 25
def cal_PMax_T(k, f):
    return k * f

# formula 27
def cal_L(L0, K, LAI, m):
    return L0 * (1 - (K * math.exp(-K * LAI)) / (1 - m))

# formula 28
def cal_k_expand(LAI, k):
    return LAI*k

# formula 29 
# PMax_T tinh theo k_expand
def cal_PMax_LT(P_MLT, PMax_T, L, L05):
    return (P_MLT * PMax_T * L) / (L + L05)


# formular 1
def dxCO2Air(CO2Air, CO2Top, i):
    # Read data from excel file
    # TODO
    data = pd.read_excel("D:/Study/Mô hình hoá/Assignment/data.xlsx")
    df = pd.DataFrame(data)
    
    ######## Calculate MCBlowAir ########
    nHeatCO2 = float(df.at[i, "nHeatCO2"])
    UBlow = float(df.at[i, "UBlow"])
    PBlow = float(df.at[i, 'PBlow'])
    AFlr = float(df.at[i, 'AFlr'])
    MCBlowAir = cal_MCBlowAir(nHeatCO2, UBlow, PBlow, AFlr)
    print(MCBlowAir)
    
    ######## Calculate MCExtAir ########
    UExtCO2 = float(df.at[i, 'UExtCO2'])
    phiExtCO2 = float(df.at[i, 'phiExtCO2'])
    MCExtAir = cal_MCExtAir(UExtCO2, phiExtCO2, AFlr)
    print(MCExtAir)

    ######## Calculate MCPadAir ########
    UPad = float(df.at[i, 'UPad'])
    phiPad = float(df.at[i, 'phiPad'])
    CO2Out = float(df.at[i, 'CO2Out'])
    MCPadAir = cal_MCPadAir_2(UPad, phiPad, AFlr, CO2Out, CO2Air)
    print(MCPadAir)

    ######## Calculate MCAirCan ########
    P = 1
    R = 0
    CBuf = float(df.at[i, 'CBuf'])
    CMax_Buf = float(df.at[i, 'CMax_Buf'])
    MCAirCan = cal_MCAirCan(P, R, CBuf, CMax_Buf)
    print(MCAirCan)

    ######## Calculate MCAirTop ########
    UThScr = float(df.at[i, 'UThScr'])
    KThScr = float(df.at[i, 'KThScr'])
    TAir = float(df.at[i, 'TAir'])
    TTop = float(df.at[i, 'TTop'])
    g = float(df.at[i, 'g'])
    pAir = float(df.at[i, 'pAir'])
    pTop = float(df.at[i, 'pTop'])
    fThScr = cal_fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop)
    MCAirTop = cal_MCAirTop(fThScr, CO2Air, CO2Top)
    print(MCAirTop)

    ######## Calculte MCAirOut ########
    # Calculate fleakage
    cleakage = float(df.at[i, 'cleakage'])
    vWind = float(df.at[i, 'vWind'])
    fleakage = cal_fleakage(cleakage, vWind) 
    
    # Calculate ppfVentRoofSide
    Cd = float(df.at[i, 'Cd'])
    URoof = float(df.at[i, 'URoof'])
    USide = float(df.at[i, 'USide'])
    ARoof = float(df.at[i, 'ARoof'])
    ASide = float(df.at[i, 'ASide'])
    hSideRoof = float(df.at[i, 'hSideRoof'])
    TOut = float(df.at[i, 'TOut'])
    Cw = float(df.at[i, 'Cw'])
    ppfVentRoofSide = cal_ppfVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind)

    # Calculate ppfVentSide
    ppfVentSide = cal_ppfVentSide(Cd, USide, ASide, vWind, AFlr, Cw)

    # Calculate fVentSide
    nSide = float(df.at[i, 'nSide'])
    nSide_Thr = float(df.at[i, 'nSide_Thr'])
    sInsScr = float(df.at[i, 'sInsScr'])
    nInsScr = cal_nInsScr(sInsScr)
    fVentSide = cal_fVentSide(nInsScr, ppfVentSide, fleakage, UThScr, ppfVentRoofSide, nSide, nSide_Thr)

    # Calculate fVentForce
    UVentForced = float(df.at[i, 'UVentForced'])
    phiVentForced = float(df.at[i, 'phiVentForced'])
    fVentForced = float(cal_fVentForced(nInsScr, UVentForced, phiVentForced, AFlr))

    MCAirOut = cal_MCAirOut(fVentSide, fVentForced, CO2Air, CO2Out)
    print(MCAirOut)

    capCO2Air = float(df.at[i, 'capCO2Air'])
    return (MCBlowAir + MCExtAir + MCPadAir - MCAirCan - MCAirTop - MCAirOut) / capCO2Air


# formula 2
def dxCO2Top(CO2Air, CO2Top, i):
    # Read data from excel file
    # TODO
    ######## Calculate MCAirTop ########
    UThScr = float(df.at[i, 'UThScr'])
    KThScr = float(df.at[i, 'KThScr'])
    TAir = float(df.at[i, 'TAir'])
    TTop = float(df.at[i, 'TTop'])
    g = float(df.at[i, 'g'])
    pAir = float(df.at[i, 'pAir'])
    pTop = float(df.at[i, 'pTop'])
    fThScr = cal_fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop)
    MCAirTop = cal_MCAirTop(fThScr, CO2Air, CO2Top)
    print(MCAirTop)


    ######## Calculate MCTopOut ########
    # Calculate ppfVentRoofSide
    AFlr = float(df.at[i, 'AFlr'])
    Cd = float(df.at[i, 'Cd'])
    URoof = float(df.at[i, 'URoof'])
    USide = float(df.at[i, 'USide'])
    ARoof = float(df.at[i, 'ARoof'])
    ASide = float(df.at[i, 'ASide'])
    hSideRoof = float(df.at[i, 'hSideRoof'])
    TOut = float(df.at[i, 'TOut'])
    Cw = float(df.at[i, 'Cw'])
    vWind = float(df.at[i, 'vWind'])
    ppfVentRoofSide = cal_ppfVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind)
    
    # Calculate ppfVentRoof
    hVent = float(df.at[i, 'hVent'])
    ppfVentRoof = cal_ppfVentRoof(Cd, URoof, ARoof, AFlr, g, hVent, TAir, TOut, Cw, vWind)
    
    # Calculate fleakage
    cleakage = float(df.at[i, 'cleakage'])
    fleakage = cal_fleakage(cleakage, vWind) 
    
    # Calculate fVentRoof
    sInsScr = float(df.at[i, 'sInsScr'])
    nInsScr = cal_nInsScr(sInsScr)
    nRoof = float(df.at[i, 'nRoof'])
    nSide = float(df.at[i, 'nSide'])
    nRoof_Thr = float(df.at[i, 'nRoof_Thr'])
    fVentRoof = cal_fVentRoof(nInsScr, fleakage, UThScr, ppfVentRoofSide, nRoof, nSide, nRoof_Thr, ppfVentRoof)

    CO2Out = float(df.at[i, 'CO2Out'])
    MCTopOut = cal_MCTopOut(fVentRoof, CO2Top, CO2Out)
    print(MCTopOut)

    capCO2Top = float(df.at[i, 'capCO2Top'])
    return (MCAirTop - MCTopOut) / capCO2Top

############## main ##############
data = pd.read_excel("D:/Study/Mô hình hoá/Assignment/data.xlsx")
df = pd.DataFrame(data)
i = 1
print(df.at[i, 'Place'])
CO2Air = float(df.at[i, 'CO2Air'])
CO2Top = float(df.at[i, 'CO2Top'])
#dxCO2Air(CO2Air, CO2Top, i)
dxCO2Top(CO2Air, CO2Top, i)