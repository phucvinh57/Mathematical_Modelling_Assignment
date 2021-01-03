import math
from math import sqrt
 
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
def cal_fThScr(UThScr, kThScr, tAir, tTop, g, PAir, PTop):
    a = UThScr * kThScr * pow(abs(tAir - tTop), 2/3)
    PMean_Air = (PAir + PTop) / 2
    b = (1 - UThScr) * pow(g * (1 - UThScr) * abs(PAir - PTop) / (2 * PMean_Air), 1/2)
    return a + b

# formula 9
def cal_MCAirOut(fVentSide, fVentForce, CO2Air, CO2Out):
    return (fVentSide + fVentForce) * (CO2Air - CO2Out)

# formula 10
def cal_fVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, TMean_Air, cW, vWind):
    a = Cd / AFlr
    b = pow(URoof * USide * ARoof * ASide, 2) / (pow(URoof * ARoof, 2) + pow(USide * ASide, 2))
    c = 2 * g * hSideRoof * (TAir - TOut) / TMean_Air
    _d = (URoof * ARoof + USide * ASide) / 2
    d = pow(_d, 2) * cW * pow(vWind, 2)
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

# formula 13
def cal_fVentSide(fleakage, UThScr, fVentRoofSide, nSide, nSide_Thr):
    # nSide_Thr la nguong Stack
    # pp_fVentSide la f"VentSide tinh bang fVentRoofSide tai ARoof = 0
    pp_fVentSide = cal_fVentRoofSide(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    nInsScr = cal_nInsScr(0)
    if nSide >= nSide_Thr:
        return nInsScr * pp_fVentSide + 0.5 * fleakage
    else:
        nInsScr * (UThScr * pp_fVentSide + (1 - UThScr) * fVentRoofSide * nSide) + 0.5 * fleakage

# formula 14
def cal_fVentForced(nInsScr, UVentForced, phiVentForced, AFlr):
    return nInsScr * UVentForced * phiVentForced / AFlr

# formula 15
def cal_MCTopOut(fVentRoof, CO2Top, CO2Out):
    return fVentRoof * (CO2Top - CO2Out)

# formula 16
def cal_fVentRoof(nInsScr, fleakage, UThScr, fVentRoofSide, nRoof, nSide, nRoof_Thr):
    # nRoof_Thr la nguong Stack
    ppfVentRoof = cal_ppfVentRoof(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    if nRoof >= nRoof_Thr:
        return nInsScr * ppfVentRoof + 0.5 * fleakage
    else:
        return nInsScr * (UThScr * ppfVentRoof + (1 - UThScr) * fVentRoofSide * nSide) + 0.5 * fleakage

# formula 17
def cal_ppfVentRoof(Cd, URoof, ARoof, AFlr, g, hRoof, TAir, TOut, TMeanAir, Cw, vWind):
    part1 = Cd * URoof * ARoof / 2 / AFlr
    part2 = g * hRoof * (TAir - TOut)  / 2 / TMeanAir + Cw * pow(vWind, 2)
    return part1 * sqrt(part2)

# formular 18 include 19
# tinh P
def cal_P(CO2Air,LAI):
    T_Can_K = 20+273
    J_POT = LAI*210*math.exp(37000*(T_Can_K-298.15)/(8.314*T_Can_K*298.15))*(1+math.exp((710*298.15-220000)/(8.314*298.15)))/(1+math.exp((710*T_Can_K-220000)/(8.314*T_Can_K)))
    J = (J_POT + 38.5 - math.sqrt(math.pow(J_POT+38.5, 2)-2.8*J_POT*38.5))/1.4
    CO2Stom = 0.67*CO2Air
    P = (J*(CO2Stom-498.1))/(4*(CO2Stom+2*498.1))
    return P
# tinh R
def cal_R(CO2Air,P):
    CO2Stom = 0.67*CO2Air
    R = P*498.1/CO2Stom
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
    # Read data from excel file
    # TODO
    
    return 0


# formula 2
def dxCO2Top(CO2Air, CO2Top):
    # Read data from excel file
    # TODO
    
    return 0


P = cal_P(440,2.5)
print(P)
R = cal_R(440,P)
print(R)
MCAirCan = cal_MCAirCan(P,R,0,1)
print(MCAirCan)