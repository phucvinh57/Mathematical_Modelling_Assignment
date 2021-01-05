import math
from math import sqrt
import pandas as pd

# formula 1
def cal_MVCanAir(VECCanAir, VPCan, VPAir):
    return VECCanAir * (VPCan - VPAir)

# formula 2
def cal_VECCanAir(pAir, c_p_Air, LAI, delta_H, y, rb, rs): #y is gamma
    return 2 * pAir * c_p_Air * LAI / (delta_H * y *(rb + rs))

# formula 3
def cal_MVPadAir(pAir, fPad, nPad, xPad, xOut):
    return pAir * fPad *(nPad * (xPad - xOut) + xOut)
def cal_fPad(UPad, phiPad, AFlr): # use for formula 3 & formula 10
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
def cal_HECAirThScr(UThScr, TAir, TThScr): # use for formula 7
    return 1.7 * UThScr * pow(abs(TAir - TThScr), 0.33)

# formula 7
def cal_MVAirTop(MWater, R, fThScr, VPAir, VPTop, TAir, TTop):
    return (MWater / R) * fThScr * (VPAir / TAir - VPTop / TTop)
def cal_fThScr(UThScr, KThScr, TAir, TTop, g, pAir, pTop): # use for formula 6
    a = UThScr * KThScr * pow(abs(TAir - TTop), 2/3)
    PMean_Air = (pAir + pTop) / 2
    b = (1 - UThScr) * pow(g * (1 - UThScr) * abs(pAir - pTop) / (2 * PMean_Air), 1/2)
    return a + b

# formula 8
def cal_MVAirOut(MWater, R, fVentSide, fVentForced, VPAir, VPTop, TAir, TTop):
    return (MWater / R) * (fVentSide + fVentForced) * (VPAir / TAir - VPTop / TTop)
def cal_fVentSide(nInsScr, ppfVentSide, fleakage, UThScr, ppfVentRoofSide, nSide, nSide_Thr):  # use for formula 8
    if nSide >= nSide_Thr:
        return nInsScr * ppfVentSide + 0.5 * fleakage
    else:
        return nInsScr * (UThScr * ppfVentSide + (1 - UThScr) * ppfVentRoofSide * nSide) + 0.5 * fleakage
def cal_ppfVentSide(Cd, USide, ASide, vWind, AFlr, Cw):  # use for fVentSide
    return Cd * USide * ASide * vWind * sqrt(Cw) / (2 * AFlr)
def cal_ppfVentRoofSide(Cd, AFlr, URoof, USide, ARoof, ASide, g, hSideRoof, TAir, TOut, Cw, vWind): # use for fVentSide
    a = Cd / AFlr
    b = pow(URoof * USide * ARoof * ASide, 2) / (pow(URoof * ARoof, 2) + pow(USide * ASide, 2))
    TMean_Air = (TAir + TOut) / 2
    c = 2 * g * hSideRoof * (TAir - TOut) / TMean_Air
    _d = (URoof * ARoof + USide * ASide) / 2
    d = pow(_d, 2) * Cw * pow(vWind, 2)
    return a * sqrt(b * c + d)
def cal_fVentForced(nInsScr, UVentForced, phiVentForced, AFlr):  # use for formula 8
    return nInsScr * UVentForced * phiVentForced / AFlr

# formula 9
def cal_MVTopOut(MWater, R, fVentRoof, VPAir, VPTop, TAir, TTop):
    return (MWater / R) * fVentRoof * (VPAir / TAir - VPTop / TTop)
def cal_fVentRoof(nInsScr, fleakage, UThScr, ppfVentRoofSide, nRoof, nSide, nRoof_Thr, ppfVentRoof): # use for formula 9
    # nRoof_Thr la nguong Stack
    if nRoof >= nRoof_Thr:
        return nInsScr * ppfVentRoof + 0.5 * fleakage
    else:
        return nInsScr * (UThScr * ppfVentRoof + (1 - UThScr) * ppfVentRoofSide * nSide) + 0.5 * fleakage

# formula 10
def cal_MVAirOut_Pad(fPad, MWater, R, VPAir, TAir):
    return fPad * MWater / R * VPAir / TAir

# formula 11
def cal_MVAirMech(HECMechAir, VPAir, VPMech):
    if VPAir <= VPMech:
        return 0
    else:
        return 6.4 * pow(10, -9) * HECMechAir * (VPAir - VPMech)
def cal_HECMechAir(UMechCool, COPMechCool, PMechCool, AFlr, TAir, TMechCool, delta_H, VPAir, VPMechCool): #use for formula 11
    A1 = UMechCool * COPMechCool * PMechCool / AFlr
    A2 = TAir - TMechCool + 6.4 * pow(10, -9) * delta_H * (VPAir -VPMechCool)
    return A1 / A2

# formula 12
def cal_MVTopCov_in(HECTopCov_in, VPTop, VPCov_in):
    if VPTop <= VPCov_in:
        return 0
    else:
        return 6.4 * pow(10, -9) * HECTopCov_in * (VPTop - VPCov_in)
def cal_HECTopCov_in(cHECin, TTop, TCov_in, ACov, AFlr): #use for formula 12
    return cHECin * pow(TTop - TCov_in, 0.33) * ACov / AFlr