from Mathematical_Modelling_Assignment.code2b import CO2Air, cal_ppfVentRoof
import math
from math import sqrt
import pandas as pd

# formula 1
def cal_MVCanAir(VECCanAir, VPCan, VPAir):
    return VECCanAir * (VPCan - VPAir)

# formula 2
def cal_VECCanAir(pAir, c_p_Air, LAI, delta_H, y, rb, rs): #y is gamma
    return 2 * pAir * c_p_Air * LAI / (delta_H * y *(rb + rs))
def cal_rs(VPCan, VPAir, R_Can):
    rs_min = 82
    c_evap1 = 4.3
    c_evap2 = 0.54
    rf_RCan = (R_Can+c_evap1)/(R_Can+c_evap2)
    R_Can_SP = 5
    S_rs = 1/(1+math.exp(-1*(R_Can-R_Can_SP)))
    c_night_evap3 = 1.1*pow(10,-11)
    c_night_evap4 = 5.2*pow(10,-6)
    c_evap3 = c_night_evap3*(1-S_rs)+c_night_evap3*S_rs
    c_evap4 = c_night_evap4*(1-S_rs)+c_night_evap4*S_rs
    n_mg_ppm = 0.554
    rf_CO2Air = 1+c_evap3*math.pow(n_mg_ppm*CO2Air-200, 2)
    rf_VPCan_VPAir = 1+c_evap4*pow(VPCan - VPAir, 2)
    return rs_min*rf_RCan*rf_CO2Air*rf_VPCan_VPAir
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
def cal_nInsScr(sInsScr):
    return sInsScr * (2 - sInsScr)
def cal_fleakage(cleakage, vWind):
    if vWind < 0.25:
        return 0.25 * cleakage
    else:
        return vWind * cleakage

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


###### Calculate dx ######
def dxVPAir(VPAir, VPTop, i):
    # TODO
    data = pd.read_excel("data_VP.xlsx")
    df = pd.DataFrame(data)
    ######## Calculate MVCanAir ########
    rb = 275
    # calculate rs
    pAir = float(df.at[i, 'pAir'])
    LAI = float(df.at[i, 'LAI'])
    VPCan = float(df.at[i, 'VPCan'])
    R_Can = float(df.at[i, 'RCan'])
    delta_H = 2450000
    y = 65.8
    c_p_Air = 1000
    rs = cal_rs(VPCan, VPAir, R_Can)
    VECCanAir = cal_VECCanAir(pAir,c_p_Air,LAI,delta_H,y,rb,rs)
    MVCanAir = cal_MVCanAir(VECCanAir, VPCan, VPAir)
    ######## Calculate MVPadAir ########
    UPad = float(df.at[i, 'UPad'])
    phiPad = float(df.at[i, 'phiPad'])
    AFlr = float(df.at[i, 'AFlr'])
    fPad = UPad*phiPad/AFlr
    nPad = float(df.at[i, 'nPad'])
    xPad = float(df.at[i, 'xPad'])
    xOut = float(df.at[i, 'xOut'])
    MVPadAir = cal_MVPadAir(pAir,fPad,nPad,xPad,xOut)
    ######## Calculate MVFogAir ########
    UFog = float(df.at[i, 'UFog'])
    phiFog = 1.39
    MVFogAir = cal_MVFogAir(UFog, phiFog, AFlr)
    ######## Calculate MVBlowAir ########
    nHeatVap = 4.43*pow(10,-8)
    UBlow = float(df.at[i, 'UBlow'])
    PBlow = float(df.at[i, 'PBlow'])
    MVBlowair = cal_MVBlowAir(nHeatVap, UBlow, PBlow, AFlr)
    ######## Calculate MVAirThScr ########
    UThScr = float(df.at[i, 'UThScr'])
    TAir = float(df.at[i, 'TAir'])
    TThScr = float(df.at[i, 'TThScr'])
    HECAirThScr = 1.7*UThScr*pow(abs(TAir - TThScr), 0.33)
    VPThScr = 
    MVAirThScr = cal_MVAirThScr(HECAirThScr, VPAir, VPThScr)
    ######## Calculate MVAirTop ########
    MWater = 18
    R = 8.314*pow(10,3)
    KThScr = float(df.at[i, 'KThScr'])
    TOut = float(df.at[i, 'TOut'])
    p_Mean_Air = float(df.at[i, 'p_Mean_Air'])
    pOut = float(df.at[i, 'pOut'])
    g = float(df.at[i, 'g'])
    fThScr = UThScr*KThScr*pow(abs(TAir - TOut), 0.66)+(1-UThScr)/p_Mean_Air*pow((0.5*p_Mean_Air*(1-UThScr)*g*abs(pAir-pOut)),0.5)
    TAir = float(df.at[i, 'TAir'])
    TTop = float(df.at[i, 'TTop'])
    MVAirTop = cal_MVAirTop(MWater,R,fThScr, VPAir, VPTop,TAir,TTop)
    ######## Calculate MVAirOut ########
    # calculate fVentSide
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
    ppfVentSide = cal_ppfVentSide(Cd,USide,ASide,vWind,AFlr,Cw)
    cleakage = float(df.at[i, 'cleakage'])
    fleakage = cal_fleakage(cleakage, vWind) 
    ppfVentRoofSide = cal_ppfVentRoofSide(Cd,AFlr,URoof,USide,ARoof,ASide,g,hSideRoof,TAir,TOut,Cw,vWind)
    nSide = float(df.at[i, 'nSide'])
    nSide_Thr = float(df.at[i, 'nSide_Thr'])
    fVentSide = cal_fVentSide(nInsScr, ppfVentSide, fleakage, UThScr, ppfVentRoofSide, nSide, nSide_Thr)
    # calculate fVentForced
    UVentForced = float(df.at[i, 'UVentForced'])
    phiVentForced = float(df.at[i, 'phiVentForced'])
    fVentForced = cal_fVentForced(nInsScr, UVentForced, phiVentForced, AFlr)
    MVAirOut = cal_MVAirOut(MWater,R,fVentSide,fVentForced,VPAir,VPTop,TAir,TTop)
    ######## Calculate MVAirOut_Pad ########
    MVAirOut_Pad = cal_MVAirOut_Pad(fPad,MWater,R,VPAir,TAir)
    ######## Calculate MVAirMech ########
    HECMechAir = 
    VPMech = 
    MVAirMech = cal_MVAirMech(HECMechAir,VPAir,VPMech)
    capVPAir = float(df.at[i, 'capVPAir'])
    return (MVCanAir+MVPadAir+MVFogAir+MVBlowair-MVAirThScr-MVAirTop-MVAirOut-MVAirOut_Pad-MVAirMech)/capVPAir

def dxVPTop(VPAir, VPTop, i):
    # TODO
    data = pd.read_excel("data_VP.xlsx")
    df = pd.DataFrame(data)
    ######## Calculate MVAirTop ########
    MWater = 18
    R = 8.314*pow(10,3)
    KThScr = float(df.at[i, 'KThScr'])
    TOut = float(df.at[i, 'TOut'])
    p_Mean_Air = float(df.at[i, 'p_Mean_Air'])
    pOut = float(df.at[i, 'pOut'])
    g = float(df.at[i, 'g'])
    UThScr = float(df.at[i, 'UThScr'])
    TAir = float(df.at[i, 'TAir'])
    pAir = float(df.at[i, 'pAir'])
    fThScr = UThScr*KThScr*pow(abs(TAir - TOut), 0.66)+(1-UThScr)/p_Mean_Air*pow((0.5*p_Mean_Air*(1-UThScr)*g*abs(pAir-pOut)),0.5)
    TAir = float(df.at[i, 'TAir'])
    TTop = float(df.at[i, 'TTop'])
    MVAirTop = cal_MVAirTop(MWater,R,fThScr, VPAir, VPTop,TAir,TTop)
    ######## Calculate MVTopCov_in ########
    cHECin = float(df.at[i, 'cHECin'])
    TCov_in = 
    ACov = float(df.at[i, 'ACov'])
    AFlr = float(df.at[i, 'AFlr'])
    MVTopCov_in = cal_HECTopCov_in(cHECin,TTop,TCov_in,ACov,AFlr)
    ######## Calculate MVTopOut ########
    # calculate fVentRoof
    nRoof = float(df.at[i, 'nRoof'])
    nRoof_Thr = float(df.at[i, 'nRoof_Thr'])
    hVent = float(df.at[i, 'hVent'])
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
    nSide = float(df.at[i, 'nSide'])
    cleakage = float(df.at[i, 'cleakage'])
    fleakage = cal_fleakage(cleakage, vWind) 
    ppfVentRoofSide = cal_ppfVentRoofSide(Cd,AFlr,URoof,USide,ARoof,ASide,g,hSideRoof,TAir,TOut,Cw,vWind)
    ppfVentRoof = cal_ppfVentRoof(Cd,URoof,ARoof,AFlr,g,hVent,TAir,TOut,Cw,vWind)
    fVentRoof = cal_fVentRoof(nInsScr,fleakage,UThScr,ppfVentRoofSide,nRoof,nSide,nRoof_Thr,ppfVentRoof)
    MVTopOut = cal_MVTopOut(MWater,R,fVentRoof,VPAir,VPTop,TAir,TTop)
    capVPTop = float(df.at[i, 'capVPTop'])
    return (MVAirTop-MVTopCov_in-MVTopOut)/capVPTop