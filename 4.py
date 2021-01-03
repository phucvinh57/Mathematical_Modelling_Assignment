from code2b import dxCO2Air, dxCO2Top
def euler(CO2Air0, CO2Top0, h, time):
    jump_step = int(time/h)
    CO2Air = CO2Air0
    CO2Top = CO2Top0

    for i in range(1, jump_step+1):
        k = h * dxCO2Air(CO2Air, CO2Top)
        t = h * dxCO2Top(CO2Air, CO2Top)

        CO2Air += k
        CO2Top += t

    return CO2Air, CO2Top


def rk4(CO2Air0, CO2Top0, h, time):
    jump_step = int(time/h)
    CO2Air = CO2Air0
    CO2Top = CO2Top0

    for i in range(1, jump_step+1):
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

    return CO2Air, CO2Top
