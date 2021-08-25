import numpy as np

def leitura(measurements,topology,shunt_bus):
    meas = np.loadtxt(measurements)
    top = np.loadtxt(topology)
    shunt = np.loadtxt(shunt_bus)
    return meas,top,shunt

def topology_and_param(topology,shunt_bus):
    top = np.loadtxt(topology)
    shunt = np.loadtxt(shunt_bus)
    return top,shunt
    
def leitura_meas_database(m,y,c):
    month=['jan', 'fev', 'mar', 'abril','maio','jun','jul','agost','set','out','nov','dez']
    year=['2019']
    # caso=['14Bus_Base']
    caso=['14Bus_CasoBase']
    # file= 'Medidas\int_5-5min\Medidas_ieee_' + month[m] + '_' + year[y] + '_' + caso[c] + 'SE.txt'
    file= 'Medidas\int_5-5min\Medidas_ieee_' + month[m] + '_' + year[y] + '_' + caso[c] + 'SESTD_artigo.txt'
    # file= 'Medidas\int_5-5min\Medidas_ieee_' + month[m] + '_' + year[y] + '_' + caso[c] + 'SEunitarioSTD.txt'
    meas = np.loadtxt(file,dtype = np.float64)
    return meas