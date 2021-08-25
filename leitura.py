import numpy as np

def leitura(measurements,topology,shunt_bus):
    meas = np.loadtxt(measurements)
    top = np.loadtxt(topology)
    shunt = np.loadtxt(shunt_bus)
    dict_sort_meds={1:1,2:3,3:2,4:4,5:5,6:6,7:7,8:8}
    dict_sort_reverse={1:1,3:2,2:3,4:4,5:5,6:6,7:7,8:8}
    num_meas=len(meas[:,0])
    for i in range(num_meas):
        meas[i,2] = dict_sort_meds[meas[i,2]]
    #Ordena de Acordo com Abur Cap 2
    meas = np.sort(meas.view('i8,i8,i8'), order=['f2'], axis=0).view(np.int)
    #Reescreve as numeracoes Originais das Medidas
    for k in range(num_meas):
        meas[i,2] = dict_sort_reverse[meas[i,2]]
        
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
    dict_sort_meds={1:1,2:3,3:2,4:4,5:5,6:6,7:7,8:8}
    dict_sort_reverse={1:1,3:2,2:3,4:4,5:5,6:6,7:7,8:8}
    num_meas=len(meas[:,0])
    for i in range(num_meas):
        meas[i,2] = dict_sort_meds[meas[i,2]]
    #Ordena de Acordo com Abur Cap 2
    meas = np.sort(meas.view('i8,i8,i8'), order=['f2'], axis=0).view(np.int)
    #Reescreve as numeracoes Originais das Medidas
    for k in range(num_meas):
        meas[i,2] = dict_sort_reverse[meas[i,2]]

    return meas