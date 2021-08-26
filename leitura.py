import numpy as np

def leitura(measurements,topology,shunt_bus):
    meas = np.loadtxt(measurements)
    top = np.loadtxt(topology)
    shunt = np.loadtxt(shunt_bus)
    meas = sort_meas_type(meas)
        
    return meas,top,shunt

def topology_and_param(topology,shunt_bus):
    top = np.loadtxt(topology)
    shunt = np.loadtxt(shunt_bus)
    return top,shunt
    
def leitura_meas_database(m,y,c):
    # month=['jan', 'fev', 'mar', 'abril','maio','jun','jul','agost','set','out','nov','dez']
    year=['2019']
    # caso=['14Bus_Base']
    caso=['14Bus_CasoBase']
    # file= 'Medidas\int_5-5min\Medidas_ieee_' + m + '_' + year[y] + '_' + caso[c] + 'SE.txt'
    file= 'Medidas\int_5-5min_std_artigo\Medidas_ieee_' + m + '_' + year[y] + '_' + caso[c] + 'SESTD_artigo.txt'
    # file= 'Medidas\int_5-5min_unitario\Medidas_ieee_' + m + '_' + year[y] + '_' + caso[c] + 'SEunitarioSTD.txt'
    meas = np.loadtxt(file,dtype = np.float64)

def sort_meas_type(meas):
    
    dict_sort_meds={1:1,2:3,3:2,4:4,5:5,6:6,7:7,8:8}
    dict_sort_reverse={1:1,3:2,2:3,4:4,5:5,6:6,7:7,8:8}
    num_meas=len(meas[:,0])
    for i in range(num_meas):
        meas[i,2] = dict_sort_meds[meas[i,2]]
    #Ordena de Acordo com Abur Cap 2
    ind = np.lexsort((meas[:,1],meas[:,0],meas[:,2]))
    meas = meas[ind]
    # meas = meas[meas[:,2].argsort()] #ordena o arquivo das medias pela coluna 3 "tipo"
    #Reescreve as numeracoes Originais das Medidas para o programa de EE
    for k in range(num_meas):
        meas[k,2] = dict_sort_reverse[meas[k,2]]

    return meas