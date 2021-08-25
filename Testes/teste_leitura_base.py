import numpy as np

month=['jan', 'fev', 'mar', 'abril','maio','jun','jul','agost','set','out','nov','dez']
year=['2019']
caso=['14Bus_Base']
file= 'Medidas\int_5-5min\Medidas_ieee_' + month[0] + '_' + year[0] + '_' + caso[0] + 'SE.txt'
meas = np.loadtxt(file)
print(meas)