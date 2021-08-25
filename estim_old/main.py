import numpy as np
import leitura as lt
import estimation as estm
import math as m
import matplotlib.pyplot as plot
#Legend- Measurements
#1 = Active Injection
#2 = Reactive Injection
#3 = Active Power Flow
#4 = Reactive Power Flow
#5 = Current(Real)
#6 = Current(Imaginary)
#7 = Voltage Magnitude
#8 = Voltage Angle
# from | t | type | meas. value | std.dev
#______________________________________________________________________________________________________________________________________







def SE():
    
    # Read Parameters and Measurements  Files
    meas,top = lt.leitura('measurements_abur_ex2_2.txt','parametros_abur_ex2_2.txt')
    z= meas[:,3]
    std_dev = meas[:,4] 
    # Residual matrix
    num_meas,col2=meas.shape
    R = estm.res_matrix(std_dev,num_meas)
    # Ybus
    Ybus=estm.Ybus(top)
    # Main Estimator
    tol=1e-4
    lin,col=Ybus.shape
    num_states=2*lin-1
    # Initialization
    x_k = np.zeros(num_states) #first angle and after voltage module
    x_k[lin-1:2*lin-1] = 1
    teta = np.zeros(lin)
    V = np.ones(lin)
    iter = 0
    num_bus=lin


    #Iterative Newton Raphson Method
    while True:
        delta_x,h = estm.evalue_state(x_k,Ybus,std_dev,meas,num_states,num_meas,R,teta,V,num_bus,z)
        if(max(abs(delta_x)) < tol):
            break
        else:
            x_k =  x_k + delta_x     
            teta[1:lin] = x_k[1:lin] #refresh angle's values but not the reference bus
            V = x_k[lin-1:2*lin-1]#refresh voltage's values
    # Print Relatories



if __name__ == '__main__':
    SE()