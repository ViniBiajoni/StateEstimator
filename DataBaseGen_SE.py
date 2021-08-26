import numpy as np
import leitura as lt
# import estimation as estm
import math as m
import matplotlib.pyplot as plt
from itertools import chain
import pandas as pd
from EE_functions import *

#-----------------------------------------------------------------State Estimation--------------------------------------------------------------
def SE():

    # Read Parameters and Measurements  Files
    #meas,top,shunt = lt.leitura('measurements_abur_ex2_2.txt','parametros_abur_ex2_2.txt','shunt_bus_test.txt')
    top,shunt = lt.topology_and_param('Parametros\IEEE14bus_cdf.txt','Parametros\IEEE14shunt_bus_cdf.txt') #the topology and measurement database is fixed
    # Ybus
    Ybus,shunt_line,t,brt= Y_bus(top,shunt)
    #Caso Escolhido
    data = pd.read_csv('meas_plan14base.txt', sep=" ", header=None)
    data.columns = ["tipo"]
    num_meas_base = len(data['tipo'].tolist())
    total_meas = int(num_meas_base*2)
    #Iterative Construction of Database
    # months=['jan', 'fev', 'mar', 'abril','maio','jun','jul','agost','set','out','nov','dez']
    months=['maio']
    #teste de 1 mês
    year=['2019']
    caso=['14Bus_Base']
    num_months= len(months)
    num_years= len(year)
    num_cases= len(caso)
    for c in range(num_cases):
        for y in range(num_years):
            for m in months:
                meas_month = lt.leitura_meas_database(m,y,c) #catch all the mensal intervals
                #DEFINIR O NUM MAX DE INSTANTES
                total_loops = int(len(meas_month[:,0])/total_meas) #Obtain the number of instants per month
                inicio=0
                fim=total_meas
                for i in range(total_loops):
                    #Obtain the meas for the instant i
                    meas= meas_month[inicio:fim,0:5] #do not use the colums 6 and 7
                    meas = lt.sort_meas_type(meas)
                    ##################################                
                    z= meas[:,3]
                    std_dev = meas[:,4] 
                    # Residual matrix
                    num_meas,col2=meas.shape
                    var=1 #1 indicates that the value is the variance, otherwise if var = 0 the value is the std
                    R = res_matrix(std_dev,num_meas,var)
                    # Main Estimator
                    tol=1e-6
                    numbus,col=Ybus.shape
                    num_states=2*numbus-1
                    # Initialization
                    x1=[0]*numbus
                    #x1[0]= algum valor -> caso o valor seja especificado na barra 1
                    x2=[1]*(2*numbus-1-(numbus-1))
                    x_k = x1+x2
                    teta =[0]*numbus #np.mat(np.zeros((numbus,1)))
                    V =[1]*numbus #np.mat(np.ones((numbus,1)))
                    Ji=[]
                    iter = 0
                    #Iterative Newton Raphson Method
                    while True:
                        iter = iter + 1
                        delta_x,h,H,J = evalue_state(x_k,Ybus,std_dev,meas,num_states,num_meas,R,teta,V,numbus,z,shunt_line,t,brt)
                        Ji.append(J)
                        delta_x=list(chain(*delta_x))
                        abs_delta =  [abs(ele) for ele in delta_x]
                        if(max(abs_delta) < tol) or iter>49:
                            break
                        else:          
                            for k in range(1,2*numbus):
                                x_k[k] = x_k[k] + delta_x[k-1]
                            teta = x_k[0:numbus]#refresh angle's values but not the reference bus
                            V = x_k[numbus:2*numbus]#refresh voltage's values
                    print("Módulo das Tensões \n")
                    print(V)
                    print('\n')
                    print("Ângulos(°) \n")
                    teta=[(teta[i]*180)/m.pi for i in range(numbus)]
                    print(teta)
                ###############Plotting J(x)################
                    # plt.plot(list(range(iter)), Ji)
                    # plt.xlabel('Iterations')
                    # plt.ylabel('J(x)')
                    # plt.yscale('log')
                    # plt.show()
                    # Print Relatories

                ##################Preparing Jacobian and h(x) to create a data_base######################
                    H_exp = (np.asarray(H)).flatten()# Makes the Jacobian vector
                    H_exp = H_exp.tolist()
                    h_exp = (np.asarray(h)).flatten()
                    h_exp = h_exp.tolist()
                    ########Especify the file name######## 
                    Hfile = 'H_file'+ month[m] + year[y] + caso[c] + '.txt'
                    hfile = 'h'+ month[m] + year[y] + caso[c] + '.txt'
                    #########Print the Data#######
                    textfile = open(Hfile, 'a') # a is used to "append" new lines and create a database of multiple H(x) and h(x)
                    with textfile as f:
                        for item in H_exp:
                            f.write("%s " % item)
                        f.write("\n")   
                    textfile.close()

                    textfile = open(hfile, 'a') # a is used to "append" new lines and create a database of multiple h(x)
                    with textfile as f:
                        for item in h_exp:
                            f.write("%s " % item)
                        f.write("\n")   
                    textfile.close()
                    #Put the cursor in the new group of meas to be evaluated by the SE
                    inicio = inicio+total_meas
                    fim = fim+total_meas


if __name__ == '__main__':
    SE() #Run the State Estimation 