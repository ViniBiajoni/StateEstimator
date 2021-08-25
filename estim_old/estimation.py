# from main import R
import numpy as np
import math as m

def Ybus(top):#Construct Ybus

    numbus = int(np.amax(top))
    Ybarra = np.zeros((numbus,numbus), dtype=complex)
    Adj = np.zeros((numbus,numbus))
    (lin,col)= top.shape
    # fb= top[:,0] #from bus
    # tb= top[:,1] #to bus
    # r= top[:,2] #resistence
    # x= top[:,3] #reactance
    # z= r + x*1j
    # y=1/z
    # b= top[:,4] # Ground Admittance, B/2
    # b= b*1j 

    for i in range(lin):
        Adj[int(top[i,0]-1),int(top[i,1]-1)] = 1 #from-to
        Adj[int(top[i,1]-1),int(top[i,0]-1)] = 1 #to-from
        Ybarra[int(top[i,0]-1),int(top[i,1]-1)] = -1*(1/complex(top[i,2],top[i,3])) #+ complex(0,top[i,4]) SEM SHUNT POR ENQUANTO
        Ybarra[int(top[i,1]-1),int(top[i,0]-1)] = -1*(1/complex(top[i,2],top[i,3])) #+ complex(0,top[i,4])
        Ybarra[int(top[i,0]-1),int(top[i,0]-1)] = Ybarra[int(top[i,0]-1),int(top[i,0]-1)] + (1/complex(top[i,2],top[i,3]))
        Ybarra[int(top[i,1]-1),int(top[i,1]-1)] = Ybarra[int(top[i,1]-1),int(top[i,1]-1)] + (1/complex(top[i,2],top[i,3]))

    print(Adj)  
    print(Ybarra) 
    return Ybarra


def res_matrix(std_dev,num_meas): #Construct R matrix
    R = np.zeros((num_meas,num_meas))

    for i in range(num_meas):

        R[i,i]=(std_dev[i])**2


    return R

def medidas_h(x_new,meas,Ybus,num_meas,num_states,teta,V,num_bus):
    h=np.zeros(num_meas)
    for i in range(num_meas):
        de=int(meas[i,0])
        para=int(meas[i,1])
        
        if meas[i,2]==1:#Pi
      
            for k in range(num_bus):
                h[i]=h[i]+V[k]*((Ybus[de-1,k].real*m.cos(teta[de-1]-teta[k]) + Ybus[de-1,k].imag*m.sin(teta[de-1]-teta[k]) ))
            h[i]= V[de-1]*h[i]
        
        if meas[i,2]==2:#Qi
          
            for k in range(num_bus):
                h[i]=h[i]+V[k]*((Ybus[de-1,k].real*m.sin(teta[de-1]-teta[k]) - Ybus[de-1,k].imag*m.cos(teta[de-1]-teta[k]) ))
            h[i]= V[de-1]*h[i]        
       
        if meas[i,2]==3:#Pij
            #desconsiderando a susceptância shunt
            h[i]=(V[de-1]**2)*(-1*(Ybus[de-1,para-1]).real) - V[de-1]*V[para-1]*(-1*(Ybus[de-1,para-1]).real*m.cos(teta[de-1]-teta[para-1]) + -1*(Ybus[de-1,para-1]).imag*m.sin(teta[de-1]-teta[para-1]))
        
        if meas[i,2]==4:#Qij
            h[i]=-1*(V[de-1]**2)*(-1*(Ybus[de-1,para-1]).imag) - V[de-1]*V[para-1]*(-1*(Ybus[de-1,para-1]).real*m.sin(teta[de-1]-teta[para-1]) - (-1*(Ybus[de-1,para-1]).imag*m.cos(teta[de-1]-teta[para-1])))
       
        if meas[i,2]==5:#Iij(real)
            h[i]=m.sqrt((((-1*Ybus[de-1,para-1]).real)**2 + ((-1*Ybus[de-1,para-1]).imag)**2)*(V[de-1]**2 + V[para-1]**2 -2*V[de-1]*V[para-1]*m.cos(teta[de-1]-teta[para-1]))) 
        
        """ if meas[i,2]==6:#Iij(imag)
            de=meas[i,0]
            h[i]= V[de-1] """
        
        if meas[i,2]==7:#V(magnitude)
            h[i]= V[de-1] 
        
        if meas[i,2]==8:#V(angle) 
            h[i]= teta[de-1]        #considering bus 1 as reference bus (not necessary if exists PMUs but it is not implemented yet)                    
    
    return h


def jacobian(x_k,Ybus,meas,num_states,num_meas,teta,V,num_bus):        
    barras_ang=np.arange(2,num_bus+1,1)#angulos
    barras_voltage=np.arange(1,num_bus+1,1)#tensoes em modulo
    H=np.zeros((num_meas,num_states))
    for i in range(num_meas):
        tipo=meas[i,2]
        de= int(meas[i,0])
        para= int(meas[i,1])

        if tipo==1: #Pi
            col=0  #contabiliza coluna
            for t in range(np.size(barras_ang, axis = 0)):
                #Teta
                if de == barras_ang[t]:
                    
                    for k in range(num_bus):
                        H[i,col]= H[i,col] + V[de-1]*V[k]*((-1*Ybus[de-1,k].real*m.sin(teta[de-1]-teta[k]) + Ybus[de-1,k].imag*m.cos(teta[de-1]-teta[k])))
                    H[i,col]= H[i,col] -(V[de-1]**2)*(Ybus[de-1,de-1].imag)
                else:
                        H[i,col]= V[de-1]*V[barras_ang[t]-1]*((Ybus[de-1,barras_ang[t]-1].real*m.sin(teta[de-1]-teta[barras_ang[t]-1]) - Ybus[de-1,barras_ang[t]-1].imag*m.cos(teta[de-1]-teta[barras_ang[t]-1])))
                col=col+1

            for v in range(np.size(barras_voltage, axis = 0)):
                 #V
                if de == barras_voltage[v]:
                    for k in range(num_bus):

                        H[i,col]= H[i,col] + V[k]*((Ybus[de-1,k].real*m.cos(teta[de-1]-teta[k]) + Ybus[de-1,k].imag*m.sin(teta[de-1]-teta[k])))     
                    H[i,col]= H[i,col] + V[de-1]*(Ybus[de-1,de-1].real)
                else:
                    H[i,col] = V[de-1]*((Ybus[de-1,barras_voltage[v]-1].real*m.cos(teta[de-1]-teta[barras_voltage[v]-1]) + Ybus[de-1,barras_voltage[v]-1].imag*m.sin(teta[de-1]-teta[barras_voltage[v]-1])))
                
                col=col+1
            
        if tipo==2: #Qi
            col=0  #contabiliza coluna
            
            for t in range(np.size(barras_ang, axis = 0)):
                #Teta
                if de == barras_ang[t]:
                    
                    for k in range(num_bus):
                        H[i,col]= H[i,col] + V[de-1]*V[k]*((Ybus[de-1,k].real*m.cos(teta[de-1]-teta[k]) + Ybus[de-1,k].imag*m.sin(teta[de-1]-teta[k]))) 
                    
                    H[i,col]= H[i,col] - (V[de-1]**2)*(Ybus[de-1,de-1].real)
                else:
                        H[i,col]= V[de-1]*V[barras_ang[t]-1]*((-1*Ybus[de-1,barras_ang[t]-1].real*m.cos(teta[de-1]-teta[barras_ang[t]-1]) - Ybus[de-1,barras_ang[t]-1].imag*m.sin(teta[de-1]-teta[barras_ang[t]-1])))
                col=col+1

            for v in range(np.size(barras_voltage, axis = 0)):
                 #V
                if de == barras_voltage[v]:
                    for k in range(num_bus):

                        H[i,col]= H[i,col] + V[k]*((Ybus[de-1,k].real*m.sin(teta[de-1]-teta[k]) - Ybus[de-1,k].imag*m.cos(teta[de-1]-teta[k])))     
                   
                    H[i,col]= H[i,col] - V[de-1]*(Ybus[de-1,de-1].imag)
                else:
                    H[i,col] = V[de-1]*((Ybus[de-1,barras_voltage[v]-1].real*m.sin(teta[de-1]-teta[barras_voltage[v]-1]) - Ybus[de-1,barras_voltage[v]-1].imag*m.cos(teta[de-1]-teta[barras_voltage[v]-1])))
                
                col=col+1

        if tipo==3: #Pij
            #teste = np.size(barras_ang, axis = 0)
            #print(teste)
            col=0  #contabiliza coluna
            for t in range(np.size(barras_ang, axis = 0)):
                #Teta
                if de == barras_ang[t]:
        
                    H[i,col]= V[de-1]*V[para-1]*((-1*Ybus[de-1,para-1].real)*m.sin(teta[de-1]-teta[para-1]) - (-1*Ybus[de-1,para-1].imag)*m.cos(teta[de-1]-teta[para-1]))

                elif para == barras_ang[t]:

                    H[i,col]= -1*V[de-1]*V[para-1]*((-1*Ybus[de-1,para-1].real)*m.sin(teta[de-1]-teta[para-1]) - (-1*Ybus[de-1,para-1].imag)*m.cos(teta[de-1]-teta[para-1]))
                
                else:
                    H[i,col]= 0
                
                col=col+1

            for v in range(np.size(barras_voltage, axis = 0)):
                 #V
                if de == barras_voltage[v]:
        
                    H[i,col]= -1*V[para-1]*((-1*Ybus[de-1,para-1].real)*m.cos(teta[de-1]-teta[para-1]) + (-1*Ybus[de-1,para-1].imag)*m.sin(teta[de-1]-teta[para-1])) + 2*(-1*Ybus[de-1,para-1].real)*V[de-1]   
                
                elif para == barras_voltage[v]:
        
                    H[i,col] = -1*V[de-1]*((-1*Ybus[de-1,para-1].real)*m.cos(teta[de-1]-teta[para-1]) + (-1*Ybus[de-1,para-1].imag)*m.sin(teta[de-1]-teta[para-1]))
                
                else:
                    H[i,col]= 0

                col=col+1

        if tipo==4: #Qij
            
            col=0  #contabiliza coluna
            
            for t in range(np.size(barras_ang, axis = 0)):
                #Teta
                if de == barras_ang[t]:
        
                    H[i,col]= -1*V[de-1]*V[para-1]*((-1*Ybus[de-1,para-1].real)*m.cos(teta[de-1]-teta[para-1]) + (-1*Ybus[de-1,para-1].imag)*m.sin(teta[de-1]-teta[para-1]))
        
                elif para == barras_ang[t]:
        
                    H[i,col]= V[de-1]*V[k]*((-1*Ybus[de-1,para-1].real)*m.cos(teta[de-1]-teta[para-1]) + (-1*Ybus[de-1,para-1].imag)*m.sin(teta[de-1]-teta[para-1]))
                
                else:

                    H[i,col]=0

                col=col+1

            for v in range(np.size(barras_voltage, axis = 0)):
                 #V
                if de == barras_voltage[v]:
        
                    H[i,col]= -1*V[para-1]*((-1*Ybus[de-1,para-1].real)*m.sin(teta[de-1]-teta[para-1]) - (-1*Ybus[de-1,para-1].imag)*m.cos(teta[de-1]-teta[para-1])) - 2*(-1*Ybus[de-1,para-1].imag)*V[de-1]   
                
                elif para == barras_voltage[v]:
        
                    H[i,col] = -1*V[de-1]*((-1*Ybus[de-1,para-1].real)*m.sin(teta[de-1]-teta[para-1]) - (-1*Ybus[de-1,para-1].imag)*m.cos(teta[de-1]-teta[para-1]))
                
                else:

                    H[i,col]=0

                col=col+1

        if tipo==5: #Iij
            col=0  #contabiliza coluna
            
            for t in range(np.size(barras_ang, axis = 0)):
                #Teta
                if de == barras_ang[t]:
        
                    H[i,col]=  ((((Ybus[de-1,].real)**2) + ((Ybus[de-1,].imag)**2))/(m.sqrt((((Ybus[de-1]).real)**2 + ((Ybus[de-1]).imag)**2)*(V[de-1]**2 + V[para-1]**2 -2*V[de-1]*V[para-1]*m.cos(teta[de-1]-teta[para-1])))))*V[de-1]*V[para-1]*m.sin(teta[de-1]-teta[para-1])
        
                elif para == barras_ang[t]:
                    
                    H[i,col]= -1*((((Ybus[de-1,].real)**2) + ((Ybus[de-1,].imag)**2))/(m.sqrt((((Ybus[de-1]).real)**2 + ((Ybus[de-1]).imag)**2)*(V[de-1]**2 + V[para-1]**2 -2*V[de-1]*V[para-1]*m.cos(teta[de-1]-teta[para-1])))))*V[de-1]*V[para-1]*m.sin(teta[de-1]-teta[para-1])

                else:

                    H[i,col]= 0
                
                col=col+1

            for v in range(np.size(barras_voltage, axis = 0)):
                 #V
                if de == barras_voltage[v]:
        
                    H[i,col]= ((((Ybus[de-1,].real)**2) + ((Ybus[de-1,].imag)**2))/(m.sqrt((((Ybus[de-1]).real)**2 + ((Ybus[de-1]).imag)**2)*(V[de-1]**2 + V[para-1]**2 -2*V[de-1]*V[para-1]*m.cos(teta[de-1]-teta[para-1])))))*(V[de-1]-V[para-1]*m.cos(teta[de-1]-teta[para-1]))
                
                elif para == barras_voltage[v]:

                    H[i,col]= ((((Ybus[de-1,].real)**2) + ((Ybus[de-1,].imag)**2))/(m.sqrt((((Ybus[de-1]).real)**2 + ((Ybus[de-1]).imag)**2)*(V[de-1]**2 + V[para-1]**2 -2*V[de-1]*V[para-1]*m.cos(teta[de-1]-teta[para-1])))))*(V[para-1]-V[de-1]*m.cos(teta[de-1]-teta[para-1]))
                
                else:
        
                    H[i,col] = 0
                
                col=col+1


        if tipo==7: #Vi
            
            col=0  #contabiliza coluna
            
            for t in range(np.size(barras_ang, axis = 0)):
                #Teta
                if de == barras_ang[t]:
        
                    H[i,col]= 0
        
                else:
        
                    H[i,col]= 0 
                
                col=col+1

            for v in range(np.size(barras_voltage, axis = 0)):
                 #V
                if de == barras_voltage[v]:
        
                    H[i,col]= 1  
                
                else:
        
                    H[i,col] = 0 
                
                col=col+1


    return H

def  evalue_state(x_k,Ybus,std_dev,meas,num_states,num_meas,R,teta,V,num_bus,z):
    x_new = np.zeros(num_states)
    H =jacobian(x_k,Ybus,meas,num_states,num_meas,teta,V,num_bus)
    h= medidas_h(x_new,meas,Ybus,num_meas,num_states,teta,V,num_bus)
    G = ((H.transpose()).dot((np.linalg.inv(R)))).dot(H)
    delta_x= (((np.linalg.inv(G)).dot((H.transpose()))).dot((np.linalg.inv(R)))).dot(np.subtract(z,h))
    return  delta_x,h

