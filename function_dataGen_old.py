def Y_bus(top,shunt_bus):#Construct Ybus

    numbus = int(np.amax(top))
    Ybarra = np.zeros((numbus,numbus), dtype=complex)
    Adj = np.zeros((numbus,numbus))
    (nbranch,col)= top.shape
    fb= top[:,0] #from bus
    tb= top[:,1] #to bus
    r= top[:,2] #resistence
    x= top[:,3] #reactance
    tap = top[:,5] # Tap value
    z= r + x*1j
    y=1/z
    b= top[:,4] 
    b= (b*1j)/2  # Ground Admittance, B/2
    shunt_line= np.zeros((numbus,numbus))
    t= np.ones((numbus,numbus)) #Tap matrix
    Bin_ref_tap = np.zeros((numbus,numbus)) # bin matrix to determine the "from" buses

    for i in range(nbranch): 
        Adj[int(top[i,0]-1),int(top[i,1]-1)] = 1 #from-to
        Adj[int(top[i,1]-1),int(top[i,0]-1)] = 1 #to-from
        
        shunt_line[int(fb[i]-1),int(tb[i]-1)] = (b[i]*1j)*-1
        shunt_line[int(tb[i]-1),int(fb[i]-1)] = shunt_line[int(fb[i]-1),int(tb[i]-1)]

        t[int(fb[i]-1),int(tb[i]-1)] = tap[i]
        t[int(tb[i]-1),int(fb[i]-1)] = t[int(fb[i]-1),int(tb[i]-1)]

        Bin_ref_tap[int(fb[i]-1),int(tb[i]-1)] = 1

    # Off Diagonal Elements    
    for k in range(nbranch):    
        # Ybarra[int(top[i,0]-1),int(top[i,1]-1)] = -1*(1/complex(top[i,2],top[i,3])) #+ complex(0,top[i,4]) SEM SHUNT POR ENQUANTO
        # Ybarra[int(top[i,1]-1),int(top[i,0]-1)] = -1*(1/complex(top[i,2],top[i,3])) #+ complex(0,top[i,4])
        # Ybarra[int(top[i,0]-1),int(top[i,0]-1)] = Ybarra[int(top[i,0]-1),int(top[i,0]-1)] + (1/complex(top[i,2],top[i,3]))
        # Ybarra[int(top[i,1]-1),int(top[i,1]-1)] = Ybarra[int(top[i,1]-1),int(top[i,1]-1)] + (1/complex(top[i,2],top[i,3]))
        Ybarra[int(fb[k]-1),int(tb[k]-1)]= Ybarra[int(fb[k]-1),int(tb[k]-1)] - (y[k]/tap[k])
        Ybarra[int(tb[k]-1),int(fb[k]-1)]= Ybarra[int(fb[k]-1),int(tb[k]-1)]

    #Diagonal Elements
    for m in range(numbus):
        for n in range(nbranch):
            if fb[n]-1 == m:
                Ybarra[m,m]= Ybarra[m,m] + (y[n]/(tap[n]**2)) + b[n]
            elif tb[n]-1 == m:
                Ybarra[m,m]= Ybarra[m,m] + y[n] + b[n]
    
    #Adding Shunt Bus in the Ybus
    for i in range(numbus):
        if shunt_bus[i]!=0:
            Ybarra[i,i]= Ybarra[i,i] + shunt_bus[i]*1j
    # print(Adj)  
    # print(Ybarra) 
    return Ybarra,shunt_line,t,Bin_ref_tap

#----------------------------------------------------------------------------------------------------------------------------------
def res_matrix(std_dev,num_meas,var): #Construct R matrix

    R = np.mat(np.zeros((num_meas,num_meas)))
    if var ==1:
        for i in range(num_meas):
            R[i,i] = std_dev[i]
    
    if var ==0:
        for i in range(num_meas):
            R[i,i] = (std_dev[i])**2


    return R
#-----------------------------------------------------------------------------------------------------------------
def lado_medida(t,brt,G,B,B_shunt,de,para): 
    tap=t.copy()
    if brt[de,para] == 1: 
        tap = 1/t[de,para]
        # g=(1-(1/t[de,para]))*G[de,para]
        # b=(1-(1/t[de,para]))*B[de,para] + B_shunt[de,para]
    else:
        tap = t[de,para]
        # g=(1-t[de,para])*G[de,para]
        # b=(1-t[de,para])*B[de,para] + B_shunt[de,para]

    return tap     
#----------------------------------------------------------------------------------------------------------------------------------
def medidas_h(x_new,meas,Ybus,num_meas,num_states,teta,V,num_bus,shunt_bus,t,brt):
    h=np.mat(np.zeros((num_meas,1)))
    #Measurement Function
    G= Ybus.real
    B= Ybus.imag
    #Here  the shunt bus is related with the lines (pi model)
    for i in range(num_meas):
        de=int(meas[i,0])
        para=int(meas[i,1])
        
        if meas[i,2]==1:#Pi
      
            for k in range(num_bus):
                h[i]=h[i]+V[k]*((G[de-1,k]*np.cos(teta[de-1]-teta[k]) + B[de-1,k]*np.sin(teta[de-1]-teta[k]) ))
            h[i]= V[de-1]*h[i]
        
        if meas[i,2]==2:#Qi
          
            for k in range(num_bus):
                h[i]=h[i]+V[k]*((G[de-1,k]*np.sin(teta[de-1]-teta[k]) - B[de-1,k]*np.cos(teta[de-1]-teta[k]) ))
            h[i]= V[de-1]*h[i]        
       
        if meas[i,2]==3:#Pij
            #desconsiderando a susceptância shunt
            tap=lado_medida(t,brt,G,B,shunt_bus,de-1,para-1)# é o tik
            h[i]= -(V[de-1]**2)*(tap*G[de-1,para-1]) - V[de-1]*V[para-1]*(-G[de-1,para-1]*np.cos(teta[de-1]-teta[para-1]) -B[de-1,para-1]*np.sin(teta[de-1]-teta[para-1]))
        
        if meas[i,2]==4:#Qij
            tap=lado_medida(t,brt,G,B,shunt_bus,de-1,para-1)# é o tik
            h[i]=-(V[de-1]**2)*(-tap*B[de-1,para-1] + shunt_bus[de-1,para-1]) - V[de-1]*V[para-1]*(-G[de-1,para-1]*np.sin(teta[de-1]-teta[para-1]) + B[de-1,para-1]*np.cos(teta[de-1]-teta[para-1]))
            
        if meas[i,2]==5:#Iij(real)
            #h[i]=m.sqrt((((-1*Ybus[de-1,para-1]).real)**2 + ((-1*Ybus[de-1,para-1]).imag)**2)*(V[de-1]**2 + V[para-1]**2 -2*V[de-1]*V[para-1]*m.cos(teta[de-1]-teta[para-1]))) 
            h[i]=np.sqrt(((-G[de-1,para-1])**2 + (-B[de-1,para-1])**2)*(V[de-1]**2 + V[para-1]**2 -2*V[de-1]*V[para-1]*np.cos(teta[de-1]-teta[para-1])))
        """ if meas[i,2]==6:#Iij(imag)
            de=meas[i,0]
            h[i]= V[de-1] """
        
        if meas[i,2]==7:#V(magnitude)
            h[i]= V[de-1] 
        
        if meas[i,2]==8:#V(angle) Randian
            h[i]= teta[de-1]        #considering bus 1 as reference bus (not necessary if exists PMUs but it is not implemented yet)                    
    
    return h

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def jacobian(x_k,Ybus,meas,num_states,num_meas,teta,V,num_bus,shunt_bus,tap,brt):        
    barras_ang=np.arange(2,num_bus+1,1)#angulos
    barras_voltage=np.arange(1,num_bus+1,1)#tensoes em modulo
    H=np.zeros((num_meas,num_states))
    H=np.mat(np.zeros((num_meas,num_states)))
    G= Ybus.real
    B= Ybus.imag
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
                        H[i,col]= H[i,col] + V[de-1]*V[k]*(-G[de-1,k]*np.sin(teta[de-1]-teta[k]) + B[de-1,k]*np.cos(teta[de-1]-teta[k]))
                    H[i,col]= H[i,col] -(V[de-1]**2)*B[de-1,de-1]
                else:
                        H[i,col]= V[de-1]*V[barras_ang[t]-1]*((G[de-1,barras_ang[t]-1]*np.sin(teta[de-1]-teta[barras_ang[t]-1]) - B[de-1,barras_ang[t]-1]*np.cos(teta[de-1]-teta[barras_ang[t]-1])))
                col=col+1

            for v in range(np.size(barras_voltage, axis = 0)):
                 #V
                if de == barras_voltage[v]:
                    for k in range(num_bus):

                        H[i,col]= H[i,col] + V[k]*((G[de-1,k]*np.cos(teta[de-1]-teta[k]) + B[de-1,k]*np.sin(teta[de-1]-teta[k])))     
                    H[i,col]= H[i,col] + V[de-1]*G[de-1,de-1]
                else:
                    H[i,col] = V[de-1]*((G[de-1,barras_voltage[v]-1]*np.cos(teta[de-1]-teta[barras_voltage[v]-1]) + B[de-1,barras_voltage[v]-1]*np.sin(teta[de-1]-teta[barras_voltage[v]-1])))
                
                col=col+1
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
        if tipo==2: #Qi
            col=0  #contabiliza coluna
            
            for t in range(np.size(barras_ang, axis = 0)):
                #Teta
                if de == barras_ang[t]:
                    
                    for k in range(num_bus):
                        H[i,col]= H[i,col] + V[de-1]*V[k]*((G[de-1,k]*np.cos(teta[de-1]-teta[k]) + B[de-1,k]*np.sin(teta[de-1]-teta[k]))) 
                    
                    H[i,col]= H[i,col] - (V[de-1]**2)*G[de-1,de-1]
                else:
                        H[i,col]= V[de-1]*V[barras_ang[t]-1]*((-G[de-1,barras_ang[t]-1]*np.cos(teta[de-1]-teta[barras_ang[t]-1]) - B[de-1,barras_ang[t]-1]*np.sin(teta[de-1]-teta[barras_ang[t]-1])))
                col=col+1

            for v in range(np.size(barras_voltage, axis = 0)):
                 #V
                if de == barras_voltage[v]:
                    for k in range(num_bus):

                        H[i,col]= H[i,col] + V[k]*((G[de-1,k]*np.sin(teta[de-1]-teta[k]) - B[de-1,k]*np.cos(teta[de-1]-teta[k])))     
                   
                    H[i,col]= H[i,col] - V[de-1]*B[de-1,de-1]
                else:
                    H[i,col] = V[de-1]*((G[de-1,barras_voltage[v]-1]*np.sin(teta[de-1]-teta[barras_voltage[v]-1]) - B[de-1,barras_voltage[v]-1]*np.cos(teta[de-1]-teta[barras_voltage[v]-1])))
                
                col=col+1
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if tipo==3: #Pij
            #teste = np.size(barras_ang, axis = 0)
            #print(teste)
            col=0  #contabiliza coluna
            for t in range(np.size(barras_ang, axis = 0)):
                #Teta
                if de == barras_ang[t]:
        
                    H[i,col]= V[de-1]*V[para-1]*(-G[de-1,para-1]*np.sin(teta[de-1]-teta[para-1]) + B[de-1,para-1]*np.cos(teta[de-1]-teta[para-1]))

                elif para == barras_ang[t]:

                    H[i,col]= -V[de-1]*V[para-1]*(-G[de-1,para-1]*np.sin(teta[de-1]-teta[para-1]) + B[de-1,para-1]*np.cos(teta[de-1]-teta[para-1]))
                
                else:
                    H[i,col]= 0
                
                col=col+1

            for v in range(np.size(barras_voltage, axis = 0)):
                 #V
                if de == barras_voltage[v]:

                    tap=lado_medida(t,brt,G,B,shunt_bus,de-1,para-1)# é o tik
                    H[i,col]= -V[para-1]*(-G[de-1,para-1]*np.cos(teta[de-1]-teta[para-1]) - B[de-1,para-1]*np.sin(teta[de-1]-teta[para-1])) - 2*(tap*G[de-1,para-1])*V[de-1]   
                
                elif para == barras_voltage[v]:
        
                    H[i,col] = -V[de-1]*(-G[de-1,para-1]*np.cos(teta[de-1]-teta[para-1]) - B[de-1,para-1]*np.sin(teta[de-1]-teta[para-1]))
                
                else:
                    H[i,col]= 0

                col=col+1
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if tipo==4: #Qij
            
            col=0  #contabiliza coluna
            
            for t in range(np.size(barras_ang, axis = 0)):
                #Teta
                if de == barras_ang[t]:
        
                    H[i,col]= -V[de-1]*V[para-1]*(-G[de-1,para-1]*np.cos(teta[de-1]-teta[para-1]) - B[de-1,para-1]*np.sin(teta[de-1]-teta[para-1]))
        
                elif para == barras_ang[t]:
        
                    H[i,col]= V[de-1]*V[para-1]*(-G[de-1,para-1]*np.cos(teta[de-1]-teta[para-1]) - B[de-1,para-1]*np.sin(teta[de-1]-teta[para-1]))
                
                else:

                    H[i,col]=0

                col=col+1

            for v in range(np.size(barras_voltage, axis = 0)):
                 #V
                if de == barras_voltage[v]:
                    
                    tap=lado_medida(t,brt,G,B,shunt_bus,de-1,para-1)# é o tik
                    H[i,col]= -V[para-1]*(-G[de-1,para-1]*np.sin(teta[de-1]-teta[para-1]) + B[de-1,para-1]*np.cos(teta[de-1]-teta[para-1])) - 2*V[de-1]*(-tap*B[de-1,para-1] + shunt_bus[de-1,para-1])  
                
                elif para == barras_voltage[v]:
        
                    H[i,col] = -V[de-1]*(-G[de-1,para-1]*np.sin(teta[de-1]-teta[para-1]) + B[de-1,para-1]*np.cos(teta[de-1]-teta[para-1]))
                
                else:

                    H[i,col]=0

                col=col+1
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if tipo==5: #Iij
            col=0  #contabiliza coluna
            
            for t in range(np.size(barras_ang, axis = 0)):
                #Teta
                if de == barras_ang[t]:
                    
                    Iij = np.sqrt(((-G[de-1,para-1])**2 + (-B[de-1,para-1])**2)*(V[de-1]**2 + V[para-1]**2 -2*V[de-1]*V[para-1]*np.cos(teta[de-1]-teta[para-1])))
                    H[i,col]=  (((G[de-1,para-1]**2) + (B[de-1,para-1]**2))*V[de-1]*V[para-1]*np.sin(teta[de-1]-teta[para-1]))/Iij
        
                elif para == barras_ang[t]:
                    Iij = np.sqrt(((-G[de-1,para-1])**2 + (-B[de-1,para-1])**2)*(V[de-1]**2 + V[para-1]**2 -2*V[de-1]*V[para-1]*np.cos(teta[de-1]-teta[para-1])))
                    H[i,col]= -1*(((G[de-1,para-1]**2) + (B[de-1,para-1]**2))*V[de-1]*V[para-1]*np.sin(teta[de-1]-teta[para-1]))/Iij

                else:

                    H[i,col]= 0
                
                col=col+1

            for v in range(np.size(barras_voltage, axis = 0)):
                 #V
                if de == barras_voltage[v]:
                    Iij = np.sqrt(((-G[de-1,para-1])**2 + (-B[de-1,para-1])**2)*(V[de-1]**2 + V[para-1]**2 -2*V[de-1]*V[para-1]*np.cos(teta[de-1]-teta[para-1])))
                    H[i,col]= (((G[de-1,para-1]**2) + (B[de-1,]**2))*(V[de-1]-V[para-1]*np.cos(teta[de-1]-teta[para-1])))/Iij
                
                elif para == barras_voltage[v]:
                    Iij = np.sqrt(((-G[de-1,para-1])**2 + (-B[de-1,para-1])**2)*(V[de-1]**2 + V[para-1]**2 -2*V[de-1]*V[para-1]*np.cos(teta[de-1]-teta[para-1])))
                    H[i,col]= (((G[de-1,para-1]**2) + (B[de-1,]**2))*(V[para-1]-V[de-1]*np.cos(teta[de-1]-teta[para-1])))/Iij
                
                else:
        
                    H[i,col] = 0
                
                col=col+1

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
def  evalue_state(x_k,Ybus,std_dev,meas,num_states,num_meas,R,teta,V,num_bus,z,shunt_line,t,brt):
    x_new = np.zeros(num_states)
    H =jacobian(x_k,Ybus,meas,num_states,num_meas,teta,V,num_bus,shunt_line,t,brt)
    h= medidas_h(x_new,meas,Ybus,num_meas,num_states,teta,V,num_bus,shunt_line,t,brt)
    #Gain Matrix
    # G = ((H.transpose()).dot((np.linalg.inv(R)))).dot(H)
    G= np.transpose(H)*np.linalg.inv(R)*H
    # delta_x= (((np.linalg.inv(G)).dot((H.transpose()))).dot((np.linalg.inv(R)))).dot(np.subtract(z,h))
    r=z-np.transpose(h) #residual
    r=np.transpose(r)
    delta_x = np.linalg.inv(G)*(np.transpose(H)*np.linalg.inv(R)*r)
    delta_x=delta_x.tolist()
    J=0
    for i in range(num_meas):
        J+=(float(r[i])**2)/float(R[i,i])

    return  delta_x,h,H,J