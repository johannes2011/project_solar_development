#importing all the necessary stuff into python
import numpy as np
import matplotlib.pyplot as plt
import math
#remember, python starts counting at 0 NOT at 1

#%%
def Dataimporter(filename):
    Tin = np.loadtxt(filename, delimiter=",", usecols=[3], encoding='utf-8')
    Tout = np.loadtxt(filename, delimiter=",", usecols=[4], encoding='utf-8')
    Tambient = np.loadtxt(filename, delimiter=",", usecols=[7], encoding='utf-8')
    
    #outside temperature top-down
    TC1 = np.loadtxt(filename, delimiter=",", usecols=[8], encoding='utf-8')    
    TC2 = np.loadtxt(filename, delimiter=",", usecols=[9], encoding='utf-8')
    TC3 = np.loadtxt(filename, delimiter=",", usecols=[10], encoding='utf-8')
    TC4 = np.loadtxt(filename, delimiter=",", usecols=[11], encoding='utf-8')
    TC5 = np.loadtxt(filename, delimiter=",", usecols=[12], encoding='utf-8')
    
    #outside temperature bottom-up
    TC7 = np.loadtxt(filename, delimiter=",", usecols=[14], encoding='utf-8')
    TC8 = np.loadtxt(filename, delimiter=",", usecols=[15], encoding='utf-8')
    TC9 = np.loadtxt(filename, delimiter=",", usecols=[16], encoding='utf-8')
    TC10 = np.loadtxt(filename, delimiter=",", usecols=[17], encoding='utf-8')
    TC11 = np.loadtxt(filename, delimiter=",", usecols=[18], encoding='utf-8')
    
    Tleft=np.array([TC1,TC2,TC3,TC4,TC5])
    Tright=np.array([TC11,TC10,TC9,TC8,TC7])
    
    F3 = np.loadtxt(filename, delimiter=",", usecols=[21], encoding='utf-8')
    F2 = np.loadtxt(filename, delimiter=",", usecols=[22], encoding='utf-8')
    #here we F2 for the flow as it is much more precise. We also took care to remove
    #the constant offset of 0.04 that is always present in this measurement.
    return (Tin, Tout,Tambient,F2-0.04,Tleft,Tright)

#%%
def c_p_h2o (T_C):
    '''
c_p_h2o (T_C)
    Return the cp of the water at the given temperature 
    it uses the formula 
    to give a result in [kJ/kg/K]
    to then convert it into [kWh/kg/K]
    '''
    T=(T_C+273.15)/1000
    A_density = -203.606
    B_density = 1523.29
    C_density = -3196.413
    D_density = 2474.455
    E_density = 3.855326
    c_p=(A_density+B_density*T+C_density*T**2+D_density*T**3+E_density/T**2)/.01801528/1000
    return c_p/3600   # to have it in kwh/kg/K instead of kJ/kg/K

#%%
def rho_calc (T_C):
    '''
rho_calc (T_C) 
    calculating the density of water acording to Kell's formula
    in [kg/L]
    '''
    rho=(999.83952+16.945176*T_C-7.9870401*10**-3*T_C**2-46.170461*10**-6*T_C**3+105.56302*10**-9*T_C**4-280.54253*10**-12*T_C**5)/(1+16.897850*10**-3*T_C)
    return (rho/1000) #to have it in kg/L and not in kg/m^3

#%%
def heat_loss_coeff(filename):
    '''
    take as input the file name of the datasheet of the heat loss
    test and give back the heat loss coefficient in [kW/K]
    '''
    global h
    (Tin, Tout,Tambient,v,Tleft,Tright)=Dataimporter(filename)
    cp_in=c_p_h2o(Tin)
    rho_in=rho_calc(Tin)
    cp_out=c_p_h2o(Tout)
    rho_out=rho_calc(Tout)
    Tstorage=sum(Tleft+Tright)/10
    L=v*(cp_in*rho_in*Tin-cp_out*rho_out*Tout)/(Tstorage-Tambient)
    t=np.arange(0,len(L),1)
    
    plt.figure(dpi=360)
    ax1=plt.subplot(211)
    ax1=plt.plot(t/60,L*60*1000,label='heat loss coefficient')
    ax1=plt.ylim([0,10])
    ax1=plt.xlabel('Time [h]')
    ax1=plt.ylabel('h  [W/K]')
    ax1=plt.grid(True)
    ax1=plt.legend()
    ax2=plt.subplot(212)
    ax2=plt.plot(t/60,Tin,label='Tin')
    ax2=plt.plot(t/60,Tout,label='Tout')
    ax2=plt.plot(t/60,Tstorage,label='Tstorage')
    ax2=plt.xlabel('Time [h]')
    ax2=plt.ylabel('T °C')
    ax2=plt.grid(True)
    ax2=plt.legend()
    plt.savefig('Heatloss.jpg')
    
    h=np.average(L[3000:3500])
    #h is in the same unit as cp*kg/min so here kWh/K/min 
    print("\nThe heat loss coefficient is : %.3f [W/K] \n" %(h*60*1000))
    return h 

h=heat_loss_coeff('Pipe-in-cylinder_test5_heat_loss.txt')
#%%
def theoreticalEofT(T0,T1,Tsupercool):
    global mPCM,LPCM,cpPCMl,cpPCMs,msteel,Cpsteel,Vwater,Tphase
    Tmin=int(round(T0,0))
    Tmax=int(round(T1,0))
    Tphase=58
    Tx=[0]*6
    Ey=[0]*6
    msteel=166      #kg
    Cpsteel=0.12944444444/1000    #kWh/K/kg
    Vwater=75       #L
    mPCM=138        #kg
    cpPCMs=2.2/3600 #kWh/kg/K
    cpPCMl=3.1/3600 #kWh/kg/K
    
    #finding the Latent heat using test 4 and 5
    T40=21.5
    T41=88.6
    T50=21.3
    T51=81.6
    E4=22.840
    E5=20.980 
    E4s=(rho_calc(Tphase)*c_p_h2o(Tphase)*Tphase-rho_calc(T40)*c_p_h2o(T40)*T40)*Vwater+(msteel*Cpsteel+mPCM*cpPCMs)*(Tphase-T40)
    E5s=(rho_calc(Tphase)*c_p_h2o(Tphase)*Tphase-rho_calc(T50)*c_p_h2o(T50)*T50)*Vwater+(msteel*Cpsteel+mPCM*cpPCMs)*(Tphase-T50)
    E4l=-(rho_calc(Tphase)*c_p_h2o(Tphase)*Tphase-rho_calc(T41)*c_p_h2o(T41)*T41)*Vwater-(msteel*Cpsteel+mPCM*cpPCMl)*(Tphase-T41)
    E5l=-(rho_calc(Tphase)*c_p_h2o(Tphase)*Tphase-rho_calc(T51)*c_p_h2o(T51)*T51)*Vwater-(msteel*Cpsteel+mPCM*cpPCMl)*(Tphase-T51) 
    L4=1/mPCM*(E4-E4s-E4l)*3600
    L5=1/mPCM*(E5-E5s-E5l)*3600  
    
    LPCM=(L5)/3600          #kWh/kg
#    LPCM=(L4+L5)/2/3600     #kWh/kg
#    LPCM=265/3600           #kWh/kg
    
    Einit=rho_calc(T0)*c_p_h2o(T0)*Vwater*T0
    Ephase=rho_calc(Tphase)*c_p_h2o(Tphase)*Vwater*Tphase
    
    Tx[0]=T0
    Tx[1]=Tphase
    Tx[2]=Tphase
    Tx[3]=T1 
    Tx[4]=Tsupercool
    Tx[5]=Tphase
    
    Ey[0]=0
    Ey[1]=Ephase-Einit+(msteel*Cpsteel+mPCM*cpPCMs)*(Tphase-T0)
    Ey[2]=Ey[1]+LPCM*mPCM
    Ey[3]=(msteel*Cpsteel+mPCM*cpPCMs)*(Tphase-T0)+(rho_calc(T1)*c_p_h2o(T1)*Vwater*T1-Einit)+(msteel*Cpsteel+mPCM*cpPCMl)*(T1-Tphase)+LPCM*mPCM    
    Ey[4]=(msteel*Cpsteel+mPCM*cpPCMs)*(Tphase-T0)+(rho_calc(Tsupercool)*c_p_h2o(Tsupercool)*Vwater*Tsupercool-Einit)+(msteel*Cpsteel+mPCM*cpPCMl)*(Tsupercool-Tphase)+LPCM*mPCM
    Ey[5]=Ey[4]
    
    return(Tx,Ey)
#%%

def Ttheosto(Econtent,change,T0,T1,T2): 
    '''calculates the theoretical temperature of the storage using the theoretical curve'''
    Tx,Ey=theoreticalEofT(T0,T1,T2)
    n=len(Econtent)
    Tstoragetrue=[0]*n
    first=True
    second=True
    third=True
    fourth=False
    fifth=False
    for i in range(n):
        if Econtent[i]<Ey[1] and first:    
            Tstoragetrue[i]=Econtent[i]/(msteel*Cpsteel+mPCM*cpPCMs+Vwater*c_p_h2o(T0)*rho_calc(T0))+T0
        elif Econtent[i]<Ey[2] and second:
            first=False
            Tstoragetrue[i]=Tphase
        elif Econtent[i]>=Ey[2] and third:
            second=False
            fourth=True
            Tstoragetrue[i]=(Econtent[i]-Ey[2])/(msteel*Cpsteel+mPCM*cpPCMl+Vwater*c_p_h2o(Tphase)*rho_calc(Tphase))+Tphase
        elif Econtent[i]>=Ey[4] and fourth:
            third=False
            fifth=True
            Tstoragetrue[i]=(Econtent[i]-Ey[2])/(msteel*Cpsteel+mPCM*cpPCMl+Vwater*c_p_h2o(Tphase)*rho_calc(Tphase))+Tphase
        elif (Econtent[i]<Ey[4] and Econtent[i]>=Ey[1] and fifth) or (Ey[4]<=Ey[1] and fifth):
            fourth=False
            first=True
            Tstoragetrue[i]=Tphase
    return(Tstoragetrue)
#%%
def HXCR(Tin,Tout,Tstorage,namefig,t,V):
    #heat exchange capacity rate
    HX=np.zeros(len(Tstorage))
    rho=rho_calc((Tout+Tin)/2)
    cp=c_p_h2o ((Tout+Tin)/2)
    value=1-(Tin-Tout)/(Tin-Tstorage)
    for i in range(len(Tstorage)):
        if value[i]>0:
            HX[i]=-V[i]*rho[i]*cp[i]*np.log(value[i])
            
    #displaying the exchange rate in time
    plt.close()
    plt.figure(dpi=500)
    plt.plot(t/60,HX,label='exchange rate')
    plt.legend(loc=0)
    plt.grid(True)
    plt.xlabel('Time [h]')
    plt.ylabel('exchange rate')
    plt.title('exchange rate of '+namefig)
    plt.savefig(namefig+'ExchangeRate.jpg')
    plt.close()
    return()
#%%
def heat_content(filename,namefig,threshold):
    (Tin, Tout,Tambient,V,Tleft,Tright)=Dataimporter(filename)
    rho_in=rho_calc(Tin)
    cp_in=c_p_h2o (Tin)
    rho_out=rho_calc(Tout)
    cp_out=c_p_h2o (Tout)
    
    #h=heat_loss_coeff('Pipe-in-cylinder_test5_heat_loss.txt')
    Pflow=V*(cp_in*rho_in*Tin-cp_out*rho_out*Tout)*60 #kWh/min*60=kWh/h=kW
    Tstorage=sum(Tleft+Tright)/10
    Plost=h*(Tstorage-Tambient)*60      #kWh/min*60=kWh/h=kW
    Elost=sum(Plost/60)                 #kWh over the whole test
    t=np.arange(0,len(Plost),1)         #time
        
    Echarge=0                           #sum of all the flow in when the flow charging
    Edischarge=0                        #sum of all the flow in when flow discharging
    Esyphon=0                           #sum of all flows when the flow is lower than threshold
    indices=[]                          #indices of time when charging and discharging (three groups)
    Echargeinteg=np.zeros(len(Plost))   #vector Energy exchanged during charging phase 0 otherwise
    Esyphoninteg=np.zeros(len(Plost))   #vector Energy lost by syphon
    Esyphonlossinteg=np.zeros(len(Plost))   #vector Energy lost by syphon and losses
    Edischargeinteg=np.zeros(len(Plost))    #vector Energy exchanged during discharging phases 0 otherwise
    Elostinteg=np.zeros(len(Plost))     #vector Energy lost by losses
    Econtent=np.zeros(len(Plost))       #vector Energy content of the heat tank

    for i in range(1,len(Pflow)):
        Econtent[i]=Econtent[i-1]+(Pflow[i]-Plost[i])/60
        Elostinteg[i]=Elostinteg[i-1]-Plost[i]/60
        
        if V[i]>threshold and Pflow[i]>0:
            Echarge+=Pflow[i]/60
            Echargeinteg[i]=Echargeinteg[i-1]+Pflow[i]/60#-Plost[i]/60
            indices=indices+[i]
            
        elif V[i]>threshold and Pflow[i]<0:
            Edischarge-=Pflow[i]/60
            Edischargeinteg[i]=Edischargeinteg[i-1]+Pflow[i]/60
            indices=indices+[i]
            
        elif V[i]<threshold:
            Esyphon-=Pflow[i]/60
            Esyphoninteg[i]=Esyphoninteg[i-1]+Pflow[i]/60
            
        Esyphonlossinteg[i]=Esyphoninteg[i]+Elostinteg[i]
        
    balance=Echarge-Edischarge-Esyphon-Elost
    
    change=[]                           #change gives the indices of changing phase
    for i in range(1,len(indices)):
        if indices[i]-indices[i-1]>10:
            change=change+[indices[i-1]]
            change=change+[indices[i]]
    change=change+[indices[len(indices)-1]]
    
    Ttheoretical=Ttheosto(Econtent,change,Tstorage[0],max(Tstorage),Tstorage[change[2]+2])
    
    #displaying a figure with the Energy contents and exchange    
    plt.close()
    plt.figure(dpi=500)
#    plt.plot(t/60,Echargeinteg,label='Echarge')
#    plt.plot(t/60,Esyphoninteg,label='Esyphon')
#    plt.plot(t/60,Elostinteg,label='Elost')
#    plt.plot(t/60,Esyphonlossinteg,label='Esyphon and losses')
#    plt.plot(t/60,Edischargeinteg,label='Edischarge')
#    plt.plot([min(t)/60,max(t)/60],[max(Econtent),max(Econtent)],'--',label='Max energy content %.3f kWh'%(max(Econtent)))
    plt.plot(t/60,Econtent,label='Energy content')
    plt.legend(loc=0)
    plt.grid(True)
    plt.xlabel('Time [h]')
    plt.ylabel('Energy [kWh]')
    plt.title('Energy content for '+namefig)
    plt.savefig(namefig+'EnergyContent.jpg')
    
    #displaying Energy in function of T
#    namefig='theor'
    Tx,Ey=theoreticalEofT(Tstorage[0],max(Tstorage),Tstorage[change[2]+2])
    plt.close()
    plt.figure(dpi=500)
#    plt.plot([min(Tstorage),max(Tstorage)],[max(Econtent),max(Econtent)],'--',label='Max energy content %.3f kWh'%(max(Econtent)))
    plt.plot(Tstorage,Econtent,label='Measured Econtent f(measured T)') 
    plt.plot(Ttheoretical,Econtent,label='Measured Econtent f(theoretical T)') 
    plt.plot(Tx,Ey,'r--',label='Theorethical Energy content')
    plt.legend(loc=0)
    plt.grid(True)
    plt.xlabel('T [°C]')
    plt.ylabel('Energy content [kWh]')
    plt.title('Energy content function of T '+namefig)
    plt.savefig(namefig+'EnergyContentT.jpg')
    plt.close()
    
    
    
    #displaying HXCR
    HXCR(Tin,Tout,Ttheoretical,namefig,t,V)
    
    #displaying T in function of time
    plt.close()
    plt.figure(dpi=500)
    plt.plot(t/60,Tin,label='Tin')
    plt.plot(t/60,Tout,label='Tout')
    plt.plot(t/60,Tstorage,label='Tstorage measured')
    plt.plot(t/60,Ttheoretical,label='Tstorage theoretical')
#    plt.xlim([0,7.5])
    plt.legend(loc=0)
    plt.grid(True)
    plt.xlabel('Time [h]')
    plt.ylabel('T [°C]')
    plt.title('Temperature '+namefig)
    plt.savefig(namefig+'temperature.jpg')
    plt.close()
    
    Tstartcharge=Tstorage[0]
    Tendcharge=Tstorage[change[0]]
    TstartdischargeL=Tstorage[change[1]]
    TenddischargeL=Tstorage[change[2]]
    TstartdischargeS=Tstorage[change[3]]
    TenddischargeS=Tstorage[change[4]]  


#charge periods
    plt.close()
    plt.figure(dpi=500)
#    plt.plot([min(Tstorage),max(Tstorage)],[max(Econtent),max(Econtent)],'--',label='Max energy content %.3f kWh'%(max(Econtent)))
    plt.plot(Echarge,Econtent,label='Measured Econtent f(measured T)') 
    plt.plot(Ttheoretical,Econtent,label='Measured Econtent f(theoretical T)') 
    plt.plot(Tx,Ey,'r--',label='Theorethical Energy content')
    plt.legend(loc=0)
    plt.grid(True)
    plt.xlabel('T [°C]')
    plt.ylabel('Energy content [kWh]')
    plt.title('Energy content function of T '+namefig)
    plt.savefig(namefig+'EnergyContentT.jpg')
    plt.close()
    
    


    
    
    #printing some tables of data
    print('\n -------------------------------------\n'+namefig+'\n')
    print("Elost | Echarge | Edischarge | Esyphon | balance \n"
          "%.3f | %.3f  | %.3f     | %.3f   | %.3f   " %(Elost,Echarge,Edischarge,Esyphon,balance))
    print("             All values in kWh   \n")
    print('State      | Tstart   | Tend    |Econtent kWh \n'
          'charge     | %.1f     | %.1f    |%.3f ' %(Tstartcharge,Tendcharge,Econtent[change[0]])+'\n'
          'syphon     | %.1f     | %.1f    |%.3f ' %(Tendcharge,TstartdischargeL,Econtent[change[1]])+'\n'
          'discharge L| %.1f     | %.1f    |%.3f ' %(TstartdischargeL,TenddischargeL,Econtent[change[2]])+'\n'
          '+losses    | %.1f     | %.1f    |%.3f ' %(TenddischargeL,TstartdischargeS,Econtent[change[3]])+'\n'
          'discharge S| %.1f     | %.1f    |%.3f ' %(TstartdischargeS,TenddischargeS,Econtent[change[4]])+'\n')
    return Elost,Echarge,Edischarge,Esyphon,indices,balance,change

#%%      
          
#for the test 1 the volume exchanged is usually 5 L/min so a threshold of half is good for discriminating
Elost1,Echarge1,Edischarge1,Esyphon1,indices1,balance1,change1=heat_content('Pipe-in-cylinder_test1.txt','test1',1)

#for the test 2 the volume exchanged is usually 2 L/min so a threshold of half is good for discriminating
Elost2,Echarge2,Edischarge2,Esyphon2,indices2,balance2,change2=heat_content('Pipe-in-cylinder_test2.txt','test2',0.5)

#for the test 4 the volume exchanged is usually 10 L/min so a threshold of half is good for discriminating
Elost4,Echarge4,Edischarge4,Esyphon4,indices4,balance4,change4=heat_content('Pipe-in-cylinder_test4.txt','test4',3)

#for the test 5 the volume exchanged is usually 1 L/min so a threshold of half is good for discriminating
Elost5,Echarge5,Edischarge5,Esyphon5,indices5,balance5,change5=heat_content('Pipe-in-cylinder_test5_heat_loss.txt','test5',0.2)

