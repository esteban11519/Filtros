import math
import numpy as np
from numpy.lib.arraypad import _view_roi
from numpy.lib.polynomial import RankWarning
import os
import subprocess
import ltspice
import matplotlib.pyplot as plt


class Tools():
    def zita_Butterworth(self,Lc):
        '''
        Lc: Pérdida por inserción en omega_c [dB]
        '''
        return 10**(Lc/10)-1
    
    def n_Butterworth(self,omega,omega_c,zita,L):
        '''
        L: Pérdida por inserción en Omega

        Omega_c= Frecuencia en donde ocurre L_c

        zita: 10**(Lc/10)-1
        '''

        a=math.log10(10**(L/10)-1)
        b=math.log10(omega/omega_c)+math.log10(zita)
        return 0.5*a/b

    def gp_Butterworth_Rin_equal_out(self,n):
        '''
        n: Orden del filtro
        '''
        n=math.ceil(n)
        gp=[]
        gp.append(1)
        for p in range(1,n+1,1):
            a=(2*p-1)*math.pi
            b=2*n
            gp.append(2*math.sin(a/b))

        gp.append(1)
        return gp

    def Butterworth_pasa_bajas(self,coeficientes,omega_c,R_in,R_out):
        '''
        coeficientes: Coeficientes de Butterworth
        omega_c: Frecuencia de corte [rad/s]
        R_in: Resistencia de entrada [Ohmios]
        R_out: Resistencia de salida [Ohmios]
        '''
        
        n=len(coeficientes)
        frequency_scaling=np.array(coeficientes)
        frequency_scaling=frequency_scaling*1/omega_c
        frequency_scaling[0]=coeficientes[0]
        frequency_scaling[n-1]=coeficientes[n-1]

        impedance_scaling=frequency_scaling
        for p in range(0,n,1):

            if p==0 or p==n-1:
                impedance_scaling[0]=coeficientes[0]*R_in
                impedance_scaling[n-1]=coeficientes[n-1]=R_out
            elif p%2==1:
                impedance_scaling[p]=impedance_scaling[p]*R_in
            else :
                impedance_scaling[p]=impedance_scaling[p]*1/R_in
        return impedance_scaling


    def m_Chebyshev(self,L,Gr,omega,omega_c):
        '''
        L: Atenuaicón en Omega |[dB]|
        Gr: Rizado en la banda pasante en |[dB]|
        Omega_c: Frecuencia de corte 
        '''
        a=math.acosh(((10**(0.1*L)-1)/(10**(0.1*Gr)-1))**0.5)
        b=math.acosh(omega/omega_c)
        return a/b
    
    def gp_Chebyshev(self,m,Gr):
        '''
        m: Orden del filtro
        Gr: Rizado en la banda de paso [dB]
        '''
        m=math.ceil(m)
        gpList=[]
        gpList.append(1) # gp_0
        xi=math.log(1/math.tanh(Gr/17.37))
        ji=math.sinh(xi/(2*m))
        an1=math.sin(math.pi/(2*m))
        bn1=ji**2+math.sin(math.pi/m)**2
        gpn1=2*an1/ji
        gpList.append(gpn1) # gp_1


        for p in range(2,m+1,1):
            ap=math.sin((2*p-1)*math.pi/(2*m))
            bp=ji**2+math.sin(p*math.pi/m)**2
            a=4*an1*ap
            b=bn1*gpn1
            gp=a/b
            
            gpList.append(gp)

            an1=ap
            bn1=bp
            gpn1=gp

        if m%2==0:
            g_m_mas_1=1/math.tanh(xi/4)**2
        else:
            g_m_mas_1=1

        gpList.append(g_m_mas_1)
        return gpList
    

    def makeNetlistPasaBajas(self,coeficientes,f_min,f_max,puntosDecada,resultados):
        '''
        coeficientes: Coeficientes para el filtro pasa bajas incluyendo las resistencias
        f_min: Frecuenci mínima para simular
        f_max: frecuencia máxima para simular
        puntosDecada: Puntos que se quieren guardar por década
        resultados: Archivo para guardar el .raw de la simulación 
        '''

        elementos=len(coeficientes) # longitud total= n+1-0+1=n+2
        netlist=".title filtro pasa bajas\n"
        netlist+="Vin in 0 dc 0 ac 1 \n"

        # El recorrido se hade desde 0 hasta n ((n+2-1)-0)-1, el for llega hasta el límite superior -1
        #Se pudo hacer hasta n+1, pero este caso se trató cuando p==0
        for p in range(0,elementos-1,1):
            if p==0:
                netlist+='Rs in n0 '+str(coeficientes[p])+'\n'
                # si tengo n elementos reactivos y el primer nodo en serie es n0
                # ceil=techo 
                # ej: con 7 elementos reactivos el nodo final es de ceil(7/2)=4
                # ej: con 8 elementos reactivos el nodo final es de ceil(8/2)=4
                # A elementos se le resta 2 para que quede de tamaño n
                # Recordar que las listas de tamaño n, empiezan en la ubicación 0 y tenminan en n-1
                netlist+='RL n'+ str(math.ceil((elementos-2)/2)) + ' 0 '+ str(coeficientes[elementos-1]) +'\n'+'\n'

            elif p%2==1:
                #Esto hace referencia a los elementos inductivos que son impares
                #Floor=piso
                #floor(p/2) y  (p+1)/2 aseguran que los nodos sean consecutivos
                netlist+='L'+str(p)+' n'+str(math.floor(p/2))+' n'+str(int((p+1)/2))+' '+str(coeficientes[p])+'\n'
            else:
                #Esto hace referencia a los elementos Capacitivos que son pares
                #De manera empírica se tiene que el nodo superior en cada elemento es de la forma p/2
                netlist+='C'+str(p)+' n'+str(int(p/2))+' 0 '+str(coeficientes[p])+'\n'
        #El último nodo se reemplaza por 'out'
        netlist=netlist.replace('n'+str(math.ceil((elementos-2)/2)),'out')
        
        #Se definen las condiciones de simulación para ngspice
        netlist+=".control\n"\
        +"ac dec "+str(puntosDecada)+" "+str(f_min)+" "+str(f_max)+" \n"\
        +"*plot (vdb(out)-vdb(in))\n"\
        +"write "+os.getcwd()+'/'+resultados+" all \n"\
        +".endc\n\n"\
        +".end"
       
        print(netlist)
        return netlist
        
    def simular(self,archivo,t):
        process = subprocess.Popen(['ngspice',archivo])
        try:
            print('Running in process', process.pid)
            process.wait(timeout=t) #En segundos
        except subprocess.TimeoutExpired:
            print('Timed out - killing', process.pid)
            process.kill()
            print("Done")

    def graficar(self,archivo):
        l = ltspice.Ltspice(os.getcwd()+'/'+archivo) 
        # Make sure that the .raw file is located in the correct path
        l.parse()
        f= l.get_frequency()*1e-6
        V_in =abs(l.get_data('V(in)'))
        V_out =abs(l.get_data('V(out)'))
        H=20*np.log10(V_out)-20*np.log10(V_in)
        
        plt.semilogx(f,H ,f,np.ones(len(H))*-3,f,np.ones(len(H))*-30)
        plt.legend(['H','-3dB','-30dB'])
        plt.ylabel('Magnitud [dB]')
        plt.xlabel('frecuencia [MHz]')
        plt.grid(True)
        #plt.xticks(np.concatenate((np.arange(0.1,1,step=0.1),np.arange(1,11,step=1))))
        plt.show()
    
    def escritura(self,data,archivo):
        '''
        data: valores a guardar
        archivo: Documento donde se guarda
        '''
        try:
            archivo= open(os.getcwd()+'/'+archivo, "w" )
            archivo.write(data)
            archivo.close()
        except IOError:
            pass

if __name__ == '__main__':
    tool=Tools()
    
    R_in=50 # Resistencia de entrada [Ohmios]
    R_out=50 # Resistencia de salida [Ohmios]
    omega_c=10e6*2*math.pi # Frecuencia de corte [rad/s]
    omega=40e6*2*math.pi  # Frecuencia con atenuación L 
    L=30 # Atenuación en [dB] en omega

    # Parámetros de simulacion
    f_min=300e3 #Límite inferior para simular
    f_max=100e6 #Límite superior para simular
    puntosDecada=10000 #Número de puntos por década

    resultados="resultados.raw"
    circuito='filtro.cir'


    # Butterworth pasa bajas
    zita=tool.zita_Butterworth(3)
    n=tool.n_Butterworth(omega,omega_c,zita,L)
    print('n=',n)
    coeficientes=tool.gp_Butterworth_Rin_equal_out(n)
    print('Coeficientes: ',coeficientes)
    pasaBajas=tool.Butterworth_pasa_bajas(coeficientes,omega_c,R_in,R_out)
    print('pasaBajas',pasaBajas)
    data=tool.makeNetlistPasaBajas(pasaBajas,f_min,f_max,puntosDecada,resultados)
    tool.escritura(data,circuito)
    tool.simular(circuito,3)
    tool.graficar(resultados)

    #print('m',tool.m_Chebyshev(15,3,1.3,1))
    
    #for i in range(1,8,1):
    #    print(tool.gp_Chebyshev(i,3))

    
    
 
