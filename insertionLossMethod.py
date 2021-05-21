import math
import numpy as np
import os
import subprocess
import ltspice
import matplotlib.pyplot as plt



class Tools():
    def zita_Butterworth(self,Lc):
        '''
        Calcula zita para el filtro Butterworth
        Lc: Pérdida por inserción en omega_c [dB]
        '''
        return 10**(Lc/10)-1
    
    def n_Butterworth(self,omega,omega_c,zita,L):
        '''
        Calcula el orden para el filtro Butterworth
        L: Pérdida por inserción en Omega

        Omega_c= Frecuencia en donde ocurre L_c

        zita: 10**(Lc/10)-1
        '''

        a=math.log10(10**(L/10)-1)-math.log10(zita)
        b=math.log10(omega/omega_c)
        return 0.5*a/b

    def gp_Butterworth_Rin_equal_out(self,n):
        '''
        Calcula las inmitancias para el filtro Butterworth
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

   

    def m_Chebyshev(self,L,Gr,omega,omega_c):
        '''
        Calcula el orden para el filtro Chebyshev
        L: Atenuaicón en Omega |[dB]|
        Gr: Rizado en la banda pasante en |[dB]|
        Omega_c: Frecuencia de corte 
        '''
        a=math.acosh(((10**(0.1*L)-1)/(10**(0.1*Gr)-1))**0.5)
        b=math.acosh(omega/omega_c)
        return a/b
    
    def gp_Chebyshev(self,m,Gr):
        '''
        Calcula las inmitancias para el filtro Chebyshev
        m: Orden del filtro
        Gr: Rizado en la banda de paso [dB]
        '''
        m=math.ceil(m)
        
        #Para garantizar simetría en la resistencia de la fuente y carga
        if m%2==0:
            m+=1
        

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
        return list(gpList)
    
    def escalamiento_pasa_bajasoAltas(self,tipo,coeficientes,omega_c,R_in,R_out):
        '''
        Escalamiento en frecuencia y resistencia para filtros pasa bajas y pasa altas

        Nota: Solo aplica cuando las resistencias de entrada y salida son iguales
        coeficientes: Son los coeficientes de orden m Butterworth o Chebyshev
        tipo: 'pasa-altas' si desea un filtro pasa-altas
              'pasa-bajas' si desea un filtro pasa-bajas

        omega_c: Frecuencia de corte [rad/s]
        R_in: Resistencia de entrada [Ohmios]
        R_out: Resistencia de salida [Ohmios]
        '''
        inductorEsImpar=1
       

        n=len(coeficientes)
        frequency_scaling=np.array(coeficientes)
        if tipo=='pasa-altas':
            inductorEsImpar=0
            frequency_scaling=frequency_scaling**(-1)
        
        frequency_scaling=frequency_scaling*1/omega_c
        frequency_scaling[0]=coeficientes[0]
        frequency_scaling[n-1]=coeficientes[n-1]

        impedance_scaling=frequency_scaling
        for p in range(0,n,1):

            if p==0 or p==n-1:
                impedance_scaling[0]=coeficientes[0]*R_in
                impedance_scaling[n-1]=coeficientes[n-1]*R_out
            elif p%2==inductorEsImpar:
                impedance_scaling[p]=impedance_scaling[p]*R_in
            else :
                impedance_scaling[p]=impedance_scaling[p]*1/R_in

        return impedance_scaling


    def makeNetlistPasaBajasOAltas(self,tipo,coeficientes,f_min,f_max,puntosDecada,resultados):
        '''
        Realiza el netlist para los filtros pasa bajas o pasa altas

        coeficientes: Coeficientes para el filtro pasa bajas o Altas incluyendo las resistencias
        f_min: Frecuenci mínima para simular
        f_max: frecuencia máxima para simular
        puntosDecada: Puntos que se quieren guardar por década
        resultados: Archivo para guardar el .raw de la simulación 
        '''

        inductorEsImpar=1
        elementoImpar='L'
        elementoPar='C'
        titulo='Filtro pasa bajas\n'
        retardoDeNodo=0

        if tipo=='pasa-altas':
            inductorEsImpar=1
            titulo='Filtro pasa altas\n'
            retardoDeNodo=0
            elementoImpar='C'
            elementoPar='L'

        elementos=len(coeficientes) # longitud total= n+1-0+1=n+2
        netlist=".title "+titulo
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
                netlist+='RL n'+ str(math.ceil((elementos-2-retardoDeNodo)/2)) + ' 0 '+ str(coeficientes[elementos-1]) +'\n'+'\n'

            elif p%2==inductorEsImpar:
                #Esto hace referencia a los elementos inductivos que son impares
                #Floor=piso
                #floor(p/2) y  (p+1)/2 aseguran que los nodos sean consecutivos
                netlist+=elementoImpar+str(p)+' n'+str(math.floor((p-retardoDeNodo)/2))+' n'+str(int((p-retardoDeNodo+1)/2))+' '+str(coeficientes[p])+'\n'
            else:
                #Esto hace referencia a los elementos Capacitivos que son pares
                #De manera empírica se tiene que el nodo superior en cada elemento es de la forma p/2
                netlist+=elementoPar+str(p)+' n'+str(math.floor(p/2))+' 0 '+str(coeficientes[p])+'\n'
        #El último nodo se reemplaza por 'out'
        netlist=netlist.replace('n'+str(math.ceil((elementos-2-retardoDeNodo)/2)),'out')
        
        #Se definen las condiciones de simulación para ngspice
        netlist+=".control\n"\
        +"ac dec "+str(puntosDecada)+" "+str(f_min)+" "+str(f_max)+" \n"\
        +"*plot (vdb(out)-vdb(in))\n"\
        +"write "+os.getcwd()+'/'+resultados+" all \n"\
        +".endc\n\n"\
        +".end"
       
        print(netlist)
        return netlist

    def makeNetlistPasaBandas(self,coeficientes,omega_u,omega_l,R_inOR_out,f_min,f_max,puntosDecada,resultados):
        '''
        Realiza el netlist para el filtro pasa bandas

        coeficientes: Coeficientes para el filtro pasa bajas incluyendo las resistencias sin escalamiento
        omega_u: Frecuencia angular inferior [rad/s]
        omega_l: Frecuencia angular superior [rad/s]
        R_inOR_out: Resistencia de entrada o de la fuente [Ohmios]
        f_min: Frecuencia mínima para simular [Hz]
        f_max: frecuencia máxima para simular [Hz]
        puntosDecada: Puntos que se quieren guardar por década
        resultados: Archivo para guardar el .raw de la simulación 
        '''

        inductorEsImpar=1
        elementoImpar='L'
        elementoPar='C'
        titulo='Filtro pasa bandas\n'
        
        elementos=len(coeficientes) # longitud total= n+1-0+1=n+2
        netlist=".title "+titulo
        netlist+="Vin in 0 dc 0 ac 1 \n"

        # El recorrido se hade desde 0 hasta n ((n+2-1)-0)-1, el for llega hasta el límite superior -1
        #Se pudo hacer hasta n+1, pero este caso se trató cuando p==0
        for p in range(0,elementos-1,1):
            if p==0:
                netlist+='Rs in n0 '+str(coeficientes[p]*R_inOR_out)+'\n'
                # Recordar que las listas de tamaño n, empiezan en la ubicación 0 y tenminan en n-1
                netlist+='RL n'+ str(math.ceil((elementos-2)/2)) + ' 0 '+ str(coeficientes[elementos-1]*R_inOR_out) +'\n\n'

            elif p%2==inductorEsImpar:
                netlist+=elementoImpar+str(p)+' n'+str(math.floor((p)/2))+' nn'+str(math.floor((p)/2))+' '+str(self.inductorPasaBandas('serie',coeficientes[p],omega_u,omega_l,R_inOR_out))+'\n'
                netlist+=elementoPar+str(p)+' nn'+str(math.floor((p)/2))+' n'+str(int((p+1)/2))+' '+str(self.capacitorPasaBandas('serie',coeficientes[p],omega_u,omega_l,R_inOR_out))+'\n\n'
            else:
                netlist+=elementoImpar+str(p)+' n'+str(math.floor(p/2))+' 0 '+str(self.inductorPasaBandas('',coeficientes[p],omega_u,omega_l,R_inOR_out))+'\n'
                netlist+=elementoPar+str(p)+' n'+str(math.floor(p/2))+' 0 '+str(self.capacitorPasaBandas('',coeficientes[p],omega_u,omega_l,R_inOR_out))+'\n\n'
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
    
    def inductorPasaBandas(self,tipo,glOgc,omega_u,omega_l,R_inOR_out):
        '''
        Función interna de makeNetlistPasaBandas, la cual se encarga de hacer el escalamiento en frecuencia y resistencia de la inductancia 
        '''

        if tipo=='serie':
            return R_inOR_out*glOgc/(omega_u-omega_l)
        else:
            return R_inOR_out*(omega_u-omega_l)/((omega_u*omega_l)*glOgc)

    def capacitorPasaBandas(self,tipo,glOgc,omega_u,omega_l,R_inOR_out):
        '''
        Función interna de makeNetlistPasaBandas, la cual se encarga de hacer el escalamiento en frecuencia y resistencia de la capacitancia 
        '''
        if tipo=='serie':
            return (omega_u-omega_l)/((omega_u*omega_l)*glOgc*R_inOR_out)
        else:
            return glOgc/((omega_u-omega_l)*R_inOR_out)

    def makeNetlistStopBandas(self,coeficientes,omega_u,omega_l,R_inOR_out,f_min,f_max,puntosDecada,resultados):
        '''
        Realiza en netlist del filtro rechaza bandas

        coeficientes: Coeficientes para el filtro pasa bajas incluyendo las resistencias, sin escalamiento
        omega_u: Frecuencia angular inferior [rad/s]
        omega_l: Frecuencia angular superior [rad/s]
        R_inOR_out: Resistencia de entrada o de la fuente [Ohmios]
        f_min: Frecuencia mínima para simular [Hz]
        f_max: frecuencia máxima para simular [Hz]
        puntosDecada: Puntos que se quieren guardar por década
        resultados: Archivo para guardar el .raw de la simulación 
        '''

        titulo='Filtro rechaza bandas\n'
        
        elementos=len(coeficientes) # longitud total= n+1-0+1=n+2
        netlist="\n\n.title "+titulo
        netlist+="Vin in 0 dc 0 ac 1 \n"

        # El recorrido se hade desde 0 hasta n ((n+2-1)-0)-1, el for llega hasta el límite superior -1
        #Se pudo hacer hasta n+1, pero este caso se trató cuando p==0
        for p in range(0,elementos-1,1):
            if p==0:
                netlist+='Rs in n0 '+str(coeficientes[p]*R_inOR_out)+'\n'
                # Recordar que las listas de tamaño n, empiezan en la ubicación 0 y tenminan en n-1
                netlist+='RL n'+ str(math.ceil((elementos-2)/2)) + ' 0 '+ str(coeficientes[elementos-1]*R_inOR_out) +'\n\n'

            elif p%2==1:
                netlist+='L'+str(p)+' n'+str(math.floor((p)/2))+' n'+str(int((p+1)/2))+' '+str(self.inductorStopBandas('serie',coeficientes[p],omega_u,omega_l,R_inOR_out))+'\n'
                netlist+='C'+str(p)+' n'+str(math.floor((p)/2))+' n'+str(int((p+1)/2))+' '+str(self.capacitorStopBandas('serie',coeficientes[p],omega_u,omega_l,R_inOR_out))+'\n\n'
            else:
                netlist+='L'+str(p)+' n'+str(math.floor(p/2))+' s'+str(math.floor(p/2))+' '+str(self.inductorStopBandas('',coeficientes[p],omega_u,omega_l,R_inOR_out))+'\n'
                netlist+='C'+str(p)+' s'+str(math.floor(p/2))+' '+' 0 '+str(self.capacitorStopBandas('',coeficientes[p],omega_u,omega_l,R_inOR_out))+'\n\n'
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
    
    def inductorStopBandas(self,tipo,glOgc,omega_u,omega_l,R_inOR_out):
        '''
        Función interna de makeNetlistStopBandas
        '''

        if tipo=='serie':
            return R_inOR_out*glOgc*(omega_u-omega_l)/(omega_u*omega_l)
        else:
            return R_inOR_out/((omega_u-omega_l)*glOgc)

    def capacitorStopBandas(self,tipo,glOgc,omega_u,omega_l,R_inOR_out):
        '''
        Función interna de makeNetlistStopBandas
        '''
        if tipo=='serie':
            return 1/((omega_u-omega_l)*glOgc*R_inOR_out)
        else:
            return (omega_u-omega_l)*glOgc/((omega_u*omega_l)*R_inOR_out)

        

    def simular(self,archivo,t):
        '''
        Esta funcíon permite simular el netlist en ngspice
        Archivo: Netlist en .cir
        t: Tiempo de simulación en segundos
        '''
        process = subprocess.Popen(['ngspice',archivo])
        try:
            print('Running in process', process.pid)
            process.wait(timeout=t) #En segundos
        except subprocess.TimeoutExpired:
            print('Timed out - killing', process.pid)
            process.kill()
            print("Done")

    def graficar(self,archivo,tipoFiltro,orden):
        '''
        Esta función permite graficar los resultados

        - Archivo: Documento .raw que se genera por el simulador
        - tipoFiltro: Tipo de filtro
        - orden: Orden del filtro
        '''
        l = ltspice.Ltspice(os.getcwd()+'/'+archivo) 
        # Make sure that the .raw file is located in the correct path
        l.parse()
        f= l.get_frequency()*1e-6
        V_in =abs(l.get_data('V(in)'))
        V_out =abs(l.get_data('V(out)'))
        H=20*np.log10(V_out)-20*np.log10(V_in)
        
        plt.figure()
        plt.semilogx(f,H ,f,np.ones(len(H))*-6,f,np.ones(len(H))*-30,f,np.ones(len(H))*-73.87)
        plt.ylabel('Magnitud [dB]')
        plt.xlabel('frecuencia [MHz]')
        
        plt.axvline(x=200, ymin=-3, ymax=3,c='m')
        plt.legend(['H','-6dB','-30dB','-73.87dB','200 MHz'])
        plt.grid(True)
        plt.title(tipoFiltro+' con orden '+str(orden))
        plt.xticks(np.concatenate((np.arange(10,100,step=10),np.arange(100,400,step=100))))
        plt.savefig(os.getcwd()+'/'+tipoFiltro.replace(' ','_').lower()+'_Orden_'+str(orden)+'.png')
        plt.show()
    
    def escritura(self,data,archivo):
        '''
        Permite escribir un archivo

        data: valores a guardar
        archivo: Documento donde se guarda
        '''
        try:
            archivo= open(os.getcwd()+'/'+archivo, "w" )
            archivo.write(data)
            archivo.close()
        except IOError:
            pass
    def BWAndOmega_02Omega_lAndOmega_u(self,BW,omega_o):
        '''
        Pasa de un ancho de banda y una frecuencia central a una frecuencia de corte inferior y otra superior
        
        - BW: Ancho de banda [rad/s]
        - omega_o: Frecuencia central [rad/s]
        
        '''
        Omega_lAndOmega_u=[];
        Omega_lAndOmega_u.append((-BW+(BW**2+4*omega_o**2)**0.5)/2)
        Omega_lAndOmega_u.append(Omega_lAndOmega_u[0]+BW)
        return Omega_lAndOmega_u

    def filtro(self,tipo,polinomio,omega_c,omega,omega_l,omega_u,L_c,L,Gr,RinORout,ordenPasaOrechazaBandas,f_min,f_max,puntosDecada,circuito,resultados,tiempoSimulacion):
        
        '''
        Permite diseñar, simular y graficar un filtro pasa bajas, pasa altas, pasa bandas y rechaza bandas.

        tipo:
            - Colocar `LP' si un pasa bajas
            - Colocar `HP' si un pasa altas
            - Colocar `BP' Si un pasa bandas
            - Colocar `' si un rechazabandas
        
        Polinomio: 
            - Colocar `B' si es un Butterworth
            - Colocar `' si es un Chebyshev

        omega_c:
            Solo aplica para pasa bajas o altas
            -Frecuencia angular [rad/s] de corte

        omega:
            Solo aplica para pasa bajas o altas
            - Frecuencia en [rad/s] donde se atenua por lo menos L[dB]

        omega_l:
            Solo aplica para pasa bandas o rechazabandas
            - Frecuencia de corte inferior [rad/s]

        omega_u:
            Solo aplica para pasa bandas o rechazabandas
            - Frecuencia de corte superior [rad/s]
        
        L_c:
            Solo aplica para pasa bajas o altas
            - Atenuación en [dB] que sucede en omega_c

        L:
            Solo aplica para pasa bajas o altas
            - Atenuación en [dB] en la frecuencia omega

        Gr: 
            Solo aplica para el tipo Chebyshev
            - Rizado en la banda pasante pico-pico en [dB]
        
        RinORout: 
            Resistencia de entrada o de la salida [Ohmios], suponiendo que son iguales
        
        ordenPasaOrechazaBandas:
            Solo aplica para pasa bandas o rechazabandas
            - Orden del filtro pasa o rechaza bandas

        Parámetros de simulación

        f_min:
            - Frecuencia mínima de simulación [Hz]
        
        f_max:
            - Frecuencia máxima de simulación [Hz]
        
        puntosDecada:
            - Número de puntos por década

        circuito:
            - Archivo para escribir el netlist, debe estár en .cir
        
        resultados:
            - Archivo para guardar los resultados, su extensión debe estar en .raw

        tiempoSimulacion:
            - Tiempo de simulación [segundos]
        '''

        resumenSimulacion="\n\nEste es la simulaci\243n que usted ingres\243\n\n"
        if tipo=='LP':
            resumenSimulacion+='Tipo: Filtro pasa bajas\n'
        elif tipo=='HP':
            resumenSimulacion+='Tipo: Filtro pasa altas\n'
        elif tipo=='BP':
            resumenSimulacion+='Tipo: Filtro pasa bandas\n'
        else:
            resumenSimulacion+='Tipo: Filtro rechaza bandas\n'

        if polinomio=='B':
            resumenSimulacion+='Polinomio: Butterworth'
        else:
            resumenSimulacion+='Polinomio: Chebyshev'

        if tipo=='LP' or tipo=='HP':
            resumenSimulacion+='\nomega_c: '+str(omega_c)+' [rad/s]'\
            +'\nomega: '+str(omega)+' [rad/s]'\
            +'\nL_c: '+str(L_c)+' [dB]'\
            +'\nL: '+str(L)+' [dB]'
        
        else:
            resumenSimulacion+='\nomega_l: '+str(omega_l)+' [rad/s]'\
            +'\nomega_u: '+str(omega_u)+' [rad/s]'\
            +'\norden del filtro pasa o rechaza bandas: '+str(ordenPasaOrechazaBandas)
        
        if polinomio=='': 
            resumenSimulacion+='\nRizado en la banda pasante: '+str(Gr)+' [dB]'

        resumenSimulacion+='\nResistencia de la fuente o de carga: '+str(RinORout)+' [Ohm]'\
        +'\nf_min: '+str(f_min)+' [Hz]'\
        +'\nf_max: '+str(f_max)+' [Hz]'\
        +'\nPuntos por decada: '+str(puntosDecada)\
        +'\ncircuito: '+str(circuito)\
        +'\nresultados: '+str(resultados)\
        +'\nTiempo de simulaci\243n: '+str(tiempoSimulacion)+'\n\n'

        print(resumenSimulacion)

        
        
        print('Prueba tipo',tipo)
        if tipo=='pasa-bajas' or tipo=='LP':
            tipoFiltro='Pasa bajas'
            print('pasa-bajas')
            if polinomio=='Butterworth' or polinomio=='B':
                zita=self.zita_Butterworth(L_c)
                orden=n=self.n_Butterworth(omega,omega_c,zita,L)
                coeficientes=self.gp_Butterworth_Rin_equal_out(n)
            else:
                orden=m=self.m_Chebyshev(L,Gr,omega,omega_c)
                coeficientes=self.gp_Chebyshev(m,Gr)

            escalamiento=self.escalamiento_pasa_bajasoAltas('pasa-bajas',coeficientes,omega_c,RinORout,RinORout)
            netlist=self.makeNetlistPasaBajasOAltas('pasa-bajas',escalamiento,f_min,f_max,puntosDecada,resultados)
            
            
        elif tipo=='pasa-altas' or tipo== 'HP':
            tipoFiltro='Pasa altas'
            print('pasa-altas')
            if polinomio=='Butterworth' or polinomio=='B':
                zita=self.zita_Butterworth(L_c)
                orden=n=self.n_Butterworth(1/omega,1/omega_c,zita,L)
                coeficientes=self.gp_Butterworth_Rin_equal_out(n)
            else:
                orden=m=self.m_Chebyshev(L,Gr,1/omega,1/omega_c)
                coeficientes=self.gp_Chebyshev(m,Gr)
                
            escalamiento=self.escalamiento_pasa_bajasoAltas('pasa-altas',coeficientes,omega_c,RinORout,RinORout)
            netlist=self.makeNetlistPasaBajasOAltas('pasa-altas',escalamiento,f_min,f_max,puntosDecada,resultados)
            
        elif tipo=='pasa-bandas' or tipo=='BP':
            tipoFiltro='Pasa bandas'
            print('pasa-bandas')
            if polinomio=='Butterworth' or polinomio=='B':
                coeficientes=self.gp_Butterworth_Rin_equal_out(ordenPasaOrechazaBandas)
            else:
                coeficientes=self.gp_Chebyshev(ordenPasaOrechazaBandas,Gr)
            
            orden=ordenPasaOrechazaBandas
            netlist=self.makeNetlistPasaBandas(coeficientes,omega_u,omega_l,RinORout,f_min,f_max,puntosDecada,resultados)

        else:
            tipoFiltro='Rechaza bandas'
            print('Stop-bandas')
            if polinomio=='Butterworth'or polinomio=='B':
                coeficientes=self.gp_Butterworth_Rin_equal_out(ordenPasaOrechazaBandas)
            else:
                coeficientes=self.gp_Chebyshev(ordenPasaOrechazaBandas,Gr)
            
            orden=ordenPasaOrechazaBandas
            netlist=self.makeNetlistStopBandas(coeficientes,omega_u,omega_l,RinORout,f_min,f_max,puntosDecada,resultados)


        self.escritura(netlist,circuito)
        self.simular(circuito,tiempoSimulacion)
        self.graficar(resultados,tipoFiltro,orden)
       
if __name__ == '__main__':
    tool=Tools()
    tipo='BP'
    polinomio=''
    R_in=75 # Resistencia de entrada [Ohmios]
    R_out=75 # Resistencia de salida [Ohmios]
    BW=10e6*2*math.pi # Ancho de banda [rad/s]
    omega_o=90e6*2*math.pi # Ancho de banda [rad/s]

    conversion=tool.BWAndOmega_02Omega_lAndOmega_u(BW,omega_o)
    omega_u=conversion[1] # Frecuencia de corte superior [rad/s]
    omega_l=conversion[0] # Frecuencia de corte inferior [rad/s]
    
    Gr=0.1 # Rizado en el banda pasante [dB]
    
    # Parámetros de simulacion
    f_min=27e6 #Límite inferior para simular
    f_max=300e6 #Límite superior para simular
    puntosDecada=10000 #Número de puntos por década

    resultados="resultados.raw"
    circuito='filtroPasaBandas.cir'
    tiempoSimulacion=3
    
    for ordenPasaOrechazaBandas in range(1,7,2):
        tool.filtro(tipo,polinomio,'omega_c','omega',omega_l,omega_u,'L_c','L',Gr,R_in,ordenPasaOrechazaBandas,f_min,f_max,puntosDecada,circuito,resultados,tiempoSimulacion)
     
