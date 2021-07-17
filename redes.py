import numpy as np
import math
import os
import subprocess
import ltspice
import matplotlib.pyplot as plt

class Tool():
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

    def RvirtualTipoT(self,Rs,RL,Q):
        '''
        Rs: Resistencia de entrada [Ohmios]
        RL: Resistencia de salida [Ohmios]
        Q: Factor de calidad inicial
        '''
        return (Q**2+1)*min(Rs,RL)

    
    def RvirtualTipoPi(self,Rs,RL,Q):
        '''
        Rs: Resistencia de entrada [Ohmios]
        RL: Resistencia de salida [Ohmios]
        Q: Factor de calidad inicial
        '''
        
        return  max(Rs,RL)/(Q**2+1)

    def RvirtualTipoTAndPi(self,tipoRed,Rs,RL,Q):
        '''
        tipoRed: 
            - 't' : Si desea una red T
            - '': Si desea una red pi
                   
        Rs: Resistencia de entrada [Ohmios]
        RL: Resistencia de salida [Ohmios]
        Q: Factor de calidad inicial
        '''
        
        if tipoRed.lower()=='t':
            return self.RvirtualTipoT(Rs,RL,Q)
        else:
            return self.RvirtualTipoPi(Rs,RL,Q)
    
    def coeficientesXsxXpxTipoT(self,Rvirtual,Rs,RL,Q):
        '''
        Rvirtual: Resistencia Virtual [Ohmios]
        Rs: Resistencia de la fuente [Ohmios]
        RL: Resistencia de la Carga [Ohmios]
        Q: Factor de calidad deseado
        '''
        coeficientes=[]
        if RL<=Rs:
            Xs2=RL*Q
            Xp2=Rvirtual/Q
            Qnuevo=(Rvirtual/Rs-1)**0.5
            Xs1=Rs*Qnuevo
            Xp1=Rvirtual/Qnuevo

        else:
            Xs1=Rs*Q
            Xp1=Rvirtual/Q
            Qnuevo=(Rvirtual/RL-1)**0.5
            Xs2=RL*Qnuevo
            Xp2=Rvirtual/Qnuevo

        coeficientes.append(Rs)
        coeficientes.append(Xs1)
        coeficientes.append(Xp1)
        coeficientes.append(Xp2)
        coeficientes.append(Xs2)
        coeficientes.append(RL)
        print("\n\treactancias red tipo T\n")
        print(f"\
        Rs={coeficientes[0]}\n\
        Xs1={coeficientes[1]}\n\
        Xp1={coeficientes[2]}\n\
        Xp2={coeficientes[3]}\n\
        Xs2={coeficientes[4]}\n\
        RL={coeficientes[5]}\n")
    
        return coeficientes

    def coeficientesXsxXpxTipoPi(self,Rvirtual,Rs,RL,Q):
        '''
        Rvirtual: Resistencia Virtual [Ohmios]
        Rs: Resistencia de la fuente [Ohmios]
        RL: Resistencia de la Carga [Ohmios]
        Q: Factor de calidad deseado
        '''
        coeficientes=[]
        if RL>=Rs:
            Xs2=Rvirtual*Q
            Xp2=RL/Q
            Qnuevo=(Rs/Rvirtual-1)**0.5
            Xs1=Rvirtual*Qnuevo
            Xp1=Rs/Qnuevo

        else:
            Xs1=Rvirtual*Q
            Xp1=Rs/Q
            Qnuevo=(RL/Rvirtual-1)**0.5
            Xs2=Rvirtual*Qnuevo
            Xp2=RL/Qnuevo

        coeficientes.append(Rs)
        coeficientes.append(Xs1)
        coeficientes.append(Xp1)
        coeficientes.append(Xp2)
        coeficientes.append(Xs2)
        coeficientes.append(RL)
        print("\n\treactancias red tipo PI\n")
        print(f"\
        Rs={coeficientes[0]}\n\
        Xp1={coeficientes[2]}\n\
        Xs1={coeficientes[1]}\n\
        Xs2={coeficientes[4]}\n\
        Xp2={coeficientes[3]}\n\
        RL={coeficientes[5]}\n")
    

        return coeficientes

    def coeficientesXsxXpxTipoTAndPi(self,tipoRed,Rvirtual,Rs,RL,Q):
        '''
        tipoRed: 
            - 't' : Si desea una red T
            - '': Si desea una red pi
                   
        Rvirtual: Resistencia Virtual [Ohmios]
        Rs: Resistencia de la fuente [Ohmios]
        RL: Resistencia de la Carga [Ohmios]
        Q: Factor de calidad deseado
        '''

        if tipoRed.lower()=='t':
            return self.coeficientesXsxXpxTipoT(Rvirtual,Rs,RL,Q)
        else:
            return self.coeficientesXsxXpxTipoPi(Rvirtual,Rs,RL,Q)


    def reactanciaParalela(self,X1,X2):

        return X1*X2/(X1+X2)
    
    def netlistTipoTAndPi(self,omega_0,coeficientes,Xs1,Xs2,tipoRed,f_min,f_max,puntosDecada,resultados):
        '''
        omega_0: Frecuencia de resonancia [rad/s] (Tipo=float)

        coeficientes: Reactancias sin normalizar (Tipo=Array)

        Xs1: -1 si Xs1 es capacitor y 1 si Xs1 es inductor (Tipo=Entero)

        Xs2: -1 si Xs2 es capacitor y 1 si Xs2 es inductor (Tipo=Entero)

        tipoRed: '' Si desea una red pi y 't' Si desea una red T (Tipo=char)

        f_min: Frecuencia mínima de simulación [Hz] (Tipo=float)

        f_max: Frecuencia máxima de simulación [Hz] (Tipo=float)
        
        puntosDecada: Puntos por década (Tipo=entero positivo)

        resultados: Archivo .raw donde se deseen guardar los resultados de simulación (Tipo=String)
        '''
        
        if tipoRed.lower()=='t':

            netlist="\n\n.title Circuito de acople T \n"
            netlist+="Vin 0 in dc 0 ac 1 \n"
            netlist+="Rs in n1 "+str(coeficientes[0])+"\n"
            netlist+="RL out 0 "+str(coeficientes[5])+"\n\n "

            if Xs1>0:
                # Las reactancias paralelas y series deben ser opuestas
                Elemento_Xs1='Ls1'
                Valor_Xs1=coeficientes[1]/omega_0
            else:
                Elemento_Xs1='Cs1'
                Valor_Xs1=1/(coeficientes[1]*omega_0)
                
            
            if Xs2>0:
                # Las reactancias paralelas y series deben ser opuestas
                Elemento_Xs2='Ls2'
                Valor_Xs2=coeficientes[4]/omega_0

            else:
                Elemento_Xs2='Cs2'
                Valor_Xs2=1/(coeficientes[4]*omega_0)
            

            # reactancia equivalente
            
            Xp_eq=self.reactanciaParalela((-Xs1)*coeficientes[2],(-Xs2)*coeficientes[3])

            if Xp_eq>0:
                Elemento_Xp_eq='Lp_eq'
                Valor_Xp_eq=Xp_eq/omega_0
            else:
                Elemento_Xp_eq='Cp_eq'
                Valor_Xp_eq=1/((-Xp_eq)*omega_0)

            
            netlist+=f"{Elemento_Xs1} n1 n2 {Valor_Xs1}\n "\
            +f"{Elemento_Xp_eq} n2 0 {Valor_Xp_eq}\n "\
            +f"{Elemento_Xs2} n2 out {Valor_Xs2}\n\n "

        else:
            netlist="\n\n.title Circuito de acople Pi \n"
            netlist+="Vin 0 in dc 0 ac 1 \n"
            netlist+="Rs in n1 "+str(coeficientes[0])+"\n"
            netlist+="RL out 0 "+str(coeficientes[5])+"\n\n "

            if Xs1>0:
                # Las reactancias paralelas y series deben ser opuestas
                Elemento_Xp1='Cp1'
                Valor_Xp1=1/(coeficientes[2]*omega_0)
            else:
                Elemento_Xp1='Lp1'
                Valor_Xp1=coeficientes[2]/omega_0
                
            
            if Xs2>0:
                # Las reactancias paralelas y series deben ser opuestas
                Elemento_Xp2='Cp2'
                Valor_Xp2=1/(coeficientes[3]*omega_0)
            else:
                Elemento_Xp2='Lp2'
                Valor_Xp2=coeficientes[3]/omega_0
       
            # reactancia equivalente

            X_s_eq=Xs1*coeficientes[1]+Xs2*coeficientes[4]

            if X_s_eq>0:
                Elemento_Xs_eq='Ls_eq'
                Valor_Xs_eq=X_s_eq/omega_0
            else:
                Elemento_Xs_eq='Cs_eq'
                Valor_Xs_eq=1/((-X_s_eq)*omega_0)

            
            netlist+=f"{Elemento_Xp1} n1 0 {Valor_Xp1}\n "\
            +f"{Elemento_Xs_eq} n1 out {Valor_Xs_eq}\n "\
            +f"{Elemento_Xp2} out 0 {Valor_Xp2}\n\n "

        netlist+=".control\n"\
        +"ac dec "+str(puntosDecada)+" "+str(f_min)+" "+str(f_max)+" \n"\
        +"*plot (vdb(out)-vdb(in))\n"\
        +"write "+os.getcwd()+'/'+resultados+" all \n"\
        +".endc\n\n"\
        +".end"

        print(netlist)
        return netlist


    def posicionMagnitudMaxima(self,y):
        nMax=0;
        fMax=y[0];

        for i in range(len(y)-1):
            if y[i+1]>fMax:
                fMax=y[i+1]
                nMax=i+1

        return nMax

    def posicionPrimerCoicidenciaSentidoPositivo(self,y,valorDeseado):
        
        for i in range(len(y)):
            if y[i] >= valorDeseado:
                return i
        
        return -1    

    def posicionPrimerCoicidenciaSentidoNegativo(self,y,valorDeseado):
        
        for i in range(len(y)-1,-1,-1):
            if y[i] >= valorDeseado:
                return i
        
        return -1    


    def graficar(self,archivo,Rs,RL):
        # Obtenención de datos
        l = ltspice.Ltspice(os.getcwd()+'/'+archivo) 
        l.parse()
        
        unidadesFrecuencia='MHz'
        unidadesFrecuenciaFactor=1e-6

        f= l.get_frequency()*unidadesFrecuenciaFactor
        v_in =l.get_data('V(in)')
        i_vin=l.get_data('i(vin)')
        v_out =l.get_data('V(out)')

         # Operaciones
        H=20*np.log10(abs(v_out/v_in))
        H_angle=(np.angle(v_out/v_in))*180/math.pi
    
        dB_corte=20*np.log10(1/(2**0.5))

        nMax=self.posicionMagnitudMaxima(H)
        f_resonancia=f[nMax]

        nOmegaLow=self.posicionPrimerCoicidenciaSentidoPositivo(H,H[nMax]+dB_corte)
        f_low=f[nOmegaLow]

        nOmegaHigh=self.posicionPrimerCoicidenciaSentidoNegativo(H,H[nMax]+dB_corte)
        f_high=f[nOmegaHigh]

        Q=f_resonancia/(f_high-f_low)

        print('f_o: ',f_resonancia,unidadesFrecuencia)
        print(f'Q={Q} {unidadesFrecuencia}/{unidadesFrecuencia}')  

        # Configuración de texto

        conf_left = {'fontsize':12,'horizontalalignment':'left','verticalalignment':'top'}
        conf_right= {'fontsize':12,'horizontalalignment':'right','verticalalignment':'top'}
        conf_center= {'fontsize':12,'horizontalalignment':'center','verticalalignment':'baseline'}


        ## Factor de calidad

        fig, axs = plt.subplots(2, 1)
        axs[0].semilogx(f,H,f,np.ones(len(H))*(H[nMax]+dB_corte))
        axs[1].semilogx(f,H_angle)
        
        axs[0].set_xticks((np.arange(70,110,step=10)))
        axs[1].set_xticks((np.arange(70,110,step=10)))
        
        axs[0].text(f[nMax],H[nMax]\
            ,f'({f[nMax].round(3)} [{unidadesFrecuencia}],{H[nMax].round(3)}[dB])',conf_center)
        
        #H_c_low=H[nOmegaLow]+(H[0]-H[nOmegaLow])/8
        axs[0].text(f[nOmegaLow],H[nOmegaLow]\
            ,f'({f[nOmegaLow].round(3)} [{unidadesFrecuencia}],{H[nOmegaLow].round(3)}[dB])',conf_right)
        #H_c_high=H[nOmegaHigh]+(H[0]-H[nOmegaHigh])/8
        axs[0].text(f[nOmegaHigh],H[nOmegaHigh]\
            ,f'({f[nOmegaHigh].round(3)} [{unidadesFrecuencia}],{H[nOmegaHigh].round(3)}[dB])',conf_left)
        H_m=H[nMax]+(H[0]-H[nMax])/2
        axs[0].text(f[nMax],H_m\
            ,f'Q={Q.round(3)}',conf_center)
        
        
        axs[1].text(f[nMax],H_angle[nMax]*0.97\
            ,f'({f[nMax].round(3)} [{unidadesFrecuencia}],{H_angle[nMax].round(3)}[$^0 $])',conf_left)
        
        axs[0].set_xlabel(f'frecuencia [{unidadesFrecuencia}]')
        axs[1].set_xlabel(f'frecuencia [{unidadesFrecuencia}]')
        
        axs[0].set_ylabel('Magnitud [dB]')
        axs[1].set_ylabel('Arg(H) [$^0 $]')

        axs[0].legend(['H','-3 [dB]'],loc="lower right")

        axs[0].grid(True)
        axs[1].grid(True)
        plt.show()

        
        
        ## Impedancia de entrada del amplificador

        Z_in_amplificador=abs(v_in/i_vin)-Rs
        Z_in_amplificador_fase=np.angle(v_in/i_vin)*180/math.pi


        fig, axs = plt.subplots(2, 1)
        axs[0].semilogx(f, Z_in_amplificador)
        axs[1].semilogx(f, Z_in_amplificador_fase)
        
        axs[0].set_xticks((np.arange(70,110,step=10)))
        axs[1].set_xticks((np.arange(70,110,step=10)))
        
        axs[0].text(f[nMax],Z_in_amplificador[nMax]*0.97\
            ,f'({f[nMax].round(3)} [{unidadesFrecuencia}],{Z_in_amplificador[nMax].round(3)}[$\Omega $])',conf_left)
        axs[1].text(f[nMax],Z_in_amplificador_fase[nMax]*0.97\
            ,f'({f[nMax].round(3)} [{unidadesFrecuencia}],{Z_in_amplificador_fase[nMax].round(3)}[$^0 $])',conf_left)
        
        axs[0].set_xlabel(f'frecuencia [{unidadesFrecuencia}]')
        axs[1].set_xlabel(f'frecuencia [{unidadesFrecuencia}]')
        
        axs[0].set_ylabel('|Z|')
        axs[1].set_ylabel('Arg(Z) [$^0 $]')

        axs[0].grid(True)
        axs[1].grid(True)
        plt.show()

        ## Potencia
        S_fuente=v_in*np.conj(i_vin) #[VA] Potencia aparente
        S_carga=v_out*np.conj(v_out/RL) #[VA] Potencia aparente

        P_fuente=np.real(S_fuente) #[W] Potencia real
        P_carga=np.real(S_carga) #[[W]] Potencia real
        
        Q_fuente=np.imag(S_fuente) #[W] Potencia real
        Q_carga=np.imag(S_carga) #[[W]] Potencia real

        H_P=P_carga/P_fuente
        H_Q=Q_carga/Q_fuente
        
        fig, axs = plt.subplots(2, 1)
        axs[0].semilogx(f, H_P)
        axs[1].semilogx(f, H_Q)
        
        axs[0].set_xticks((np.arange(70,110,step=10)))
        axs[1].set_xticks((np.arange(70,110,step=10)))
        
        axs[0].text(f[nMax],H_P[nMax]\
            ,f'({f[nMax].round(3)} [{unidadesFrecuencia}],{H_P[nMax].round(3)} [W/W])',conf_left)
        axs[1].text(f[nMax],H_Q[nMax]\
            ,f'({f[nMax].round(3)} [{unidadesFrecuencia}],{H_Q[nMax].round(3)} [VAR/VAR])',conf_left)
        
        axs[0].set_xlabel(f'frecuencia [{unidadesFrecuencia}]')
        axs[1].set_xlabel(f'frecuencia [{unidadesFrecuencia}]')
        
        axs[0].set_ylabel('Pout/Pin [W/W]')
        axs[1].set_ylabel('Pout/Pin [VAR/VAR]')

        axs[0].grid(True)
        axs[1].grid(True)
        plt.show()


if __name__ == '__main__':
    tool=Tool()
    
    #Condiciones de diseño

    RL=220 # Resistencia de carga [Ohmios]
    Rs=50 # Resistencia de la fuente [Ohmios]

    tipoRed='T' # tipoRed='': red tipo pi y tipoRed='T' red tipo T 
    Xs1_m=1   #Xs1: -1 si Xs1 es capacitor y 1 si Xs1 es inductor (Tipo=Entero)
    Xs2_m=1    #Xs2: -1 si Xs2 es capacitor y 1 si Xs2 es inductor (Tipo=Entero)


    f_o=103.9e6 # Frecuencuencia central [Hz]
    Q=f_o/200e3 # Factor de calidad
    #Q=10

    ## Condiciones de simulación

    f_min=70e6 # Frecuencia mínima de simulación [Hz] 
    f_max=f_o**2/f_min
    puntosDecada=100000 # Puntos por década
    resultados='resultados.raw'
    circuito='circuito.cir'
    tiempoSimulacion=5 # Tiempo de simulación [s]
    
    
    # Se realizan los cálculos (No modificar)
    omega_o=f_o*2*math.pi
    
    Rv=tool.RvirtualTipoTAndPi(tipoRed,Rs,RL,Q)
    Coe=tool.coeficientesXsxXpxTipoTAndPi(tipoRed,Rv,Rs,RL,Q)
    netlist=tool.netlistTipoTAndPi(omega_o,Coe,Xs1_m,Xs2_m,tipoRed,f_min,f_max,puntosDecada,resultados)
    tool.escritura(netlist,circuito)
    
    tool.simular(circuito,tiempoSimulacion)
    tool.graficar(resultados,Rs,RL)