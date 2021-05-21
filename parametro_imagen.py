import os  # Para saber la dirección absoluta
import math # Para cálculos matématicos
import subprocess # Para crear subprocesos
import ltspice # Extrae parámetros de archivos .raw
import matplotlib.pyplot as plt # Grafica en python de manera parecida a MatLab
import numpy as np # Permite hacer operaciones entre matrices o vectores


class Tools():

    def paralelo(self,a,b):
        '''
        Calcula la impedancia en paralelo
        '''
        return a*b/(a+b)

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

    def circuitoPasaBajas(self):
        print("El circuito pasabajas tiene la siguiente forma: ")
        c="\n        ---in--Lpi1_in---Lk1--nk--Lk2---Lm1--nm1--Lm2---Lpi1_out--out--\n\
           |                  |               |                    | \n\
           |                  |               |                    | \n\
        Lpi2_in               |              Lm3                Lpi2_out\n\
           |                  |               |                    | \n\
        np_in                Ck1             nm2                np_out \n\
           |                  |               |                    | \n\
        Cpi1_in               |              Cm1               Cpi1_out\n\
           |                  |               |                    | \n\
           |                  |               |                    | \n\
        ----------------------------------------------------------------gnd"

        print(c)
    
    def circuitoPasaAltas(self):
        print("El circuito pasa-altas tiene la siguiente forma: ")
        c="\n        ---in--Cpi1_in---Ck1--nk--Ck2---Cm1--nm1--Cm2---Cpi1_out--out--\n\
           |                  |               |                    | \n\
           |                  |               |                    | \n\
        Lpi2_in               |              Lm1                Lpi2_out\n\
           |                  |               |                    | \n\
        np_in                Lk1             nm2                np_out \n\
           |                  |               |                    | \n\
        Cpi2_in               |              Cm3               Cpi2_out\n\
           |                  |               |                    | \n\
           |                  |               |                    | \n\
        ----------------------------------------------------------------gnd"

        print(c)

    def makeNetlistPasaBajas(self,f_c,f_infty,z_o,f_min,f_max,n,guardarDatosEn):
        # Especificacion de parametros para filtro pasa-bajos por metodo de imagenes
        omega_c=f_c*2*math.pi # [rad/s]
        omega_infty=f_infty*2*math.pi # [rad/s]
        Z_N=z_o # Debe ser igual a la impedancia caracteristica


        # Desarrollo de calculos
        m=(1-(omega_c/omega_infty)**2)**(0.5) # Parametro m
        z_0=Z_N
            # Parametros de filtros-k
        C=2/(z_0*omega_c)
        L=z_0**2*C

        Ck1=C
        Lk1=Lk2=L/2
        print('Lk1= ',Lk1)
     

        # Circuito equivalente m
        Lm1=Lm2=m*L/2
        Lm3=((1-m**2)/(4*m))*L
        Cm1=m*C
        print('Lm1= ',Lm1)
        # Circuito compuesto

        mc=0.6 # Por el comportamiento de la impedancia normalizada de arreglo pi
        Lpi1_out=Lpi1_in=mc*L/2
        Lpi2_out=Lpi2_in=(1-mc**2)/(2*mc)*L
        Cpi1_out=Cpi1_in=mc*C/2
        print('Lpi1_out= ',Lpi1_out)


        netlist=".title filtro pasa bajas mediante parametro de imagen\n"\
        +"Vin in_1 0 dc 0 ac 1 "+"\n"\
        +"Rs in_1 in "+str(z_0)+"\n"\
        +"Lpi2_in in np_in "+str(Lpi2_in) +"\n"\
        +"Cpi1_in np_in 0 "+str(Cpi1_in) +"\n"\
        +"Lpi1_in_p_Lk1 in nk "+str(Lpi1_in+Lk1) +"\n"\
        +"Ck1 nk 0 "+str(Ck1)+ "\n"\
        +"Lk2_p_Lm1 nk nm1 "+str(Lk2+Lm1)+ "\n"\
        +"Lm3 nm1 nm2 "+str(Lm3)+ "\n"\
        +"Cm1 nm2 0 "+str(Cm1)+ "\n"\
        +"Lm2_p_Lpi1_out nm1 out "+str(Lm2+Lpi1_out) +"\n"\
        +"Lpi2_out out np_out "+str(Lpi2_out)+ "\n"\
        +"Cpi1_out np_out 0 "+str(Cpi1_out)+ "\n"\
        +"RL out 0 "+str(z_0)+"\n\n"\
        +".control\n"\
        +"ac dec "+str(n)+" "+str(f_min)+" "+str(f_max)+" \n"\
        +"*plot (vdb(out)-vdb(in))\n"\
        +"write "+os.getcwd()+'/'+guardarDatosEn+" all \n"\
        +".endc\n\n"\
        +".end"
        
        print(netlist)
        return netlist

    def makeNetlistPasaAltas(self,f_c,f_infty,z_o,f_min,f_max,n,guardarDatosEn):
        # Especificación de parámetros para filtro pasa-altas por método de imágenes
        omega_c=f_c*2*math.pi # [rad/s]
        omega=f_infty*2*math.pi # [rad/s]
        Z_N=z_o # Debe ser igual a la impedancia característica


        # Desarrollo de cálculos
        z_0=Z_N
            # Parámetros de filtros-k
        C=0.5/(z_0*omega_c)
        L=0.5*z_0/omega_c

        Ck2=Ck1=2*C
        Lk1=L
        print('C',C)
        print('Ck',Ck2)


        # Circuito equivalente m
        m=(1-(omega/omega_c)**2)**(0.5) # Parámetro m
        Cm1=Cm2=2*C/m
        Lm1=L/m
        Cm3=4*m*C/(1-m**2)
        print('Cm',Cm2)

        # Circuito compuesto

        mc=0.6 # Por el comportamiento de la impedancia normalizada de arreglo pi
        Cpi1_in=Cpi1_out=2*C/mc
        Lpi2_in=Lpi2_out=2*L/mc
        Cpi2_in=Cpi2_out=2*mc*C/(1-mc**2)
        print('c_pi',Cpi1_out)

        netlist=".title filtro pasa altas mediante parametro de imagen\n"\
        +"Vin in_1 0 dc 0 ac 1 "+"\n"\
        +"Rs in_1 in "+str(z_0)+"\n"\
        +"Lpi2_in in np_in "+str(Lpi2_in) +"\n"\
        +"Cpi2_in np_in 0 "+str(Cpi2_in) +"\n"\
        +"Cpi1_in_p_Ck1 in nk "+str(self.paralelo(Cpi1_in,Ck1)) +"\n"\
        +"Lk1 nk 0 "+str(Lk1)+ "\n"\
        +"Ck2_p_Cm1 nk nm1 "+str(self.paralelo(Ck2,Cm1))+ "\n"\
        +"Lm1 nm1 nm2 "+str(Lm1)+ "\n"\
        +"Cm3 nm2 0 "+str(Cm3)+ "\n"\
        +"Cm2_p_Cpi1_out nm1 out "+str(self.paralelo(Cm2,Cpi1_out)) +"\n"\
        +"Lpi2_out out np_out "+str(Lpi2_out)+ "\n"\
        +"Cpi2_out np_out 0 "+str(Cpi2_out)+ "\n"\
        +"RL out 0 "+str(z_0)+"\n\n"\
        +".control\n"\
        +"ac dec "+str(n)+" "+str(f_min)+" "+str(f_max)+" \n"\
        +"*plot (vdb(out)-vdb(in))\n"\
        +"write "+os.getcwd()+'/'+guardarDatosEn+" all \n"\
        +".endc\n\n"\
        +".end"
        
        print(netlist)
        return netlist

    def simular(self,archivo,t):
        '''
        Archivo: Netlist en .cir
        t: Tiempo de simulación en segundos
        '''
        process = subprocess.Popen(['ngspice',os.getcwd()+'/'+archivo])
        try:
            print('Running in process', process.pid)
            process.wait(timeout=t) #En segundos
        except subprocess.TimeoutExpired:
            print('Timed out - killing', process.pid)
            process.kill()
            print("Done")

    def graficar(self,archivo,filtro,f_infty):
        if filtro=='Circuito_Pasa_Bajas'.lower():
            titulo='Circuito pasa bajas con frecuencia de cero dB en '
        else:
            titulo='Circuito pasa altas con frecuencia de cero dB en '

        l = ltspice.Ltspice(os.getcwd()+'/'+archivo) 
        # Make sure that the .raw file is located in the correct path
        l.parse()
        f= l.get_frequency()
        V_in =abs(l.get_data('V(in)'))
        V_out = abs( (l.get_data('V(out)')))
        H=20*np.log10(V_out)-20*np.log10(V_in)
        f=np.multiply(f,1e-6)
        
        #plt.figure()
        plt.semilogx(f, H,f,np.ones(len(H))*-3,f,np.ones(len(H))*-10.3,f,np.ones(len(H))*-30)
        plt.legend(['H','-3dB','-10,3 dB','-30dB'])
        plt.ylabel('Magnitud [dB]')
        plt.xlabel('frecuencia [MHz]')
        #plt.title(titulo+str(f_infty)+' Hz')
        #plt.xlim([60,100])
        #plt.ylim([-2,2])
        plt.grid(True)
        plt.xticks(np.concatenate((np.arange(10,100,step=10),np.arange(100,300,step=100))))
        plt.savefig(os.getcwd()+'/'+filtro+'_'+str(f_infty)+'_Hz'+'.png')
        #plt.show()
        

    def tipoFiltro(self,tipo,f_c,f_infty,z_o,f_min,f_max,n,guardarDatosEn,circuito,tiempoSimulacion):
        '''
        tipo:
            - Escribir: 'Circuito_Pasa_Bajas' si es un pasa bajas
            - Escribir: 'Circuito_Pasa_Altas' si es un pasa altas
        f_c: frecuencia de corte
        f_infty: frecuencia donde la constante de propagación en una sección T tiende a infinito.
        z_o: Impedancia nominal
        f_min: frecuencia mínima
        f_max: frecuencia máxima
        n: número de punto por década
        guardarDatosEn: Archivo para guardar los datos
        '''
        if tipo.lower()== 'Circuito_Pasa_Bajas'.lower():
            self.circuitoPasaBajas()
            netlist= self.makeNetlistPasaBajas(f_c,f_infty,z_o,f_min,f_max,n,guardarDatosEn)
        elif tipo.lower()== 'Circuito_Pasa_Altas'.lower():
            self.circuitoPasaAltas()
            netlist= self.makeNetlistPasaAltas(f_c,f_infty,z_o,f_min,f_max,n,guardarDatosEn)
        else: print('Elija otra opci\243n')

        self.escritura(netlist,circuito)
        self.simular(circuito,tiempoSimulacion)
        self.graficar(guardarDatosEn,filtro,f_infty)
        

if __name__ == '__main__':
    for step in range(91,100,1):
        
        # Se ingreasan los parámetros
        filtro='Circuito_Pasa_Bajas'
        f_c=90E6 #Frecuencia de corte [Hz]
        #f_infty=(step)*1E6 #frecuencia [Hz] donde la constante de propagación en una sección T tiende a infinito.
        f_infty=(112.5)*1E6 
        z_0=75 #Impedancia caracteristica [Omh]
        f_min=10e6 #Límite inferior para simular [Hz]
        f_max=300e6 #Límite superior para simular [Hz]
        n=10000 #Número de puntos por década
        tiempoSimulacion=3 # Tiempo en el que dura toda la simulación [s]

        guardarDatosEn="resultados.raw" #Archivo .raw
        circuito='pasa_bajas.cir' #Netlist para ngspice
        
        # No se modifica
        tools=Tools()
        tools.tipoFiltro(filtro,f_c,f_infty,z_0,f_min,f_max,n,guardarDatosEn,circuito,tiempoSimulacion)
        
    


  




