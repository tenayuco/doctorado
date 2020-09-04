#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Apr 10 19:58:16 2020

@author: emilio

En este codigo esta todas las funciones para hacer pruebas simples con la ecuacion mono
y la ecuacion logistica. 
Sobre este codigo puedo hacer las pruebas temporales.
"""

from operator import add
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

"""
Las primeras tres funciones estan muy relacionadas pero las puse por separado
para poder llamar a algunas por separadao y tambien para dejar mayor claridad
de lo que estoy haciendo

"""
    


"""

esta funcion recibe los parametros (condicion inicial Zpre, los pasos a correr 
pasos, los parametros iniciales, el paso de integracion y el modo)
Llama a la funcion euler_red 

Y arroja una matriz con los valores actualizados de la ecuación para cada paso (va poniendo
los Z_t0, Z_t1, Zt2, etc )

"""
def creadorMatriz_sim(zpre, steps_, para_, size_step, mode_, dur_lluvia, dur_sol):
    pasos = steps_
    
    matrix_zvalues = np.zeros((pasos,len(zpre)))
    matrix_zvalues[0,:]= zpre
    D = 0
    
    t1 = np.arange(0, pasos, dur_lluvia+dur_sol)
    t2 = np.arange(0+dur_lluvia, pasos, dur_lluvia+dur_sol)
    transiciones = np.concatenate((t1, t2))
    for n in range(pasos-1):
        for k in transiciones:
            if n==k:
                D= abs(D-1)
            else:
                pass
        
        z1= euler_red(matrix_zvalues[n], para_, size_step, mode_, D)
        matrix_zvalues[n+1]= z1
    return(matrix_zvalues)
        
"""
esta funcion recibe los parametros, Z_t, el paso, el modo y genera un renglón de Z_t+h
Z_t+h = f(s,i,p)*h + Z_t
"""  

def euler_red(z_pre, para_, size_step, mode_, delta_): #z es para referir el vector sip que le metes
    h= size_step
    suma= [i * h for i in SIX_eu(z_pre, para_, mode_, delta_)]#aqui hago el calcylo
    znu= list(map(add,z_pre, suma))     # esto es la z más la suma de la diferencial pero por si quiero ver varios pasos en algun momento
    return(znu) #renglon


"""
esta funcion recibe el paso Z_t (z_pre), los parametros, el modo y genera (I_t+paso - I_t)/paso
la funcion aqui decide si es log, mono la funcion del crecimiento
"""
def SIX_eu(z_pre, para_, mode_, delta_):
   #parametros del modelo
    
    D = delta_
    para= para_
    z= z_pre
    
    r= para[0]
    b1=para[1]
    b2 = para[2]
    g= para[3]
    a= para[4]
    l= para[5]

    s= z[0]
    i= z[1]
    p= z[2]
    f= z[3]
    
    #ecuaciones
    if str(mode_)== "log":
        dsdt = r*s*(1-(s+i))-b1*s*p-b2*s*i  
    
    elif str(mode_) == "mono":
        dsdt = r*(1-(s+i))-b1*s*p-b2*s*i
    
    else: 
        print ("ERROR DE MODO")
    
    didt = b1*s*p+b2*s*i-g*i
    
    dpdt= a*i + l*( D*f-(1-D)*p)-p ## el menos P es esencial
    dfdt = l*((1-D)*p-D*f)#
    
    dzdt = [dsdt, didt, dpdt, dfdt]    
    return dzdt

"""
Esta funcion recibe los parametros, los pasos, el modo, el tamanio de integracion y los pasos
y hace un diagrama de fase con el campo vectorial (emplea el SIX) y una soluciones con diferentes 
iniciales 
"""



def matrizFases(parametros_ini, N, mode_, steps_, size_step):

    ###con esta primera parte se genera el campo vectorial
    ### se hace la rejilla en donde se va a resolver la ecuacion
    y1 = np.linspace(-0.2, 1.2, 20) 
    y2 = np.linspace(-0.2, 1.2, 20)
    Y1, Y2 = np.meshgrid(y1, y2)
    u, v = np.zeros(Y1.shape), np.zeros(Y2.shape)
    
    
    NI, NJ = Y1.shape
    
    
    r=parametros_ini[0]
    b1= parametros_ini[1]
    b2 = parametros_ini[2]
    g= parametros_ini[3]
    a= parametros_ini[4]
    
    

    for i in range(NI):
        for j in range(NJ):
            x = Y1[i, j]
            y = Y2[i, j]
            z_pre = [x, y, 0]
            
            yprime= SIX_eu(z_pre, parametros_ini, mode_) # se llama directo la ecuacion diferencial porque las flechas representan la derivada
           
            u[i,j] = yprime[0] # y aqui se van 
            v[i,j] = yprime[1]  


    fig = plt.figure(figsize=(10,10))
    Q = plt.quiver(Y1, Y2, u, v, color='r')

    plt.xlabel('$S$')
    plt.ylabel('$I$')
    plt.xlim([-0.5, 1.5])
    plt.ylim([-0.5, 1.5])
    plt.axhline()
    plt.axvline()
   
    """
    En esta parte se calcula algunas soluciones
    """
    pasos= steps_
        
    for p1 in [0.1, 0.5, 1]:
            z0 = [1, 0, p1]
            matriz = creadorMatriz_sim(z0, pasos, parametros_ini, size_step, mode_)
            plt.plot(matriz[:,0], matriz[:,1], 'c-')
            plt.plot(matriz[0,0], matriz[0,1], 'o')
            plt.plot(matriz[-1,0], matriz[-1,1], 's')   
        
    ro= (b2+b1*a)/g
    
    plt.title("ro= "+str(ro))
    #plt.show()
    plt.savefig('/home/emilio/python/figurasFase/pruebas/fig'+str(N)+'.png')
    plt.close()
    

"""
el graficador recibe una matriz (en el tiempo) de soluciones de S, I, P y arroja el valor de estas
en el tiempo
"""    
    
def graficador(matriz, dur_lluvia, dur_sol): #meterle un tuple
    X= matriz
    t = range(len(matriz))
    pasos= len(matriz)
# plot results
    
    fig, axs= plt.subplot(3, 1, figsize= (15, 3), facecolor= "w", edgecolor='k', sharey=True)
    fig.subplots_adjust(hspace = .5, wspace=.01)
    

    plt.plot(t,X[:,0],'g-',label="S", linewidth=4)
    plt.plot(t,X[:,1],'r-',label="I",linewidth=4)
    plt.plot(t,X[:,0]+X[:,1],'k-',label="T",linewidth=3)

    
    tran = np.arange(0, pasos, dur_lluvia+dur_sol)  #estos marcan los inicios de las lluvias (después de 1 ciclo de lluvoa y otro de sol)
    #t2 = np.arange(0+dur_lluvia, pasos, dur_lluvia+dur_sol)
        
    
    for k in tran:
        plt.axvspan(k, k+dur_lluvia, alpha=0.2, color='blue')

    plt.ylabel('Cantidad')
    plt.xlabel('Tiempo')

    plt.legend(loc='best')
    plt1= plt.show()
    
   
    plt.plot(t,X[:,2],"k-", label="Xact",linewidth=4)

    
    plt.plot(t,X[:,3],'y-',label="Einac",linewidth=4)
    for k in tran:
        plt.axvspan(k, k+dur_lluvia, alpha=0.2, color='blue')
    
    plt.legend(loc='best')
    plt.show()

    
    
    plt.plot(X[:,0],X[:,1],'g--',label="SvsI")
    plt.ylabel('I')
    plt.xlabel('S')
    plt.legend(loc='best')
    
    plt3= plt.show()
    
    return[plt1, plt3]
        
    
    """
    el graficador Recursivo hace diferentes graficas y las guarda
    """


def graficador2(matriz, dur_lluvia, dur_sol): #meterle un tuple
    X= matriz
    t = range(len(matriz))
    pasos= len(matriz)
# plot results
    
    fig, axs= plt.subplots(1, 3, figsize= (15, 4), facecolor= "w", edgecolor='k', sharey=False)
    fig.subplots_adjust(hspace = .05, wspace=.2)
    

    axs[0].plot(t,X[:,0],'g-',label="S", linewidth=4)
    axs[0].plot(t,X[:,1],'r-',label="I",linewidth=4)
    axs[0].plot(t,X[:,0]+X[:,1],'k-',label="T",linewidth=3)

    
    tran = np.arange(0, pasos, dur_lluvia+dur_sol)  #estos marcan los inicios de las lluvias (después de 1 ciclo de lluvoa y otro de sol)
    #t2 = np.arange(0+dur_lluvia, pasos, dur_lluvia+dur_sol)
        
    
    for k in tran:
        axs[0].axvspan(k, k+dur_lluvia, alpha=0.2, color='blue')

    axs[0].set_ylabel('Cantidad')
    axs[0].set_xlabel('Tiempo')
    axs[0].legend(loc="upper right")
    #axs[0].set_ylim([0,1])

    
    axs[1].plot(t,X[:,2],"k-", label="Xact",linewidth=4)
    axs[1].plot(t,X[:,3],'y-',label="Einac",linewidth=4)
    axs[1].set_ylabel('Cantidad')
    axs[1].set_xlabel('Tiempo')
    axs[1].legend(loc="upper right")
    #axs[1].set_ylim([0,4])


    for k in tran:
        axs[1].axvspan(k, k+dur_lluvia, alpha=0.2, color='blue')
   
    axs[2].plot(X[:,0],X[:,1],'g--',label="SvsI")
    axs[2].set_ylabel('I')
    axs[2].set_xlabel('S')
    axs[2].legend(loc="upper right")
    #axs[2].set_ylim([0,1])


    fig.savefig('/home/emilio/Desktop/figurasSimulaciones/estacionalidad/fig'+str(dur_lluvia)+"_"+str(dur_sol)+'.png')
    plt.close(fig)
   




def graficadorRec(matriz, N, ro_): #meterle un tuple
    X= matriz
    t = range(len(matriz))
    
    ro= ro_
# plot results
    plt.plot(t,X[:,0],'g-',label="S", linewidth=3)
    plt.plot(t,X[:,1],'r-',label="I",linewidth=3)
    
    
    plt.ylabel('Cantidad')
    plt.xlabel('Tiempo')
    plt.legend(loc='best')
    
    plt.title("ro= "+str(ro))
    plt.show()
    #plt.savefig('/home/emilio/Desktop/figurasSimulaciones/figurasTemporalesMono/fig'+str(N)+'.png')
    #plt.close()




"""
AQUI PARA CORRER TODO
"""
###script para correr la funcion una vez

"""
Primero para poder graficar las soluciones en el tiempo
"""

### el chiste es controlar los parametros como si fuera un año
#############CARACRTETRISTICAS GENERALES
steps = 400
sizeStep = 0.01
mode_= "mono"
lluvia = 30 #aquí recuerda que representa el optimo. 
sol= 100 #FALTA PODER JUGAR CON LA LLUVIA MAS Y CON LOS PARAMETROS #aqui puse el lag
#PARAMETROS 
b=1
b1=2
b2=6
a= 3
g= 6
l= 100
    
condIni = [1,0,0,4]

parametros =[b, b1,b2, g, a, l]


matrizSim = creadorMatriz_sim(condIni, steps, parametros, sizeStep, mode_, lluvia, sol)

ro= (b2+(b1*a))/(g) 
print ("ro=",ro)
graficador(matrizSim, lluvia, sol)


############3

for llu in np.arange(1, 151, 10):
    for sec in np.arange(1, 151, 10):
        #print("llu=",llu, "sec=", sec)
        matrizSim = creadorMatriz_sim(condIni, steps, parametros, sizeStep, mode_, llu, sec)
        graficador2(matrizSim, llu, sec)
    
    
    
    




#    
#
