# -*- coding: utf-8 -*-


# SERIES DE FOURIER DISCRETAS
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

#Señal de ecg phystonet
file = open('ecg.txt','r')

openedFile = file.read().replace('\n',' ')

data_str = openedFile.split(' ')
data_num = list(map(int,data_str))

#Tamaño del registro
TR = len(data_num)
raw_signal = np.zeros(int(TR/2))
filt_signal = np.zeros(int(TR/2))

#separar señales a partir de los datos del registro
i = 0
j = 0
while (j< TR):
  raw_signal[i] = data_num[j]
  filt_signal[i] = data_num[j+1]
  j += 2
  i += 1
  
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(raw_signal)
plt.title("Señal Original ")
plt.xlabel("Muestras ")
plt.subplot(2,1,2)
plt.plot(filt_signal)
plt.title("Señal Filtrada ")
plt.xlabel("Muestras ")

#extraer un pulso de la señal
#extraer un latido (ciclo) del registro
signal = raw_signal[5800:6200]
#tamaño en muestras del pulso
M = len(signal)
#frecuencia de muestreo
fs = 500
#periodo de muestreo
ts = 1/fs
#vector de tiempos
t = np.arange(0,M*ts,ts)
#condición de periodicidad
#signal[0] = signal[-1]
signal[-1] = signal[0]

#grafica del pulso
plt.figure(2)
plt.plot(t,signal)
plt.title("Señal de ECG")
plt.xlabel("Tiempo (seg)")

### metodo de interpolación
#numero de nuevos puntos de la señal (numero impar)
Nn = 401
#Numero de trapezoides
Nt = Nn-1

#vector de tiempos
tr = np.linspace(0,t[-1],Nn)
#Aplica metodo de interpolación
inter = interpolate.splrep(t,signal)
yr = interpolate.splev(tr,inter)

#grafica
plt.figure(2)
plt.plot(tr,yr,'.r')

#Analisis de frecuencia
#Calculo del incremento del ángulo
S = (2*np.pi)/Nt
#Factor o escala de cada coeficiente
CTE = (1/np.pi)*S
#coeficiente a0 = offset
a0 = CTE*np.sum(yr[0:Nt])

#Número de componentes de frecuencia
N = int((Nn-1)/2)
#Frecuencia fundamental
fo = 1/t[-1]
#vector de angulos
tn = np.linspace(0,2*np.pi-S,Nn-1)
#ciclo para calcular an,bn,cn,pn,fn
an = np.zeros(N)
bn = np.zeros(N)
cn = np.zeros(N)
pn = np.zeros(N)
fn = np.zeros(N)
for n in range(1,N):
  an[n] = CTE*np.sum(yr[0:Nt]*np.cos(n*tn))
  bn[n] = CTE*np.sum(yr[0:Nt]*np.sin(n*tn))
  cn[n] = np.sqrt(an[n]**2+bn[n]**2)
  fn[n]= n*fo
  an[40] = 0
  cn[40] = 0
  
  #primer cuadrante
  if an[n] >= 0 and bn[n] >= 0:
      pn[n] = np.arctan(bn[n]/an[n])
   #primer cuadrante
  if an[n] < 0 and bn[n] >= 0:
      pn[n] = np.pi - np.arctan(bn[n]/an[n])

  #primer cuadrante
  if an[n] < 0 and bn[n] < 0:
      pn[n] = np.pi + np.arctan(bn[n]/an[n])
      
  #primer cuadrante
  if an[n] >= 0 and bn[n] < 0:
      pn[n] = -np.arctan(bn[n]/an[n])
 
#Asigna el coeficiente a0 a an y cn
an[0] = a0    
cn[0] = np.abs(a0)
#cn[40] = 0
#an[40] = 0
#bn[40] = 0
   
plt.figure(3)
plt.subplot(2,2,1)
plt.stem(fn,an)
plt.title("Espectro del coeficiente an ")
plt.xlabel("Frecuencia en hertz ")

plt.subplot(2,2,2)
plt.stem(fn,bn)
plt.title("Espectro del coeficiente bn ")
plt.xlabel("Frecuencia en hertz ")
  
plt.subplot(2,2,3)
plt.stem(fn,cn)
plt.title("Espectro de Magnitud ")
plt.xlabel("Frecuencia en hertz ")
  
plt.subplot(2,2,4)
plt.stem(fn,pn)
plt.title("Espectro de Fase ")
plt.xlabel("Frecuencia en hertz ")


#Tarea3 hacer la reconstrucción de la señal usando 
#f(t) = a0/2 + sum(ancos(nt) + bnsen(nt)) desde n=1 a inf
#donde t = np.linspace(0,2*pi-S,400) -> esto ya es igual a tn
#Reconstrucción de la señal
tt = np.linspace(0,2*np.pi-S,400)
#Arreglo para almacenar la señal reconstruida
ft = np.zeros(400)
#ft[0] = a0/2
#ciclo para evaluar f(t) = a0/2 + sum(ancos(nt) + bnsen(nt)) desde n=1 a inf
for n in range(1,N):
    if ft[n] == 0 :
        ft = a0/2
    ft = ft + (an[n]*np.cos(n*tt)+bn[n]*np.sin(n*tt))
plt.figure(4)
plt.title("Señal de ECG Reconstruida")
plt.xlabel("Tiempo (seg)")
plt.plot(t,signal)
plt.plot(t,ft,'g')    

#flitrar el ruido de 50Hz y hacer la reconstrucción
