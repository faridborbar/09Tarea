
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import (leastsq, curve_fit)

def func_modelo(params, x):
    '''
    Modelo que utilizó Hubble: V = Ho * D (Caso a)
    d = (1 / Ho) * v (Caso b)
    En este caso "x" será D
    '''
    Ho = params
    return Ho * x



def func_a_minimizar(params, x_data, y_data):
    '''
    Función chi-cuadrado o función error = v - Ho * d (Caso a)
    Función chi-cuadrado o función error = d - v / Ho (Caso b)
    Donde X e Y ocuparan los roles de D y V, o V y D, para el caso a y b respectivamente. Despues habrá que invertir (1/Ho) para el caso b.
    '''
    return (y_data - func_modelo(params, x_data))

def bootstrap(data):
    '''
    Simulacion de bootstrap para encontrar el intervalo de confianza al 95%. Toma como argumento un archivo de datos de dos columnas.
    '''

    N = data.shape[0]
    N_bootstrap = 10000
    H = np.zeros(N_bootstrap)
    for i in range(N_bootstrap):
        s = np.random.randint(low=0, high=N, size=N)
        datos_falsos = data[s][s]
        x = datos_falsos[:, 0]
        y = datos_falsos[:, 1]
        H1, S1 = leastsq(func_a_minimizar, a0, args=(x, y))
        casi_H2, S2 = leastsq(func_a_minimizar, a0, args=(y, x))
        H2 = (1 / casi_H2)
        Hpromedio = (H1 + H2) / 2
        #print Hpromedio

        #Segun esta distribucion, obtenemos todos los H_o_promedio y los almacenamos
        H[i] = Hpromedio
    H = np.sort(H)
    lim_inf = H[int(N_bootstrap * 0.025)]
    lim_sup = H[int(N_bootstrap * 0.975)]
    print "El intervalo de confianza al 95% es: [{}:{}]".format(lim_inf, lim_sup)


'''
Implementamos el codigo que carga los datos, calcula la optimizacion para H_o y posteriormente muestra los resultados
'''

datos = np.loadtxt("data/SNIa.dat", usecols=(1,2))
d = datos[:, 1] # Distancia [Mpc]
v = datos[:, 0] # Velocidad [km/s]

# Adivinanza para el valor de Ho, caso a
a0 = 1
# Minimizacion del chi-cuadrado caso a
resultado_caso_a = leastsq(func_a_minimizar, a0, args=(d, v))
print "Status para a: ", resultado_caso_a[1]
print "mejor fit para Ho, caso a: ", resultado_caso_a[0]
Ho_a = resultado_caso_a[0]

#Adivinanza para el valor de 1/Ho, caso b
a1 = 100
# Minimizacion del chi-cuadrado caso b
resultado_caso_b = leastsq(func_a_minimizar, a1, args=(v, d))
print "Status para b: ", resultado_caso_b[1]
# En este caso el parámetro que optimizamos fue 1/Ho, por lo tanto lo invertimos para obtener Ho
Ho_b = 1 / resultado_caso_b[0]
print "mejor fit para Ho, caso b: ", Ho_b
print

#para la alternativa simetrica promediamos los valores
Ho_promedio = (Ho_a + Ho_b) / 2
print Ho_promedio


# Se grafican los valores experimentales y el ajuste con la minimización
# de la función chi-cuadrado (usando el H_o optimo)
fig = plt.figure()
fig.clf()
ax1 = fig.add_subplot(111)

ax1.plot(d, v, '*', label="Datos experimentales")
ax1.plot(d, Ho_promedio * d, label="Ajuste para $H_o$ optimo")

#ax1.set_xlim([-0.5, 2.5])
ax1.set_xlabel("Distancia $[Mpc]$")
ax1.set_ylabel("Velocidad $[km/s]$")
ax1.set_title("Grafico de distancia v/s velocidad")

plt.legend(loc='lower right')
plt.draw()
plt.show()

# Intervalo de confianza
interv_confianza = bootstrap(datos)
