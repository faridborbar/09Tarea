import numpy as np
import matplotlib.pyplot as plt


''' utilizamos el metodo de montecarlo para calcular intervalor de confianza, para este debemos calcular el ajuste lineal
tal que y = Mx + N'''


def montecarlo(flujo_i, error_i, flujo_z, error_z):

    N_mc = 10000
    a = np.zeros(N_mc)
    b = np.zeros(N_mc)
    for i in range(N_mc):
        r = np.random.normal(0, 1, size=len(flujo_i))
        muestra_i = flujo_i + error_i * r
        muestra_z = flujo_z + error_z * r
        a[i], b[i] = np.polyfit(muestra_i, muestra_z, 1)
    a = np.sort(a)
    b = np.sort(b)
    limite_bajo = a[int(N_mc * 0.025)]
    limite_alto = a[int(N_mc * 0.975)]
    limite_bajo_2 = b[int(N_mc * 0.025)]
    limite_alto_2 = b[int(N_mc * 0.975)]
    print "El intervalo de confianza al 95% de M es: [{}:{}]".format(limite_bajo,
                                                                limite_alto)
    print "El intervalo de confianza al 95% de N es: [{}:{}]".format(limite_bajo_2,
                                                                limite_alto_2)

# main

datos = np.loadtxt("data/DR9Q.dat", usecols=(80, 81, 82, 83))
#print datos
#pasamos los datos a unidad de "nmaggie"
flujo_i = datos[:, 0] * 3.631
error_i = datos[:, 1] * 3.631
flujo_z = datos[:, 2] * 3.631
error_z = datos[:, 3] * 3.631

#ahroa indentificamos los elementos del ajuste de grado 1: M y N
M, N = np.polyfit(flujo_i, flujo_z, 1)
print "ecuacion ajustada : y = {}x + {}".format(M, N)

# plot
plt.clf()
x = np.linspace(0, 500, len(flujo_i))
plt.plot(x, M*x + N, color='y', label='Ajuste')
plt.errorbar(flujo_i, flujo_z, xerr=error_i, yerr=error_z, fmt='*',
             label='Datos Observados')
plt.xlabel("Flujo Banda_I [10^{-6} Jy]")
plt.ylabel("Flujo Banda_Z [10^{-6} Jy]")
plt.title('Flujo Banda_I vs Banda_Z')
#plt.grid(True)

plt.savefig('p3.png')
plt.legend(loc=2)
plt.show()

# intervalo de confianza
intervalo = montecarlo(flujo_i, error_i, flujo_z, error_z)
