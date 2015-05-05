#Para entender el codigo comleto de todo el punto 2 lo mejor es leer el makefile que da una vision global de lo que hace 
#cada parte del programa

# <nbformat>3.0</nbformat>

# <codecell>
from __future__ import division 
import pyfits
import numpy as np
import matplotlib.pyplot as plt

s = "hmi.m_45s.magnetogram.subregion_x1y1.fits"
n_iterations = 20000 #this is the number of iterations I want to make

# <codecell>

pyfits.info(s)

# <codecell>

 #Se obtienen los datos
data = pyfits.getdata(s,0)
data.shape

# <codecell>

#Se grafica una  linea de tiempo
test = data[:,138,241]
x = np.linspace(0, 397, 398)
y = np.linspace(0, 198, 199)
z = np.linspace(0, 206, 207)
#plt.plot(z, test)
#plt.show()

# <codecell>

#Se lee el archivo tiempos.csv que trae la distribucion de los tiempos
s = 'tiempos.csv'
infile = open(s, 'r')
tiempos = []
text = infile.readlines()
for x in text[1:]:
    #parte el texto x segun ','
    a = x.split(',')
    a = a[1].split('\n')
    tiempos.append(float(a[0]))
tiempos = np.array(tiempos)

#Se convierte el tiempo a minutos y se fija el tiempo inical en 0
tiempos = (tiempos-tiempos[0])*525600


# <codecell>

#plt.plot(tiempos,data[1:,3,3])

#plt.show()

# <codecell>

#Se define el likelihood
def likelihood(y_obs, y_model):
    chi_squared = (1.0/2.0)*sum((y_obs-y_model)**2)*10**(-5)
    return np.exp(-chi_squared)

# <codecell>

#Se define la funcion lineal al igual que en el enunciado
def lin(x_obs, a, b):
    return x_obs*b + a

# <codecell>

#Se define la funcion exponencial al igual que en el enunciado
def expo(x_obs, c, d, nu , sigma):
    if sigma <= 0:
        result = -10000
    else:    
        result = c +  x_obs*d + (1/((sigma*np.pi)**(0.5)))*np.exp(-((x_obs -nu)/sigma)**2 /2)
    return result

# <codecell>

#Se define la funcion tangente al igual que en el enunciado
def tan(x_obs, f, g,h,n,t0):
    return f+ x_obs*g + h*(1+(2/np.pi)*np.arctan(n*(x_obs -t0)))

# <codecell>

#El siguiente codigo realiza el proceso de cadenas de markov para encontrar el mejor ajuste al caso lineal.
#El codigo usado es parecido al dado como muestra en clase  a excepcion de unas pequenas modificaciones

#input yobs : datos del punto a observar
def findLin(y_obs):
    a_walk = np.empty((0)) #this is an empty list to keep all the steps
    b_walk = np.empty((0))
    l_walk = np.empty((0))
    
    a_walk = np.append(a_walk, np.random.random()*100 -50)
    b_walk = np.append(b_walk, np.random.random()*1 -0.5)
    
    y_init = lin(tiempos, a_walk[0], b_walk[0])
    l_walk = np.append(l_walk, likelihood(y_obs, y_init))
    
    
    for i in range(n_iterations):
        a_prime = np.random.normal(a_walk[-1], 1) 
        b_prime = np.random.normal(b_walk[-1], 0.1)
    
        y_init = lin(tiempos, a_walk[-1], b_walk[-1])
        y_prime = lin(tiempos, a_prime, b_prime)
        
        l_prime = likelihood(y_obs, y_prime)
        l_init = likelihood(y_obs, y_init)
        
        alpha = l_prime/l_init
        if(alpha>=1.0):
            a_walk  = np.append(a_walk,a_prime)
            b_walk  = np.append(b_walk,b_prime)
            l_walk = np.append(l_walk, l_prime)
            
            #Eliminamos los datos que ya no necesitamos para mejorar la efectividad. Note que en este caso, el likelihood
            #del anterior es menor que el del siguiente paso por lo que podemos eliminarlo
            a_walk = np.delete(a_walk, -2)
            b_walk = np.delete(b_walk, -2)
            l_walk = np.delete(l_walk, -2)
            
        else:
            beta = np.random.random()
            if(beta<=alpha):
                a_walk = np.append(a_walk,a_prime)
                b_walk = np.append(b_walk,b_prime)
                l_walk = np.append(l_walk, l_prime)
                
    return a_walk, b_walk , l_walk

# <codecell>


#Los comentarios que se muestran abajo son codigos usados para graficar los datos y verificar su veracidad, esa era mi form
#De comprobar si la actualizacion realizada tenia sentido

#La siguiente funcion invoca findLin y halla el mejor a y b.
#Recibe los datos del punto a observar
def hallarMins(y_obs):
    a_walk, b_walk, l_walk = findLin(y_obs)
    
    #plt.scatter(a_walk, b_walk)    
    #plt.show()
    
    #plt.scatter(a_walk, -np.log(l_walk))
    #plt.show()
    
    max_likelihood = np.argmax(l_walk)
    best_l = l_walk[max_likelihood]
    best_a = a_walk[max_likelihood]
    best_b = b_walk[max_likelihood]
    
    #plt.plot(tiempos,y_obs,'b' , tiempos, lin(tiempos, best_a, best_b),'r'  )
    #plt.show()
    return best_a, best_b, best_l

# <codecell>

hallarMins(data[1:,3,3])

# <codecell>

#Imprime el resultado sobre los pixeles de 10*10
fh = open('lineal.dat', 'w')
fh.write('X'+ '\t'+ 'Y' + '\t' + ' a'+'\t'+' b'+'\t'+ 'Likelihood' + '\n')
for i in range(0,10):
    for j in range(0,10):
        besta, bestb, bestl = hallarMins(data[1:,i,j])
        fh.write(repr(i)+ '\t'+ repr(j) + '\t' + repr(besta)+'\t'+repr(bestb)+'\t'+repr(bestl) + '\n')
fh.close()       



