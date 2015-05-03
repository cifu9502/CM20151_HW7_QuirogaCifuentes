# -*- coding: utf-8 -*-
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


#Es el analogo a findlin para tan
def findTan(y_obs):
    f_walk = np.empty((0)) #this is an empty list to keep all the steps
    g_walk = np.empty((0))
    h_walk = np.empty((0))
    n_walk = np.empty((0))
    t0_walk = np.empty((0))
    l_walk = np.empty((0))
    
    f_walk = np.append(f_walk, np.random.random()*100 -50)
    g_walk = np.append(g_walk, np.random.random()*1 -0.5)
    h_walk = np.append(h_walk, np.random.random()*200 - 100)
    n_walk = np.append(n_walk, np.random.random()*100)
    t0_walk = np.append(t0_walk, np.random.random()*300 +80)
    
    y_init = tan(tiempos, f_walk[0], g_walk[0], h_walk[0],n_walk[0],t0_walk[0])
    l_walk = np.append(l_walk, likelihood(y_obs, y_init))
    
    
    for i in range(n_iterations):
        f_prime = np.random.normal(f_walk[-1], 1) 
        g_prime = np.random.normal(g_walk[-1], 0.1)
        h_prime = np.random.normal(h_walk[-1], 2) 
        n_prime = np.random.normal(n_walk[-1], 1)
        t0_prime = np.random.normal(t0_walk[-1], 5)
    
        y_init =  tan(tiempos, f_walk[-1], g_walk[-1],h_walk[-1] ,n_walk[-1],t0_walk[-1])
        y_prime = tan(tiempos, f_prime, g_prime, h_prime, n_prime,t0_prime)
        
        l_prime = likelihood(y_obs, y_prime)
        l_init = likelihood(y_obs, y_init)
        
        alpha = l_prime/l_init
        if(alpha>=1.0):
            f_walk  = np.append(f_walk,f_prime)
            g_walk  = np.append(g_walk,g_prime)
            h_walk  = np.append(h_walk,h_prime)
            n_walk  = np.append(n_walk,n_prime)
            t0_walk  = np.append(t0_walk,t0_prime)
            
            l_walk = np.append(l_walk, l_prime)
            
            
            #Eliminamos los datos que ya no necesitamos para mejorar la efectividad. Note que en este caso, el likelihood
            #del anterior es menor que el del siguiente paso por lo que podemos eliminarlos
            f_walk = np.delete(f_walk, -2)
            g_walk = np.delete(g_walk, -2)
            h_walk = np.delete(h_walk, -2)
            n_walk = np.delete(n_walk, -2)
            t0_walk = np.delete(t0_walk, -2)
            
            l_walk = np.delete(l_walk, -2)
            
        else:
            beta = np.random.random()
            if(beta<=alpha):
                f_walk = np.append(f_walk,f_prime)
                g_walk = np.append(g_walk,g_prime)
                h_walk = np.append(h_walk,h_prime)
                n_walk = np.append(n_walk,n_prime)
                t0_walk = np.append(t0_walk,t0_prime)
                
                l_walk = np.append(l_walk, l_prime)
    #print l_walk            
                
    return f_walk, g_walk ,h_walk, n_walk , t0_walk, l_walk

# <codecell>


#Los comentarios que se muestran abajo son codigos usados para graficar los datos y verificar su veracidad

#La siguiente funcion invoca findLin y halla el mejor a y b.
#Recibe los datos del punto a observar
def hallarMinTan(y_obs):
    f_walk, g_walk, h_walk, n_walk, t0_walk, l_walk = findTan(y_obs)
    
    #plt.scatter(f_walk, g_walk)    
    #plt.show()
    
    #plt.scatter(f_walk, -np.log(l_walk))
    #plt.show()
    
    max_likelihood = np.argmax(l_walk)
    best_l = l_walk[max_likelihood]
    best_f = f_walk[max_likelihood]
    best_g = g_walk[max_likelihood]
    best_h = h_walk[max_likelihood]
    best_n = n_walk[max_likelihood]
    best_t0 = t0_walk[max_likelihood]
    
    #plt.plot(tiempos,y_obs,'b' , tiempos, tan(tiempos, best_f, best_g,best_h, best_n, best_t0),'r'  )
    #plt.show()
    return best_f, best_g, best_h, best_n, best_t0 , best_l

# <codecell>

hallarMinTan(data[1:,3,3])

# <codecell>

fh = open('tangente.dat', 'w')
fh.write('X'+ '\t'+ 'Y' + '\t' + ' f'+'\t'+' g'+'\t'+ 'h'+'\t'+'n'+'\t'+'t0'+'\t'+'Likelihood' + '\n')
for i in range(0,10):
    for j in range(0,10):
        bestf, bestg, besth, bestn, bestt0, bestl = hallarMinTan(data[1:,i,j])
        fh.write(repr(i)+ '\t'+ repr(j) + '\t' + repr(bestf)+'\t'+repr(bestg)+'\t'+repr(besth)+'\t'+repr(bestn)+'\t'+repr(bestt0)+'\t'+repr(bestl) + '\n')
fh.close()    

# <codecell>


