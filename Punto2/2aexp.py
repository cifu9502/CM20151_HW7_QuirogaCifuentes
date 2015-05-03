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

# <codecell>

#El siguiente codigo realiza el proceso de cadenas de markov para encontrar el mejor ajuste al caso lineal.
#El codigo usado es parecido al dado como muestra en clase , a excepcion de unas pequeñas modificaciones

#input y_obs -- datos del punto a observar


# <codecell>

#Es el analogo a la funcion findLin para la exponencial
def findExp(y_obs):
    c_walk = np.empty((0)) #this is an empty list to keep all the steps
    d_walk = np.empty((0))
    nu_walk = np.empty((0))
    sigma_walk = np.empty((0))
    
    l_walk = np.empty((0))
    
    c_walk = np.append(c_walk, np.random.random()*100 -50)
    d_walk = np.append(d_walk, np.random.random()*1 -0.5)
    nu_walk = np.append(nu_walk, np.random.random() + 100)
    sigma_walk = np.append(sigma_walk, np.random.random() +150)
    
    y_init = expo(tiempos, c_walk[0], d_walk[0], nu_walk[0],sigma_walk[0])
    l_walk = np.append(l_walk, likelihood(y_obs, y_init))
    
    
    for i in range(n_iterations):
        c_prime = np.random.normal(c_walk[-1], 1) 
        d_prime = np.random.normal(d_walk[-1], 0.1)
        nu_prime = np.random.normal(nu_walk[-1], 2) 
        sigma_prime = np.random.normal(sigma_walk[-1], 1)
    
        y_init =  expo(tiempos, c_walk[-1], d_walk[-1],nu_walk[-1] ,sigma_walk[-1])
        y_prime = expo(tiempos, c_prime, d_prime, nu_prime, sigma_prime)
        
        l_prime = likelihood(y_obs, y_prime)
        l_init = likelihood(y_obs, y_init)
        
        alpha = l_prime/l_init
        if(alpha>=1.0):
            c_walk  = np.append(c_walk,c_prime)
            d_walk  = np.append(d_walk,d_prime)
            nu_walk  = np.append(nu_walk,nu_prime)
            sigma_walk  = np.append(sigma_walk,sigma_prime)
            
            l_walk = np.append(l_walk, l_prime)
            
            
            #Eliminamos los datos que ya no necesitamos para mejorar la efectividad. Note que en este caso, el likelihood
            #del anterior es menor que el del siguiente paso por lo que podemos eliminarlos
            c_walk = np.delete(c_walk, -2)
            d_walk = np.delete(d_walk, -2)
            nu_walk = np.delete(nu_walk, -2)
            sigma_walk = np.delete(sigma_walk, -2)
            
            l_walk = np.delete(l_walk, -2)
            
        else:
            beta = np.random.random()
            if(beta<=alpha):
                c_walk = np.append(c_walk,c_prime)
                d_walk = np.append(d_walk,d_prime)
                nu_walk = np.append(nu_walk,nu_prime)
                sigma_walk = np.append(sigma_walk,sigma_prime)
                
                l_walk = np.append(l_walk, l_prime)
    #print l_walk            
                
    return c_walk, d_walk ,nu_walk, sigma_walk , l_walk

# <codecell>


#Los comentarios que se muestran abajo son codigos usados para graficar los datos y verificar su veracidad

#La siguiente funcion invoca findLin y halla el mejor a y b.
#Recibe los datos del punto a observar
def hallarMinExp(y_obs):
    c_walk, d_walk, nu_walk, sigma_walk, l_walk = findExp(y_obs)
    
    #plt.scatter(c_walk, sigma_walk)    
    #plt.show()
    
    #plt.scatter(c_walk, -np.log(l_walk))
    #plt.show()
    
    max_likelihood = np.argmax(l_walk)
    best_l = l_walk[max_likelihood]
    best_c = c_walk[max_likelihood]
    best_d = d_walk[max_likelihood]
    best_nu = nu_walk[max_likelihood]
    best_sigma = sigma_walk[max_likelihood]
    
    
    #plt.plot(tiempos,y_obs,'b' , tiempos, expo(tiempos, best_c, best_d,best_nu, best_sigma),'r'  )
    #plt.show()
    return best_c, best_d, best_nu, best_sigma, best_l

# <codecell>

hallarMinExp(data[1:,3,3])

# <codecell>

fh = open('exponencial.dat', 'w')
fh.write('X'+ '\t'+ 'Y' + '\t' + ' c'+'\t'+' c'+'\t'+ 'nu'+'\t'+'sigma'+'Likelihood' + '\n')
for i in range(0,10):
    for j in range(0,10):
        bestc, bestd, bestnu, bestsigma, bestl = hallarMinExp(data[1:,i,j])
        fh.write(repr(i)+ '\t'+ repr(j) + '\t' + repr(bestc)+'\t'+repr(bestd)+'\t'+repr(bestnu)+'\t'+repr(bestsigma)+'\t'+repr(bestl) + '\n')
fh.close()    

# <codecell>



