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
#Se obtienen los datose
data = pyfits.getdata(s,0)

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

#Se convierte el tiempo a minutos
tiempos = (tiempos-tiempos[0])*525600

# <codecell>


def lin(x_obs, a, b):
    return x_obs*b + a

def expo(x_obs, c, d, nu , sigma): 
    result = c +  x_obs*d + (1/(sigma*(np.pi)**(0.5)))*np.exp(-((x_obs -nu)/sigma)**2 /2)
    return result

def tan(x_obs, f, g,h,n,t0):
    return f+ x_obs*g + h*(1+(2/np.pi)*np.arctan(n*(x_obs -t0)))

# <codecell>

slin = 'lineal.dat'
sexp = 'exponencial.dat'
stan = 'tangente.dat'

infile1 = open(slin, 'r')
textlin = infile1.readlines()

infile2 = open(sexp, 'r')
textexp = infile2.readlines()
    
infile3 = open(stan, 'r')
texttan = infile3.readlines()

# <codecell>

fh = open('bestmodels.txt', 'w')

#Se itera sobre los 100 puntos y se grafica la recta correspondiente
for k in range(1,len(textlin)):
    lin1 = textlin[k].split('\t')
    exp1 = textexp[k].split('\t')
    tan1 = texttan[k].split('\t')
    
    lin1[-1] = lin1[-1].split('\n')[0]
    exp1[-1] = exp1[-1].split('\n')[0]
    tan1[-1] = tan1[-1].split('\n')[0]
          
    lin1[2:]= map(float, lin1[2:])
    exp1[2:]= map(float, exp1[2:])
    tan1[2:]= map(float, tan1[2:])
          
    lin1[0:2]= map(int, lin1[0:2])
    exp1[0:2]= map(int, exp1[0:2])
    tan1[0:2]= map(int, tan1[0:2])
    if((lin1[-1] >= exp1[-1])& (lin1[-1] >= tan1[-1]) ):
            plt.clf()
            plt.plot(tiempos,data[1:,lin1[0],lin1[1]],'b' , tiempos, lin(tiempos, lin1[2], lin1[3]),'r'  )
            plt.title('Lineal ('+repr(lin1[0])+','+repr(lin1[1])+')')
            plt.ylabel('B')
            plt.xlabel('t')
            plt.savefig('Lineal ('+repr(lin1[0])+','+repr(lin1[1])+')')
            fh.write(repr(lin1[0])+ ','+ repr(lin1[1]) + '\t' +'Lineal' +'\t'+repr(lin1[2])+'\t'+repr(lin1[3]) + '\n')
                
    elif((lin1[-1] <= exp1[-1])& (exp1[-1] >= tan1[-1]) ):
            plt.clf()
            plt.plot(tiempos,data[1:,lin1[0],lin1[1]],'b' , tiempos, expo(tiempos, exp1[2], exp1[3],exp1[4],exp1[5]),'r'  )
            plt.title('Exponencial ('+repr(lin1[0])+','+repr(lin1[1])+')')
            plt.ylabel('B')
            plt.xlabel('t')
            plt.savefig('Exponencial ('+repr(lin1[0])+','+repr(lin1[1])+')')
            fh.write(repr(lin1[0])+ ','+ repr(lin1[1]) + '\t' +'Exponencial' +'\t'+repr(exp1[2])+'\t'+repr(exp1[3])+ '\t'+repr(exp1[4])+ '\t'+repr(exp1[5]) + '\n')
    else:
            plt.clf()
            plt.plot(tiempos,data[1:,lin1[0],lin1[1]],'b' , tiempos, tan(tiempos, tan1[2], tan1[3],tan1[4],tan1[5],tan1[6]),'r'  )
            plt.title('Tangente ('+repr(lin1[0])+','+repr(lin1[1])+')')
            plt.ylabel('B')
            plt.xlabel('t')
            plt.savefig('Tangente ('+repr(lin1[0])+','+repr(lin1[1])+')')
            fh.write(repr(lin1[0])+ ','+ repr(lin1[1]) + '\t' +'Exponencial' +'\t'+repr(tan1[2])+'\t'+repr(tan1[3])+ '\t'+repr(tan1[4])+ '\t'+repr(tan1[5])+ repr(tan1[6]) + '\n')

fh.close()              
          

# <codecell>


