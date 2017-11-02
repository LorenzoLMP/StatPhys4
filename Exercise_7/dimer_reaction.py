# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:05:56 2017

@author: LorenzoLMP
"""

from pylab import *  
from scipy import * 
from scipy import special
from scipy import misc
from mpl_toolkits.mplot3d import Axes3D
#import scipy.signal as sig

def G(z, alpha, beta):
    p1 = ((z+1)/2)**((1-2*beta)/2)
    p2 = special.iv(2*beta-1,sqrt(8*alpha*(z+1)))
    p3 = special.iv(2*beta-1,sqrt(16*alpha))
    return p1*p2/p3

#def H(z, ratio):
#    l = len(ratio);
#    alpha = 1*ones(l);
#    beta = alpha*ratio;
#    p1 = ((z+1)/2)**((1-2*ones(l)*beta)/2)
#    p2 = special.iv(2*ones(l)*beta-1,sqrt(8*ones(l)*alpha*(z+1)))
#    p3 = special.iv(2*ones(l)*beta-1,sqrt(16*ones(l)*alpha))
#    return p1*p2/p3

def H(z, alpha, beta):
    l = len(alpha);
#    alpha = 1*ones(l);
#    beta = alpha*ratio;
    p1 = ((z+1)/2)**((1-2*beta)/2)
    p2 = special.iv(2*beta-1,sqrt(8*ones(l)*alpha*(z+1)))
    p3 = special.iv(2*beta-1,sqrt(16*ones(l)*alpha))
    return p1*p2/p3

N = 1000
alpha = logspace(-3,1,N)
beta = logspace(-3,1,N)    
    
#mean_val = zeros((N,N))
#variance = zeros((N,N))
#for i in range(N):
#    mean_val[i][:] = misc.derivative(H, 1, dx=1e-10, n=1, args=(alpha,beta[i]))
#    variance[i][:] = sqrt(misc.derivative(H, 1, dx=1e-6, n=2, args=(alpha,beta[i])) + mean_val[i][:] - mean_val[i][:]**2)
#
#print(mean_val)
#
#
#fluct = divide(variance, mean_val)
#
#print('nonzero', nonzero(fluct[where( fluct<0.1 )]))
#print(fluct[where( fluct<0.1 )])

def Hmean(alpha, beta):
    return misc.derivative(G, 1, dx=1e-6, n=1, args=(alpha,beta))

def Hstd(alpha, beta):
    return misc.derivative(G, 1, dx=1e-6, n=2, args=(alpha,beta)) + Hmean(alpha, beta) -  Hmean(alpha, beta)**2  

def H_fluct(alpha, beta):
    return sqrt(Hstd(alpha, beta))/Hmean(alpha, beta)
    
#alpha, beta = meshgrid(alpha,beta)
#
#Axes3D.plot_surface(alpha, beta, fluct)
#def GG(z):
#    alpha = 1
#    beta = 1
#    return G(z,alpha, beta)
#
#def FF(z, *p):
#    return H(z,p[0])
#
#par = [1];
##misc.derivative(FF, 1, dx=1e-6, n=1)    
#    
#def Fmean():
#    return misc.derivative(H, 1, dx=1e-6, n=1)
#    
#def Fstd():
#    return misc.derivative(H, 1, dx=1e-6, n=2) + Gmean() -  Gmean()**2    
#    
#def Gmean():
#    return misc.derivative(GG, 1, dx=1e-6, n=1)
#    
#def Gstd():
#    return misc.derivative(GG, 1, dx=1e-6, n=2) + Gmean() -  Gmean()**2
#    
#def fluctuation():
#    return sqrt(Gstd())/Gmean()
#
#print('The mean of the distribution is: ', Gmean())
#print('The std of the distribution is: ', sqrt(Gstd()))
#print('The fluctuation is: ', fluctuation())
#
fig = figure()
ax = fig.add_subplot(111, projection='3d')
#ax = fig
X,Y = meshgrid(alpha, beta)
zs = array([H_fluct(x,y) for x,y in zip(ravel(X), ravel(Y))])
Z = zs.reshape(X.shape)

ax.plot_surface(log10(X),log10(Y),log10(Z), cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

ax.set_xlabel('log10(alpha)')
ax.set_ylabel('log10(beta)')
ax.set_zlabel('log10(fluctuation)')
ax.set_title(r'Surface plot of the fluctuations as function of $\alpha$, $\beta$')

savefig('fluctuations.png', dpi=600)
show()