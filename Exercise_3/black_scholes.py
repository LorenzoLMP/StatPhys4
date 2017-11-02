# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 15:16:20 2017

@author: LorenzoLMP
"""
from pylab import *  
from scipy import * 
from scipy.stats import norm

def d1(x, c, s, r, tau, t):
    return (log(x/c) + (r + 0.5*s**2)*(tau-t))/(s*sqrt(tau-t))

def d2(x, c, s, r, tau, t):
    return (log(x/c) + (r - 0.5*s**2)*(tau-t))/(s*sqrt(tau-t))
    
def w(x, c, r, t, tau, d1, d2):
    return x*norm.cdf(d1) - c*exp(r*(t-tau))*norm.cdf(d2)


            

#DEFINITION OF PARAMETERS
c = 20
s = sqrt(0.2)
r = 1
tau = 1
t1 = 0.1
t2 = 0.5
t3 = 0.8
t4 = 0.99
    
rc('font', size=14)
xlabel(r' $X(t)$ [$ \$ $]')
ylabel(r'Option price $w(X,t)$ [$ \$ $]')
minorticks_on()
grid(which='major')
#yscale('log')
#xscale('log')
title("Option Price vs Stock Price")
xdata = linspace(0.1, 40, 30)
xxdata = linspace(0.1, 40, 200)
xxxdata = linspace(20, 40, 100)
plot(xdata, w(xdata, c, r, t1, tau, d1(xdata, c, s, r, tau, t1), d2(xdata, c, s, r, tau, t1)), linestyle='None', marker= 'o', color="red", label='t1=0.1')
plot(xdata, w(xdata, c, r, t2, tau, d1(xdata, c, s, r, tau, t2), d2(xdata, c, s, r, tau, t2)), linestyle='None', marker= '^', color="b",label='t2=0.5')
plot(xdata, w(xdata, c, r, t4, tau, d1(xdata, c, s, r, tau, t4), d2(xdata, c, s, r, tau, t4)), linestyle='None', marker= 'd', color="orange", label='t3=0.8')
plot(xdata, w(xdata, c, r, t3, tau, d1(xdata, c, s, r, tau, t3), d2(xdata, c, s, r, tau, t3)), linestyle='None', marker= 's', color="g",label='t4=0.99')


plot(xxdata, xxdata, linestyle='-',color="k", label='max. option')

plot(xxxdata, xxxdata-c, linestyle='--',color="k", label='min. option')

legend(loc='upper left')
savefig('black_scholes.png', dpi=600)
show()