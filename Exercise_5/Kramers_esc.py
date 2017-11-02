import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import scipy.integrate as integrate
from scipy import interpolate
from scipy.misc import derivative
from scipy import constants

#constants
c =299792458.  #m / s
m =1.10
pi = constants.pi
n_silica=1.46
n_water=1.33
P=0.05 #W
lam=532.0e-9 #m
w_0= 1.e-6 # 2*pi*c/lam m
R=50.0e-9 #m
k_b = 1.38064852e-23 #J/K

def omega(z):
	return w_0*np.sqrt( 1+((lam/n_water*z)/(pi*w_0**2))**2 )
	
def I(x,y,z):
	return 2*P/(pi*omega(z)**2) * np.exp(-2*(x**2+y**2)/omega(z)**2)

def V(x,y,z):
	return -2*pi*n_water*R**3*((m**2-1)/(m**2+2))*I(x,y,z)/c

def Vz(z):
	return -2*pi*n_water*R**3*((m**2-1)/(m**2+2))*I(0,0,z)/c

def V_norm(x,y,z):
	return V(x,y,z)/k_b
 
print(V_norm(0,0,0))
print(V_norm(0,0,5e-5))
print(derivative(Vz,0,n=2))
print(derivative(Vz,5e-5,n=2))

d_silica=2329.0
mass= d_silica*4*0.33*pi*R**3 
gamma = 0.8509e-3 #Pa/s
T = 300.0

r = 1/(2*pi)*np.sqrt(derivative(Vz,0,n=2)*derivative(Vz,5e-5,n=2))*np.exp( (V(0,0,5e-5)-V(0,0,0))/(k_b*T))/(gamma*mass)
print('r=',r)
 


print(V_norm(0,0,0))
#x = np.linspace(-3.e-6, 3.e-6, 100, endpoint=True)
#plt.figure(1)  
#plt.plot(x, V_norm(x,0,0))
#plt.xlabel('x (micrometer)')
#plt.ylabel('V(x)')
#plt.title('V(x)')
#plt.grid(True)
#plt.savefig("V(x).png")

rc('font', size=14)
xlabel(r' $\delta z$ [$m$]', size=16)
ylabel(r' Normalized V')
minorticks_on()
grid(which='major')
z = np.linspace(-50.e-6, 50.e-6, 100, endpoint=True)
plt.figure(1)  
plt.plot(z, V_norm(0,0,z))
#plt.xlabel('z (micrometer)')
#plt.ylabel('V(z)')
plt.title('Potential Trap V(z)')
plt.grid(True)
plt.savefig("V(z).pdf")


