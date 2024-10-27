#--------------------------------------------------
# Tool 5: Differentiation (Euler)
# Author: Nicolas Beaudoin
# October 20, 2023
#--------------------------------------------------

import math
import matplotlib.pyplot as plt

# funnel problem -------------------------

def dh(h):
    g = 32.2
    alpha = 2/5
    Ae = 1/4
    return - math.sqrt(2*g)/(math.pi*alpha**2) * Ae * h**(-3/2)

def nexth(hn, delt):
    return hn + delt*dh(hn)


# MAIN
delt = 0.1
h1 = 20
h2 = 0

hn = h1
t = 0

hs = [h1]
ts = [0]

while hn > h2:    #for i in range(0, n):
    hn = nexth(hn, delt)
    t+= delt
    hs.append(hn)
    ts.append(t)

print("t: " + str(t))
print("h: " + str(hn))

plt.plot(ts, hs,".")
plt.show()


# -----------------------------------------




# inverted pendulum -----------------------
# from moderncontrolengineering_ogata

M = 2
m = 0.1
l = 0.5
g = 9.81


def u(t,theta, thetadot):   #control force
    #return 1
    #return 25*theta                # P controller
    return 27.892*theta + thetadot  # PD controller
    #I to avoid deviation on long term # PID controller

    #Note:
    #frequency of correction depends on time step, thus changes D param
    #program could be made differet frquency of u than timestep


def x1dot(x2):
    return x2

def x2dot(x1,u):
    return ((M+m)*g*x1 -u)/(M*l)

def x3dot(x4):
    return x4

def x4dot(x1,u):
    return (-m*g*x1 +u)/M 


   
# MAIN

#params of numerical steps
n = 400         #nbr of iterations
delt = 0.1      #time step 

#starting values
theta0 = 0.1
thetadot0 = 0

#initializing storage vars
t = 0
thetan = theta0
thetadotn = thetadot0
x3n = 0
x4n = 0

#lists for graph
ts = []
thetas = []
xs = []
thetadots = []
xdots = []

#calculation loop
for i in range(0, n):
    #stepping
    thetadotn += delt * x2dot(thetan,u(t,thetan, thetadotn))
    thetan += delt * x1dot(thetadotn)
    x4n += delt * x4dot(thetan,u(t,thetan,thetadotn))
    x3n += delt * x3dot(x4n)
    t += delt

    #adding step to list for graph
    thetas.append(thetan)
    thetadots.append(thetadotn)
    xs.append(x3n)
    xdots.append(x4n)
    ts.append(t)
    
#output
print("theta: " + str(thetan))
print("thetadot: " + str(thetadotn))
print("x: " + str(x3n))
print("xdot: " + str(x4n))

plt.plot(ts, thetas,"y.")
#plt.show()
plt.plot(ts, thetadots,"g.")
#plt.show()
plt.plot(ts, xs,"b.")
#plt.show()
plt.plot(ts, xdots,"k.")
plt.show()

#------------------------------------------

