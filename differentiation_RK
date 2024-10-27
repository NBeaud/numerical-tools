#--------------------------------------------------
# Tool 6: Runge-Kutta Differentiation
#
# Author: Nicolas Beaudoin
# October 19, 2023
#--------------------------------------------------

import math
import matplotlib.pyplot as plt

# funnel problem --------

def dh(t,h):
    g = 32.2
    alpha = 2/5
    Ae = 1/4
    return - math.sqrt(2*g)/(math.pi*alpha**2) * Ae * h**(-3/2)

def dt(h,t):
    g = 32.2
    alpha = 2/5
    Ae = 1/4
    return (math.pi*alpha**2)/(math.sqrt(2*g)*Ae) * h**(3/2)
 

# -----------------------

#callable function for R-K 4th order
#def RK4(t0,h0,delt,n):
    #tn = t0
    #hn = h0



# main -------------------
delt = 0.001
h1 = 20
h2 = 0

hn = h1
tn = 0

hs = [h1]
ts = [0]

        
while hn > h2:
    k1 = dh(tn,hn)
    k2 = dh(tn+delt/2, hn + abs(k1)*delt/2)
    k3 = dh(tn+delt/2, hn + abs(k2)*delt/2)
    k4 = dh(tn+delt, hn + abs(k3)*delt)

    hn += (k1+2*k2+2*k3+k4)/6 * delt
    tn += delt
    hs.append(hn)
    ts.append(tn)

print("t: " + str(tn))
print("h: " + str(hn))
plt.plot(ts, hs,".")
plt.show()


#inverse eq

delt = 0.001
h1 = 0
h2 = 179.28

hn = h1
tn = 0

hs = [h1]
ts = [0]

        
while hn < h2:
    k1 = dt(tn,hn)
    k2 = dt(tn+delt/2, hn + abs(k1)*delt/2)
    k3 = dt(tn+delt/2, hn + abs(k2)*delt/2)
    k4 = dt(tn+delt, hn + abs(k3)*delt)

    hn += (k1+2*k2+2*k3+k4)/6 * delt
    tn += delt
    hs.append(hn)
    ts.append(tn)

print("\nBelow obtained with dt/dh equation: ")
print("h: " + str(tn))
print("t: " + str(hn))
plt.plot(ts, hs,".")
plt.show()
