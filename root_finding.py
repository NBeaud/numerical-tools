#--------------------------------------------------
# Assignment 3: Root finding
# Author: Nicolas Beaudoin
# September 28, 2023
#--------------------------------------------------

import matplotlib.pyplot as plt 
import math


TrueAns = 0.16672356243778486 #for comparison


#params
k1 = 40000
k2 = 40
m = 95
g = 9.81
h = 0.43
    
def f(d):       #returns 0 for proper value of d
    return (2/5)*k2*d**(5/2) + 0.5*k1*d**2 - m*g*d - m*g*h


def df1(d):      #first derivative
    return k2* d**(3/2) + k1*d -m*g


def df2(d):     #second derivative
    return (3/2)*k2*math.sqrt(d) + k1


# bisection method
def byBisect(xl,xu,error, n):
    plotPoints = []
    
    e = error+0.1     #have e greater than error just to go through loop
    xrold = 0
    count = 0
    
    while abs(e) > error and count < n:  #or... and... ?
        xr = (xl + xu)/2
        test = f(xl) * f(xr)

        if  test < 0 :  
            xu = xr

        elif test > 0 :
            xl = xr

        else:
            e = 0   #e just get rewriten, need break statement

        e = (xr - xrold)/xr #*100
        xrold = xr
        count +=1
        plotPoints.append(abs(e))

    plt.plot(plotPoints,"b")
    return xr    


# false position (regula falsi) method
def byFalsePos(xl,xu, error, n):
    plotPoints = []
    
    e = error +0.1     #have e greater than error just to go through loop
    xrold = 0
    count = 0
    
    while abs(e) > error and count < n:  
        xr = xu - f(xu)*(xl-xu)/(f(xl)-f(xu))
        test = f(xl) * f(xr)    #f(xl)*f(xu) ??? or just - f(xr)?

        if  test < 0 :
            xu = xr

        elif test > 0 :
            xl = xr

        else:
            e = 0

        e = (xr - xrold)/xr #*100
        xrold = xr
        count +=1
        plotPoints.append(abs(e))

    plt.plot(plotPoints,"r")
    return xr 


# Newton-Raphson method
def byNR(xi,ei, error, n):
    plotPoints = []
    
    e = ei     
    count = 0
    xr = xi
    
    while abs(e) > error and count < n:  
        xr = xr - f(xr)/df1(xr)   
        e = -df2(xr)/(2*df1(xr)) * e**2

        if df2(xr) == 0:
            print("f''(x) was equal zero. Error calculation may be affected")

        count +=1
        plotPoints.append(abs(e))
    
    plt.plot(plotPoints,"g")
    return xr 


# secant method
def bySecant(xi, xiold,error, n):
    plotPoints = []
       
    count = 0
    xr = xi
    xrold = xiold
    e = error +0.1
    
    while abs(e) > error and count < n:
        xrtemp = xr

        #if throws div/0 error, n is too large
        xr = xr - f(xr)*abs(xrold - xr)/abs(f(xr)-f(xrold))
        xrold = xrtemp
        e = xr - TrueAns

        count +=1
        plotPoints.append(abs(e))
    
    plt.plot(plotPoints,"k")
    return xr



#main ------------------------------------

roots = [byBisect(0,5,0.0,57), byFalsePos(0,5,0.0,546), byNR(5,2.5,0.0,10), bySecant(0,5,0.0,13)]
#root = byBisect(0,5,0.001,100)
#root = byFalsePos(0,5,0,100)
#root = byNR(5,2.5,0,10)
#root = bySecant(0,20,12)

for root in roots:
    print("d = " + str(root))
    print ("This should be ~0 if the root is correct: " + str(f(root))+"\n")
plt.show()
