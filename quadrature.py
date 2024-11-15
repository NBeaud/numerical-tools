#--------------------------------------------------
# Tool 4: Numerical integration (Quadrature)
# Author: Nicolas Beaudoin
# October 5, 2023
#--------------------------------------------------


#import numpy as np
import math
import matplotlib.pyplot as plt

GaussW = [
    [],                                                             #0
    [0],                                                            #1
    [1,1],                                                          #2
    [0.5555556,0.8888889,0.5555556],                                #3
    [0.3478548,0.6521452,0.6521452,0.3478548],                      #4
    [0.2369269,0.4786287,0.5688889,0.4786287,0.2369269],            #5
    [0.1713245,0.3607616,0.4679139,0.4679139,0.3607616,0.1713245]   #6
    ]

GaussX = [
    [],
    [0],
    [-0.577350269,0.577350269],
    [-0.774596669,0,0.774596669],
    [-0.861136312,-0.339981044,0.339981044,0.861136312],
    [-0.906179846,-0.538469310,0,0.538469310,0.906179846],
    [-0.932469514,-0.661209386,-0.238619186,0.238619186,0.661209386,0.932469514]
    ]

def C(T):
    # Compute the heat required to raise 1kg of this material from -100 to 200C.
    return (0.132 + 1.56*10**(-4)* T + 2.64*10**(-7)* T**2) #Cp

def f(z):   #force function to integrate, for 0 to 30 range
    return 200*(z/(5+z))* math.e**(-2*z/30)

def M(z):   #moment function to integrate, for 0 to 30 range
    return 200*(z/(5+z))* math.e**(-2*z/30) *z


def Trapezoid(a,b,n):
    print("Using Trapezoid method...")
    step = (b-a)/n
    
    x1 = a
    x2 = a + step
    I = 0
    
    while x2 <= b:
        I += (x2-x1)* (C(x1) + C(x2))/2
        x1 = x2
        x2 += step

    return I


def Simpson(a,b,n,o=2):
    h = (b-a)/n
    I = 0
    
    if o == 3:
        # 3/8 rule
        print("Using Simpson 3/8 method...")
        x0 = a
        x1 = a + h
        x2 = a + 2*h
        x3 = a + 3*h

        for i in range(0,n,o):
            I += 3*h/8 *(C(x0) + 3*C(x1) + 3*C(x2) + C(x3))
            x0 = x3
            x1 += 3*h
            x2 += 3*h
            x3 += 3*h

    else:
        # 1/3 rule
        print("Using Simpson 1/3 method...")
        x0 = a
        x1 = a + h
        x2 = a + 2*h
        
        for i in range(0,n,2):
            I += h/3* (C(x0) + 4*C(x1) + C(x2))
            x0 = x2
            x1 += 2*h
            x2 += 2*h

    return I


def Gauss(a,b,n,o):
    print("Using Gauss quadrature...")
    
    if o<=1 or o>6:
        print("o is out of tabulated range")
        return 0
    
    h = (b-a)/n
    x1 = a
    x2 = x1 + h
    I=0

    for i in range(0,n):
        #variable change
        a0 = (x2+x1)/2
        a1 = (x2-x1)/2
        # new f(x) = f(a0+a1*xd) * a1
        
        #loop sum with tabulated vals       
        for j in range(0,o):
            #       Cj        *      f(xj)
            I += GaussW[o][j] * M(a0 + a1*GaussX[o][j])*a1

        x1=x2
        x2 += h    
        
    return I


# main ------------

#plot function
x = []
y = []
for p in range (0,30):
    x.append(p)
    y.append(M(p))

plt.plot(x,y)
plt.show()
    

# choose how to decunstruct f for easy computation
# diff segment diff methods
I = 0
I += Simpson(-100,200,300,3)
#Q = 1000* I     # m x integral

print("Numerically calculated value: " + str(I))
#print analytical integral when known for comparison
#print("Expected value (analytical): " + "42.732 kCal")
