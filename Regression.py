#--------------------------------------------------
# Tool 2: Regression by Least Squares Best Fit
# Author: Nicolas Beaudoin
# September 26, 2023
#--------------------------------------------------

import matplotlib.pyplot as plt
import math
import numpy as np


#linear function f of order o
def f(o, X, x): #order o, param array X, at value x
    sum = 0
    for i in range(0,len(X)):
        sum += X[i] * x**(o-i)
    return sum


def fit(x,y, o = 1, ForceOrigin =False):    #done with normal equations
    print(" Fit of order " + str(o) + " with ForceOrigin = " + str(ForceOrigin))
    n = len(x)
    if len(y) != n:
        print("error: y and x arrays are of different size")
        return

    # AX=B              
    A = []              # matrix of data xs of mixed order
    # X = [a,b,c,...]   # unknown param a,b,c,etc.
    B = y               # vector of ys of data

    #populate A
    for i in range(0,n):
        #populate row of A
        A.append([])
        if ForceOrigin:
            for j in range(0,o):
                A[i].append(strain[i]**(o-j))
        else:        
            for j in range(0,o+1):
                A[i].append(strain[i]**(o-j))

    #convert to np.array object        
    A = np.array(A)
    B = np.array(B)
    
    X = np.linalg.solve( np.matmul(A.T,A) , np.matmul(A.T,B) )
    print(X)   

    #-----
    # fit outputs below
    #-----
    
    #error sqared sum
    Sr=0
    for i in range(0,n):
        Sr+=stress[i]-(f(o,X,x[i]))
    print("Sr = "+ str(Sr))

    #y average
    yavg=0
    for val in y:
        yavg+=val
    yavg = yavg/n

    #St
    St = 0
    for i in range(0,n):
        St+=(y[i]-yavg)**2

    #standard error
    Sxy = math.sqrt(abs(Sr)/(n-2))
    print("Sx/y = " + str(Sxy))

    #Standard deviation
    Sy = math.sqrt(St/(n-1))
    print("Sy = " + str(Sy))

    #quality of fit
    r = math.sqrt((St-Sr)/St)
    print("r = " + str(r)+ "\n")

    #plot f()
    plotX = []
    plotY = []
    for i in range(0,100):  # can modify point spacing and range for different data
        plotX.append(i*0.0001)
        plotY.append( f(o, X, i*0.0001) )

    plt.plot(plotX,plotY)

    

#dataset, applied to stress strain curve's linear portion
#strain = [0.0020,   0.0045, 0.0060, 0.0013, 0.0085, 0.0005]   #x
#stress = [4965,     5172,   5517,   3586,   6896,   1241]     #y
strain = [4,8,5,39,19,17,21]
stress = [0.8,2,1.2,10,4.2,3.7,5.2]

#main
plt.plot(strain,stress,'.')
fit(strain,stress,2)
fit(strain,stress,2,True)
plt.show()






'''
# simpleton linear fit (avoids linear algebra)

#number of data points
n = len(strain)

#linear fit
Sumx=0
Sumy=0
Sumxy=0
Sumx2=0

for i in range(0,n):
    Sumx += strain[i]
    Sumy += stress[i]
    Sumxy += strain[i]*stress[i]
    Sumx2 += strain[i]**2


#from minimisation formulas derived in class
a = (n*Sumxy-Sumx*Sumy)/(n*Sumx2-Sumx**2)
b = (Sumy - a*Sumx)/n
print("Best fit equation: "+ str(a)+ 'x + '+ str(b))

'''







