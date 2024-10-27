#---------------------
# Tool 1: interpolation
# Author: Nicolas Beaudoin
# 
#---------------------

import math
import matplotlib.pyplot as plt

#chose example function : cos(x) + sin(2x)

#differentiation of example function (for taylor series)
def df(x, n):
    #could use symbolic package, but for trig functions it is simple to do this way
    if (n%4 == 0):
        return math.cos(x) + 2**n * math.sin(2*x)    # df0 = cosx + sin(2x)  ,  df4 = cos(x) + 16sin(2x)
    elif (n%4 == 1):
        return -math.sin(x) + 2**n * math.cos(2*x)   # df1 = -sin(x) + 2cos(2x)  ,  ...
    elif (n%4 == 2):            
        return -math.cos(x) - 2**n * math.sin(2*x)   # df2 = -cos(x) - 4sin(2x) 
    elif (n%4 == 3):
        return math.sin(x) - 2**n * math.cos(2*x)    # df3 = sin(x) - 8cos(2x)
    else:
        print("you messed up in df()")


# Taylor series interpolation
def taylor(x,x0,n):
    prodT = 1
    facT = 1
    sumT = 0
    
    for i in range(1,n+1):  
        sumT += prodT*df(x0,i-1)/facT
        prodT *=(x-x0)
        facT *=i   

    return sumT


# Lagrange interpolation
def theBarn(xjs,yjs,x):
    sumL = 0
  
    for i in range(0,len(yjs)):
        prodL = 1
        
        for j in range(0,len(xjs)):
            if i!=j:
                prodL *=(x-xjs[j])/(xjs[i]-xjs[j])   
                
        sumL += prodL*yjs[i]
        
    return sumL


#-----------------------
# main
#-----------------------

# Taylor single point test
n = 13     # chosen order of polynomial
x = 1.75   #chosen new point x    
x0 = 0     #chosen reference point x0

print("taylor calculation: " + str(taylor(x,x0,n)))
print("actual value: "+ str(df(x,0)))


# Lagrange single point test
x = 1.75     #
xjs = []     #data points of Xj
yjs = []     #

for p in range(0,41):  #interval of 0.5 for 40 points
    xjs.append(0.5*p)   #populate xjs
    
for p in xjs:
    yjs.append(df(p,0))  #populate yjs


print("lagrange calculation: " + str(theBarn(xjs,yjs,x)))
print("actual value: "+ str(df(x,0)))


# plot
xs = []   # points for calculation
Tys = []  # Taylor points
Lys = []  # Lagrange points
Fys = []  # Actual function points

for i in range(0,200):          
    xs.append(0.1*i)        #use consistent interval of 0.1 for points
    Tys.append(taylor(0.1*i,i//10,n))   
    Lys.append(theBarn(xjs,yjs,0.1*i))
    Fys.append(df(0.1*i, 0))  #df with n=0 is the function itself

plt.plot(xs,Tys,"g-")
plt.show()               #comment out to see all on same plot
plt.plot(xs,Lys,"r-")
#plt.show()               #comment out to see all on same plot
plt.plot(xs,Fys,"b.")
plt.show() 


