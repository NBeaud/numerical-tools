#--------------------------------------------------
# Tool 9: Finite element method
# Author: Nicolas Beaudoin
# Nov 24, 2023
#--------------------------------------------------

import numpy as np

#utility
def S(x=0, t=0):     #heat source
    return 100*x#10

def ExtEffect1(x, x1, x2):
    N1 = (x2-x)/(x2-x1)
    # N2 = (x-x1)/(x2-x1)   #same thing so will ignore    
    return S(x)*N1

def Integral(a,b,n):
    #Using Trapezoid method
    step = (b-a)/n
    
    x1 = a
    x2 = a + step
    I = 0
    
    while x2 <= b:
        I += (x2-x1)* (ExtEffect1(x1,a,b) + ExtEffect1(x2,a,b))/2
        x1 = x2
        x2 += step

    return I


# Problem solver
def HeatTransfer1D(n,L,T0,Tn,): #add k param
    K = []
    F = []
    x0 = 0

    l= L/n
    N = n+1     #number of points vs number of elements

    a = 1/l                     # 0.4 in example 31.2 of textbook
    b = Integral(0,l,10)       #12.5 in example, 1 subdivision 


    #populat matrices
    for i in range(N):   #row
        K.append([])
        for j in range(N):   #column
            if (i==0 and j==0):
                K[i].append(1)
            elif (i==N-1 and j==N-1):
                K[i].append(-1)
            elif (j==0 or j==N-1):
                K[i].append(0)
            elif i == j:
                K[i].append(a+a)
            elif j==i-1 or j==i+1:
                K[i].append(-a)
            else:
                K[i].append(0)
        '''
        due to boundary conditions
        solution T will feature dt/dx at begin and end
        '''
        #populate F 
        if i==0:                
            F.append(b-a*T0)      # because switched with dt/dx
        elif i==N-1:    
            F.append(b-a*Tn)      # because switched with dt/dx
        elif i==1:
            F.append(2*b + a*T0)
        elif i ==N-2:
            F.append(2*b + a*Tn)
        else:
            F.append(2*b)   #connectivity

    #print(K)
    #print(F)
    return np.linalg.solve(K,F)



# MAIN
#T = HeatTransfer1D(4,10,40,200)    # 4 elems, 10cm, k=1, 40C, 200C
T = HeatTransfer1D(100,1,-20,100)    # 100 elems, 1m, k=1, -20C, 100C
print(T)

