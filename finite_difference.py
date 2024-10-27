#--------------------------------------------------
# Tool 8: Finite difference method
# Author: Nicolas Beaudoin
# Nov 9, 2023
#--------------------------------------------------

import numpy as np
import math


def S(x=0, y=0):
    #function for heat source. returns 0 if none.
    return 1000#100*x + 100*y


def HeatTransfer1D(n,x0,L,k,T0,TL):
    #rod heat transfer by finite difference, non insulated end
    delx = L/n

    #AT=b
    A =[]
    #T=[T1,T2,T3, ..., Tn]
    b=[]
    
    for i in range(0,n+2):
        # since T[i-1] -2T[i] +T[i+1] = -delx**2/k * S(x[i])

        #add to A
        A.append([])
        for j in range(0,n+2):
            if i==0:
                A[i] = [1]+(n+1)*[0]
                break
            elif i == n+1:
                A[i] = (n+1)*[0] + [1]
                break
                
            elif j == i-1:
                A[i].append(1)
            elif j == i:
                A[i].append(-2)
            elif j == i+1:
                A[i].append(1)
            else:
                A[i].append(0)

        #add to b        
        if i ==0:
            b.append(T0)    #boundary condition
        elif i == n+1:
            b.append(TL)    #boundary condition
        else:
            b.append(-delx**2/k*S(x0))

        x0 +=delx
        
    
    return np.linalg.solve(A,b)    #solv and return T[]
    
    

def HeatTransfer2D(step, L,H,k,T0j,Ti0,Tnj,Tin):
    #2D heat transfer on square plate
    # T0j,Ti0,Tnj,Tin are boundary conditions on 4 sides

    #
    # should have gone with iterative method instead...
    #

    x0 = 0
    y0 = 0

    nx = round(L/step)  #floor or round?
    ny = round(H/step)
    n = ny*nx

    
    #AT=b
    A =[]
    #T=[T11,T21,T31, ..., Tnxny]
    b=[]
    
    
    for i in range(n): #kinda confusing ij notation
        A.append([])
        for j in range(n):
            '''
            if j==i or j==i-2 or j==i-4 or j== i-6:
                A[i].append(1)
            elif j==i-3:
                A[i].append(-4)
            else:
                A[i].append(0)
            '''
            if j==i-3 or j==i-1 or j==i+1 or j==i+3:
                A[i].append(1)
            elif j==i:
                A[i].append(-4)
            else:
                A[i].append(0)
            
                

        #add to b
        # BC outside edge of plate at 20C
        '''
        if i==0:
            b.append(Ti0)
        elif i ==n+5:
            b.append(Tnj)
        else:
            b.append(-step**2/k*S(x0,y0))
        '''

        b.append(-step**2/k*S(x0,y0))
        
        '''
        if i%nx ==0:
            b.append(T0j)    #boundary condition
        elif i%nx == nx-1:
            b.append(Tnj)
        elif i<nx:      #<=???
            b.append(Ti0)
        elif i%nx == ny-1:
            b.append(Tin)    #boundary condition
        else:
            b.append(-step**2/k*S(x0,y0))
        '''
        
        x0 +=step
        if i%nx ==0:   #end of row, next line
            x0 = 0
            y0 +=step
        #print (x0,y0)
    #print(A)
    #print(b)
    return np.linalg.solve(A,b)    #solve and return T[]


# MAIN
#T = HeatTransfer1D(10,0,1,200,100,20)
T = HeatTransfer2D(2,10,10,100,20,100,20,100)

print(T)
