#---------------------
# Project: airfoil flow computation (incomplete)
#
# Author: Nicolas Beaudoin
# Date: December 16, 2023
#---------------------

import math
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------------
# Utility functions
#-------------------------------------------------------------

# Tridiagonal solver (used in cubic spline interpolation)
def triDiagSolve(A,b):
    #solves Ax=b
    for i in range(1, len(A)):  #len(A)-1?
        #decomposition
        A[i][0] /= A[i-1][1]
        A[i][1] -= A[i][0] * A[i-1][2]

        #forward substitution
        b[i] -= A[i][0]*b[i-1]
        #print(b[i])

    #back substitution
    x = [None]* len(A)    # initialize return vector
    
    x[len(A)-1] = b[len(A)-1] / A[len(A)-1][1]
    for i in range(len(A)-2,-1,-1):
        x[i] = (b[i] - A[i][2]*x[i+1]) / A[i][1]    

    return x



#Cubic spline interpolation for airfoil profile, determines second derivatives
def cubicSpline(xs,ys):
    #returns df2[] for known xs[], ys[]
    A = []
    b = []
    for i in range(len(xs)-2):
        #populate A row
        A.append([])
        for j in range(len(xs)-2):
            if (i==0 and j==0):  #first row, to keep 3 wide diagonal
                A[i].append(0)
            if j==i-1:
                A[i].append(xs[i+1]-xs[i])
            elif j==i:
                A[i].append(2*(xs[i+2]-xs[i]))
            elif j==i+1:
                A[i].append(xs[i+2]-xs[i+1])
            if(i==len(xs)-3 and j==len(xs)-3):  #last row, to keep 3 wide
                A[i].append(0)

        #populate b
        b.append(6*(ys[i+2]-ys[i+1])/(xs[i+2]-xs[i+1])  
                 + 6*(ys[i]-ys[i+1])/(xs[i+1]-xs[i]))   

    df2 = [0] + triDiagSolve(A,b) + [0]     #add zero at ends because f"(x)=0 at those points
    return df2                              #vector of second derivatives is returned



#quadratic spline profile for known second derivatives
def profile(x, df2, xs, ys):
    # the loop check if inputx>xi, then get proper fi() from found index i, then return eq
    i=0
    while i<len(df2)-1:
        i+=1
        if x < xs[i]:  
            break       
    
    return (df2[i-1]/(6*(xs[i]-xs[i-1])) *(xs[i]-x)**3
            + df2[i]/(6*(xs[i]-xs[i-1])) *(x-xs[i-1])**3
            + (ys[i-1]/(xs[i]-xs[i-1]) - df2[i-1]*(xs[i]-xs[i-1])/6) *(xs[i]-x)
            + (ys[i]/(xs[i]-xs[i-1]) - df2[i]*(xs[i]-xs[i-1])/6) *(x-xs[i-1]))



def profileByBisect(y,df2,xs,ys,xl,xu,error):
    #uses root finding to find interpolated value for  given y value
    
    #nested function for the equation to solve
    def f(x):
        return profile(x,df2,xs,ys)-y       #profile(x) = y
    
    e = error+1     #have e greater than error just to go through loop
    xrold = 0
    count = 0
    
    while abs(e) > error:  
        xr = (xl + xu)/2
        test = f(xl) * f(xr)
        #print(profile(xr,df2,xs,ys))
        #print(xr)
        if  test < 0 :
            xu = xr
        elif test > 0 :
            xl = xr
        else:
            e = 0
            break

        e = (xr - xrold)/xr #*100
        xrold = xr
        
    return round(xr,12)+10**-11   # managing rounding errors. not elegant                                    



def discreteProfile(xs,ys,df2,step):
    #returns points at intersections with the grid steps
    error = 0.000000000001*step          #root finding error relative to step
    pxs = []
    pys = []
    x=0
    y=0
    xr=0
    while x < xs[len(xs)-1]:
        pxs.append(x)
        pys.append(profile(x,df2,xs,ys))        #add x intercept
        
        yl = profile(x,df2,xs,ys)
        yr = profile(x+step,df2,xs,ys)
        if yl < yr:
            y = (yl//step) *step
        elif yl > yr:
            y = math.ceil(yl/step)*step
        while xr < x+step:
            if y < yr:
                y += step
            elif y >yr:
                y -= step

            xr = profileByBisect(y,df2,xs,ys,x,x+step,error)
            if xr>x+step: break
            pxs.append(xr)    
            pys.append(profile(xr,df2,xs,ys))   #add y intercepts

        x+=step
            
    return pxs,pys









#=====================================
#object for each square volume cell of the grid.

class Cell:
    #magnitude var, so dont have to cal every time? depends on algo

    def __innit__(self):  #i,j (grid pos) params/var?
        #self.rho = rho
        self.P = 0
        self.u = {"N":0,"S":0,"E":0,"W":0}
        self.v = {"N":0,"S":0,"E":0,"W":0}

    def solve():
        #might need adjacent cells, so make this ext funct
        pass

    def norm():    #magnitude of velocity vector
        return math.sqrt(self.x**2 + self.y**2)

    def __str__(self):
        return str(norm()) + f"({self.x}, {self.y}"

    # index pos in grid, tensor/all side info, what about partial elems on arfoil (mult sides)?
    # airfoil boundary elems just hold resultant, no sides.



def generateGrid(x,y,step):
    #creates X*Y matrix, divided in step*step sized elems. ceil X,Y to next step
    nx = math.ceil(x/step)     
    ny = math.ceil(y/step)
    #grid = [[None]*nx]*ny   #innitialize empty matrix of grid dimensions
    grid = []
    for j in range(ny):
        grid.append([])
        for i in range(nx):
            grid[j].append(Cell())   #create nx*ny grid filled with cell obj
    return grid


# boundary conditions
def setBC(grid,u,P):
    for j in range(len(grid)):
        grid[j][0].u["W"] = u                       #left
        grid[j][0].P = P #["W"] = P
        
        grid[j][len(grid)-1].u["E"] = u             #right
        grid[j][len(grid)-1].P =P #["E"] = P

        grid[j][0].P = P #["N"] = P                 #top
        grid[j][len(grid)-1].P = P#["S"] = P        #bottom

    return grid

'''
def populateGrid(rho,v,P,grid,xs,Uys,Lys,df2,step): #just be main at this point
    x0 = math.ceil(0.5*xs[len(xs)-1]/step)-1      #index of origin cell
    y0 = round((0.5*max(Uys)-min(Lys))/step)-1    #index of row of origin
    xptr = x0
    
    i=0
    j=y0
    while j>0:
        while i<xptr and i< len(grid[0]):
            if i==0:
                grid[j][i] = Cell(v,0,rho,P)      #source flow
                i+=1
                continue

            if i==xptr:
                #grid[j][i] = Cell(,rho)
                j+=1
                #if profile(step*(xptr-x0),df2,xs,Uys) < step*(i+2):
                break
                    
            else:
                lv = grid[j][i-1].v
                bv = grid[j-1][1].v
                grid[j][i] = Cell(lv[0]+bv[0],lv[1]+bv[1],rho)

    # yeah... you have to rethink this bud
                
'''

def populateGrid(grid,pxs,Upys,Lpys):
    #airfoil is starts at (xptr,yptr) grid index
    x0 = len(grid[0])//4        #pos x of airfoil head
    y0 = len(grid)//2           #pos y of airfoil head
    pptr = 1                    #position pointer in profile point list
    xptr = 0                    #ptr for x in grid
    yptr = y0                   #ptr for y in grid
    
    #for i in range(len(grid[0])):
        #for j in range(len(grid)//2):   #may leave last row unpopulated/out of range
            #if within profile space:
                   #continue    #cell stays 0
            #calcCell grid[y0+j][i]
            #calcCell grid[y0-j][i]



#---------------------------------------------------
# MAIN
#---------------------------------------------------

#airfoil profile datapoints
#NACA0012 ,alpha =0, p.330,417
xs = [0, 1.25, 2.5, 5.0, 7.5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 95, 100]
ys = [0, 1.894, 2.615, 3.555, 4.200, 4.683, 5.345, 5.737, 5.941, 6.002,
    5.803, 5.294, 4.563, 3.664, 2.623, 1.448, 0.807, 0] #0.126] need to force last point to 0   

#NACA0010-66 ,alpha =0
#xs = [0, 1.25, 2.5, 5.0, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100]
#ys = [0,1.489,2.011,2.656,3.089,3.400,3.856,4.178,4.578,
#      4.822,4.956,5.0,4.889,4.3,2.833,1.656,0]

#xs = [i * 1 for i in xs]       #xs with scaling factor
Uys = [i * 1 for i in ys]       #upper half
Lys = [i * -1 for i in ys]      #lower half

step = 0.1    #0.1% of c


Udf2 = cubicSpline(xs,Uys)
Ldf2 = cubicSpline(xs,Lys)

gridHeight = math.ceil(6*max(ys))               #3* thickness of airfoil
gridLength = math.ceil(2*xs[len(xs)-1])         #2* cord
grid = generateGrid(gridLength,gridHeight,step)    

Upx, Upy = discreteProfile(xs,Uys,Udf2,step)
Lpx, Lpy = discreteProfile(xs,Lys,Ldf2,step)    #Lpx = Upx

#       12
#50     100     50
#       12


#plotting grid point profile
plt.plot(Upx,Upy,".")
plt.plot(Lpx,Lpy,".")
plt.grid()
plt.show()


'''
#plotting interpolated profile
pxs = []
pUys = []
pLys = []
for x in range(0,1000):
    pUys.append(profile(0.1*x,Udf2,xs,Uys))
    pLys.append(profile(0.1*x,Ldf2,xs,Lys))
    pxs.append(0.1*x)

plt.plot(xs,ys,".")
plt.plot(pxs,pUys,".")
plt.plot(pxs,pLys,".")
plt.show()
'''
