#--------------------------------------------------
# Tool 7: Tridiagonal Matrix Solver (Thomas algorithm)
# Author: Nicolas Beaudoin
# Nov 3, 2023
#--------------------------------------------------

import numpy as np #for comparison of results
import copy
import time


def TriDiagSolve(A,b):
    A = copy.deepcopy(A)
    b = copy.deepcopy(b)
    for i in range(1, len(A)):  #len(A)-1?
        #decomposition
        A[i][0] /= A[i-1][1]
        A[i][1] -= A[i][0] * A[i-1][2]

        #forward substitution
        b[i] -= A[i][0]*b[i-1]
        #print(b[i])

    #back substitution
    x = [None]* len(A)    # return vector
    
    x[len(A)-1] = round(b[len(A)-1] / A[len(A)-1][1],10)    #round to 10 decimals
    for i in range(len(A)-2,-1,-1):
        x[i] = round((b[i] - A[i][2]*x[i+1]) / A[i][1], 10) #round to 10 decimals

    return x

def findEigenVal(A,n, bn):  # bn = initial guess
     # limitations:
     #              bn not orthogonal to A

    #print("\nEigen values by power method : ")
    for i in range(n):
        bni = np.dot(A,bn)
        norm = np.linalg.norm(bni)
        eigval = norm   #/np.linalg.norm(bn)
        #print("Lambda " + str(i)+": "+ str(eigval))
        bn = bni/norm

    return eigval, bn



def populateA(n,e,f,g):
    #copies e,f,g for an nxn sized matrix
    A = []
    
    for i in range(0,n):
        if i==0:
            A.append([0,f,g])
        elif i==n-1:
            A.append([e,f,0])
        else:
            A.append([e,f,g])
    return A

def populateb(n):
    #creates vector of size n of [1,2,3,...,n]
    b=[]
    for i in range(0,n):
        b.append(i+1)
        
    return b


    
# MAIN ---------------------------------------------------------

#equation solver
'''
A = [[0,2,-1],
     [-1,2,-1],  #only store diagonals of matrix
     [-1,2,-1],  #omits 0s for memory saving 
     [-1,2,0]]   #represents tri diagonals [e,f,g]
'''

n = 4   #size of matrix

A = populateA(n,-1,2,-1)
b = populateb(n)

begin = time.time()
x = TriDiagSolve(A,b)
end = time.time()

print("Input tridiagonal A: ")
print(A)
print("\nInput vector b: ")
print(b)
print("\nSolution x: ")
print(x)
print("Elapsed time for tridiagsolve = " + str(end-begin))



#comparison with numpy.linalg.solve() function below

'''
otherA =[[2,-1,0,0],    #numpy takes normal matrix format (with all zeros)
         [-1,2,-1,0],      
         [0,-1,2,-1],
         [0,0,-1,2]]
'''

#repopulate A in normal matrix format (add zeros) 
newA=[]
for i in range(0,len(A)):
    newA.append([])
    for j in range(0,len(A)):
        if j==i-1:
            newA[i].append(A[i][0])
        elif j == i:
            newA[i].append(A[i][1])
        elif j == i+1:
            newA[i].append(A[i][2])
        else:
            newA[i].append(0)
A = newA 
being = time.time()    
verif = np.linalg.solve(A,b)  #takes normal matrix format
end=time.time()

print("\nExpected solution x: ")
print(verif)
print("Elapsed time for linalg.solve = " + str(end-begin))
print()


#Eigenvalues

                                # a) 
#A = [[7,4,1],[4,4,4],[1,4,7]]  # b)
#A = [[2,1,1],[4,-6,0],[-2,7,2]]
bn = np.random.rand(np.array(A).shape[1])  #initial random guess

eigval,eigvect = findEigenVal(A,10,bn)
eigvalmin,eigvectmin = findEigenVal(np.linalg.inv(A),10,bn)
eigvalmin = 1/eigvalmin #need inverse of returned eigenvalue of A^-1

expeigvals, expeigvects = np.linalg.eig(A)  #expected vals for comparison


print("\nInput matrix for eigenvalue calculation: ")
print(A)
print("\nEigenvector max by power method: ")
print(eigvect)
print("\nLambda max by power method: ")
print(eigval)
print("\nEigenvector min by inverse iteration: ")
print(eigvectmin)
print("\nLambda min by inverse iteration: ")
print(eigvalmin)
print("\nCondition number: ")
print(eigval/eigvalmin)
print("\nExpected eigenvalues:")
print(expeigvals)
print("\nExpected eigenvectors:")
print(expeigvects)






