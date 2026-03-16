from sympy import *

numberOfDimensions = 2
numberOfXi = 2
numberOfNodes = 4

dx, dy = symbols('dx,dy')
xi1, xi2  = symbols('xi1,xi2')
k = symbols('k')

C = zeros(numberOfNodes,numberOfNodes)
CE = zeros(numberOfNodes,numberOfNodes)

def Evaluate(Object):
    return Object.evalf(subs={dx:3/4,dy:1/4,k:1})

def PrintMatrix(Mat,Name):
    print("")
    for iIdx in range(0,numberOfNodes):
        for jIdx in range(0,numberOfNodes):
            print(" %s[%2d,%2d] = %s" % (Name,iIdx+1,jIdx+1,simplify(Mat[iIdx,jIdx])))
            
Phi = Matrix([[ (1-xi1)*(1-xi2) ],
                  [ xi1*(1-xi2) ],
                  [ (1-xi1)*xi2 ], 
                  [ xi1*xi2 ]])
J = dx*dy

for mIdx in range(0,numberOfNodes):
     for nIdx in range(0,numberOfNodes):
         integrand = k*Phi[mIdx]*Phi[nIdx]*J
         integral1 = integrate(integrand,(xi2,0,1))
         integral2 = simplify(integrate(integral1,(xi1,0,1)))
         C[mIdx,nIdx] = integral2
                 
CE = Evaluate(C)

PrintMatrix(C,'C')
PrintMatrix(CE,'CE')
