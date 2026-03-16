from sympy import *

numberOfDimensions = 2
numberOfXi = 2
numberOfNodes = 4

dx, dy = symbols('dx,dy')
xi1, xi2  = symbols('xi1,xi2')
k = symbols('k')

K = zeros(numberOfNodes,numberOfNodes)
KE = zeros(numberOfNodes,numberOfNodes)

def Evaluate(Object):
    return Object.evalf(subs={dx:3/4,dy:1/4,k:-0.01})

def PrintMatrix(Mat,Name):
    print("")
    for iIdx in range(0,numberOfNodes):
        for jIdx in range(0,numberOfNodes):
            print(" %s[%2d,%2d] = %s" % (Name,iIdx+1,jIdx+1,simplify(Mat[iIdx,jIdx])))
            
dPhidxi = Matrix([[ -(1-xi2), -(1-xi1) ],
                  [   1-xi2 ,    -xi1  ],
                  [    -xi2 ,   1-xi1  ], 
                  [     xi2 ,     xi1  ]])
dPhimdx = Matrix([ 0, 0 ])
dPhindx = Matrix([ 0, 0 ])
gu = Matrix([[ 0, 0 ],
             [ 0, 0 ]])
J = dx*dy

dxidx = Matrix([[  1/dx,     0 ],
                [     0,  1/dy ]])

for xiIdx1 in range(0,numberOfXi):
    for xiIdx2 in range(0,numberOfXi):
        gu[xiIdx1,xiIdx2]=0
        for coordIdx in range(0,numberOfDimensions):
            gu[xiIdx1,xiIdx2] = gu[xiIdx1,xiIdx2] + dxidx[xiIdx1,coordIdx]*dxidx[xiIdx2,coordIdx]
for mIdx in range(0,numberOfNodes):
     for nIdx in range(0,numberOfNodes):
         integrand = 0
         for xiIdx1 in range(0,numberOfXi):
             for xiIdx2 in range(0,numberOfXi):
                 integrand = integrand + dPhidxi[mIdx,xiIdx1]*dPhidxi[nIdx,xiIdx2]*gu[xiIdx1,xiIdx2]
         integrand = k*integrand*J
         integral1 = integrate(integrand,(xi2,0,1))
         integral2 = simplify(integrate(integral1,(xi1,0,1)))
         K[mIdx,nIdx] = integral2
         
                 
KE = Evaluate(K)

PrintMatrix(K,'K')
PrintMatrix(KE,'KE')
                 

