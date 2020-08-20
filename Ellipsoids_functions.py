"""
    HMATbook is a set of jupyter notebooks to calculate the homogeneous properties of heterogenous materials 
    Copyright (C) 2020  by Karim AÏT AMMAR, Enguerrand LUCAS, Julie DIANI 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Please if use or improve the software contact julie.diani@gmail.com

"""
"""
Useful functions for ellipsoid calculation
Uses fortran_tools which must have been compiled with f2py previously
"""
import numpy as np
from numpy import pi,cos,sin,arccos,arcsin
from numpy.random import random_sample
from numpy.linalg import inv
from numpy import dot
import random
import math
from fortran_tools import *

# Useful global values, must match the global values in classes_v*.py
nu_max = 0.49999999

def Comp3333_to_66 (G) : 
    "Transforms a 3x3x3x3 behaviour tensor G  into a 6x6 behaviour matrix F"
    F=np.zeros((6,6))
    for i in range(3):
        for j in range(3):
            F[i,j] = G[i,i,j,j]
            
        F[i,5]=(G[i,i,0,1]+G[i,i,1,0])/2.
        F[i,3]=(G[i,i,1,2]+G[i,i,2,1])/2. 
        F[i,4]=(G[i,i,2,0]+G[i,i,0,2])/2. 
        F[3,i]=(G[1,2,i,i]+G[2,1,i,i])/2. 
        F[4,i]=(G[0,2,i,i]+G[2,0,i,i])/2.
        F[5,i]=(G[0,1,i,i]+G[1,0,i,i])/2.

    F[4,4]=(G[0,2,0,2]+G[2,0,0,2]+G[0,2,2,0]+G[2,0,2,0])/4. 
    F[3,3]=(G[1,2,1,2]+G[2,1,1,2]+G[1,2,2,1]+G[2,1,2,1])/4.  
    F[5,5]=(G[0,1,0,1]+G[1,0,0,1]+G[0,1,1,0]+G[1,0,1,0])/4.  
    F[4,3]=(G[0,2,1,2]+G[2,0,1,2]+G[0,2,2,1]+G[2,0,2,1])/4.  
    F[4,5]=(G[0,2,1,0]+G[2,0,1,0]+G[0,2,0,1]+G[2,0,0,1])/4.  
    F[3,4]=(G[1,2,0,2]+G[2,1,0,2]+G[1,2,2,0]+G[2,1,2,0])/4.  
    F[5,4]=(G[0,1,0,2]+G[1,0,0,2]+G[0,1,2,0]+G[1,0,2,0])/4.  
    F[3,5]=(G[1,2,1,0]+G[2,1,1,0]+G[1,2,0,1]+G[2,1,0,1])/4.   
    F[5,3]=(G[0,1,1,2]+G[1,0,1,2]+G[0,1,2,1]+G[1,0,2,1])/4. 
    
    return F

def Comp66_to_3333(F) : 
    "Transforms a 6x6 matrix F into a tensor 3x3x3x3 G"
    G = np.zeros((3,3,3,3))
    for i in range(3) :
        for j in range(3) :
            G[i,i,j,j]=F[i,j]
       
        G[i,i,0,1]=F[i,5]
        G[i,i,1,2]=F[i,3]
        G[i,i,2,0]=F[i,4]
        G[0,2,i,i]=F[4,i]
        G[1,2,i,i]=F[3,i]
        G[0,1,i,i]=F[5,i]
        G[i,i,1,0]=F[i,5]
        G[i,i,2,1]=F[i,3]
        G[i,i,0,2]=F[i,4]
        G[2,0,i,i]=F[4,i]
        G[2,1,i,i]=F[3,i]
        G[1,0,i,i]=F[5,i]
        
    G[0,1,0,1]=F[5,5]
    G[0,1,0,2]=F[5,4]
    G[0,1,1,0]=F[5,5]
    G[0,1,1,2]=F[5,3] 
    G[0,1,2,0]=F[5,4]
    G[0,1,2,1]=F[5,3]

    G[0,2,0,1]=F[4,5]
    G[0,2,0,2]=F[4,4]
    G[0,2,1,0]=F[4,5]
    G[0,2,1,2]=F[4,3] 
    G[0,2,2,0]=F[4,4]
    G[0,2,2,1]=F[4,3]

    G[1,0,0,1]=F[5,5]
    G[1,0,0,2]=F[5,4]
    G[1,0,1,0]=F[5,5]
    G[1,0,1,2]=F[5,3] 
    G[1,0,2,0]=F[5,4]
    G[1,0,2,1]=F[5,3]

    G[1,2,0,2]=F[3,4]
    G[1,2,1,0]=F[3,5]
    G[1,2,1,2]=F[3,3] 
    G[1,2,2,0]=F[3,4]
    G[1,2,2,1]=F[3,3]

    G[2,0,0,1]=F[4,5]
    G[2,0,0,2]=F[4,4]
    G[2,0,1,0]=F[4,5]
    G[2,0,1,2]=F[4,3] 
    G[2,0,2,0]=F[4,4]
    G[2,0,2,1]=F[4,3]

    G[2,1,0,1]=F[3,5]
    G[2,1,0,2]=F[3,4]
    G[2,1,1,0]=F[3,5]
    G[2,1,1,2]=F[3,3] 
    G[2,1,2,0]=F[3,4]
    G[2,1,2,1]=F[3,3]
 
    return G 

def Rotation_matrices(n) : 
    'Create a matrix nx3x3 composed of n rotation matrix with two angles and one cos(angle) such as orientation are evenly distributed '
    Q = np.zeros((n,3,3))
    for i in range(n) : 
        psi = random.random()*2*np.pi
        theta = random.random()*2*np.pi
        phi = math.acos(random.random())
        Q[i,0,0]=cos(psi)*cos(theta)-cos(phi)*sin(theta)*sin(psi)
        Q[i,0,1]=sin(theta)*cos(psi)+cos(phi)*sin(psi)*cos(theta)
        Q[i,0,2]=sin(phi)*sin(psi)
        Q[i,1,0]=-sin(psi)*cos(theta)-sin(theta)*cos(phi)*cos(psi)
        Q[i,1,1]=cos(psi)*cos(phi)*cos(theta)-sin(theta)*sin(psi)
        Q[i,1,2]=cos(psi)*sin(phi)
        Q[i,2,0]=sin(phi)*sin(theta)
        Q[i,2,1]=-sin(phi)*cos(theta)
        Q[i,2,2]=cos(phi)
    
        for j in range(3) : 
            for k in range(3):
                if (abs(Q[i,j,k]) < 10**-6 ) :
                    Q[i,j,k] = 0.0
            
    return Q

def Isotropic_Stiffness_Matrix(K,G) :
    'Returns the stiffness matrix of an isotropic material'
    C = np.zeros((6,6))
    lamb = K - (2*G/3.)
    C[0,0]=2.*G+lamb
    C[1,1]=2.*G+lamb
    C[2,2]=2.*G+lamb

    C[3,3]=G
    C[4,4]=G
    C[5,5]=G

    C[0,1]=lamb
    C[0,2]=lamb
    C[1,2]=lamb
    C[1,0]=lamb
    C[2,1]=lamb
    C[2,0]=lamb
    
    return C

def Isotropic_Compliance_Matrix(E,nu) :
    'Returns the compliance matrix of an isotropic material'
    S = np.zeros((6,6))
    S[0,0]=1./E
    S[1,1]=1./E
    S[2,2]=1./E

    S[3,3]=2.*(1+nu)/E
    S[4,4]=2.*(1+nu)/E
    S[5,5]=2.*(1+nu)/E

    S[0,1]=-nu/E
    S[0,2]=-nu/E
    S[1,2]=-nu/E
    S[1,0]=-nu/E
    S[2,1]=-nu/E
    S[2,0]=-nu/E
    
    return S
    

def Isotropic_S(S) : 
    Eh = (1./3.) * (1/S[0,0]+1/S[1,1]+1/S[2,2]) 
    nuh = (min(- (1/6)*Eh*(S[0,1] + S[0,2] + S[1,2] + S[1,0] + S[2,0] + S[2,1]),nu_max))
    return Eh,nuh    

def Isotropic_C(C) : 
    lam=(1./6.)*(C[0,1]+C[0,2]+C[1,2]+C[1,0]+C[2,0]+C[2,1])
    G = (1./3.) *(C[3,3]+C[4,4]+C[5,5])
    Eh = G * (3*lam+2*G)/(lam+G)
    nuh = min(lam / (2* (lam+G)),nu_max)
    return Eh,nuh
             
def Clear_matrix2 (C) : 
    n = C.shape[0]
    for i in range(n) : 
        for j in range(n) :
            if C[i,j]<10**-6 : 
                C[i,j] = 0
                

def Eshelby_tensor(Axis,Cm,Sm) : 
    
    Cm3 = Comp66_to_3333(Cm)
    a0,a1,a2 = Axis
    IJV = np.array([[0,0],[1,1],[2,2],[1,2],[0,2],[0,1]])
    Nit = 40
    Ntop = Nit
    Mtop = Nit
    dphi = pi/(Ntop-1)
    dtheta = pi/(Ntop-1)
    A = np.zeros((6,6))
    B = np.zeros((6,6,Mtop))
    G = np.zeros((6,6,Ntop))
    E = np.zeros((6,6))
    
    # Green function integration for half of the ellipsoid
    for m in range(Mtop) : 
        phi = m*dphi
        for n in range(Ntop) : 
            theta = n*dtheta
            X = np.array([sin(theta)*cos(phi)/a0 , sin(theta)*sin(phi)/a1 , cos(theta)/a2])
            CXX = np.zeros((3,3))
            for i in range(3) :
                for j in range(3) :
                    for k in range(3) : 
                        for l in range(3) :
                            CXX[i,k] += Cm3[i,j,k,l]*X[j]*X[l]
            CXX = inv(CXX)
            for i in range(6) :
                for j in range(6) :                     
                    I1 = IJV[i,0]
                    J1 = IJV[j,0]
                    I2 = IJV[i,1]
                    J2 = IJV[j,1]
                    G[i,j,n] = 0.5 * sin(theta) * (CXX[I1,J1]*X[I2]*X[J2] + CXX[I2,J1]*X[I1]*X[J2] + CXX[I1,J2]*X[I2]*X[J1] + CXX[I2,J2]*X[I1]*X[J1])
        
        
        B[:,:,m] = 0.5 * dtheta * (G[:,:,0]+G[:,:,Ntop-1])
        for i in range(1,Ntop-1) : 
            B[:,:,m] +=  dtheta * G[:,:,i]

    A = 0.5*(B[:,:,0]+B[:,:,Ntop-1])* dphi/(4*pi)
    for i in range(1,Ntop-1) : 
         A += B[:,:,i]* dphi/(4*pi)  
    
    for i in range(6) : 
        for j in range(6) : 
            E[i,j]=A[i,0]*Cm[0,j]+A[i,1]*Cm[1,j]+A[i,2]*Cm[2,j] + 4* (A[i,3]*Cm[3,j]+A[i,4]*Cm[4,j]+A[i,5]*Cm[5,j]) 
    
    return E

def Fast_Eshelby_tensor(Axis,Cm,Sm) : 
    
    Cm3 = Comp66_to_3333(Cm)
    a0,a1,a2 = Axis
    IJV = np.array([[0,0],[1,1],[2,2],[1,2],[0,2],[0,1]])
    Nit = 40
    Ntop = Nit
    Mtop = Nit
    dphi = pi/(Ntop-1)
    dtheta = pi/(Ntop-1)
    A = np.zeros((6,6))
    B = np.zeros((6,6,Mtop))
    G = np.zeros((6,6,Ntop))
    E = np.zeros((6,6))
    
    # Integrate Green function on half of the ellipsoid
    for m in range(Mtop) : 
        phi = m*dphi
        for n in range(Ntop) : 
            theta = n*dtheta
            X = np.array([sin(theta)*cos(phi)/a0 , sin(theta)*sin(phi)/a1 , cos(theta)/a2])
            TCXX = np.zeros((3,3))
            for i in range(3) :
                for j in range(3) :
                    for k in range(3) : 
                        for l in range(3) :
                            TCXX[i,k] += Cm3[k,j,i,l]*X[j]*X[l] ## Transpose of CXX
             
            CXX = inversiont(TCXX)  ## Calculate the inverse of TCXX with fortran
            
            for i in range(6) :
                for j in range(6) :                     
                    I1 = IJV[i,0]
                    J1 = IJV[j,0]
                    I2 = IJV[i,1]
                    J2 = IJV[j,1]
                    G[i,j,n] = 0.5 * sin(theta) * (CXX[I1,J1]*X[I2]*X[J2] + CXX[I2,J1]*X[I1]*X[J2] + CXX[I1,J2]*X[I2]*X[J1] + CXX[I2,J2]*X[I1]*X[J1])
        
        
        B[:,:,m] = 0.5 * dtheta * (G[:,:,0]+G[:,:,Ntop-1])
        for i in range(1,Ntop-1) : 
            B[:,:,m] +=  dtheta * G[:,:,i]

    A = 0.5*(B[:,:,0]+B[:,:,Ntop-1])* dphi/(4*pi)
    for i in range(1,Ntop-1) : 
         A += B[:,:,i]* dphi/(4*pi)  
    
    for i in range(6) : 
        for j in range(6) : 
            E[i,j]=A[i,0]*Cm[0,j]+A[i,1]*Cm[1,j]+A[i,2]*Cm[2,j] + 4* (A[i,3]*Cm[3,j]+A[i,4]*Cm[4,j]+A[i,5]*Cm[5,j]) 
    
    return E 

def Matrix_to_vecteur(A) : 
    L = np.zeros(81)
    for i in range(3) : 
        for j in range(3) : 
            for k in range(3):
                for l in range(3) : 
                    L[27*i + 9*j + 3*k + l] = A[i,j,k,l]
    return L

def Vecteur_to_matrix(L) : 
    A = np.zeros((3,3,3,3))
    for i in range(3) : 
        for j in range(3) : 
            for k in range(3):
                for l in range(3) : 
                    A[i,j,k,l] = L[27*i + 9*j + 3*k + l]
    return A