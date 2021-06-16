# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 14:40:41 2021

@author: zolabar
"""

from sympy import *
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
sns.set()


# Geometrical constants and prescribed velocity
relativeEccentrcity=0.5
R2_=7.6
R1_=5
shift=(R2_-R1_)*relativeEccentrcity
u_R_=0.4
dp_=10
l_=1.55
mu_=10.11


R1, R2,gamma, xi, eta, x,y,b, A, B, C, epsilon, kappa, alpha, beta, c, M, F, Psi, u_R, mu, l, dp =symbols('R1 R2 gamma xi eta x y b A B C epsilon kappa alpha beta c M F Psi u_R mu l dp', real=True)
k, m, n = symbols('k m n', integer=True)


u=(dp/(mu*l))*M**2*(Psi+A*eta+B-(cosh(eta)-cos(xi))/(4*(cosh(eta)+cos(xi))))+(u_R/(beta-alpha))*(eta-alpha)


A_=(coth(alpha)-coth(beta))/(2*(alpha-beta))
B_=(beta*(1-2*coth(alpha))-alpha*(1-2*coth(beta)))/(4*(alpha-beta))

F_=(R2**2-R1**2+b**2)/(2*b)
M_=sqrt(F**2-R2**2)

alpha_=0.5*log((F+M)/(F-M))
beta_=0.5*log((F-b+M)/(F-b-M))


Psi_=Sum((-1)**n*(cos(n*xi)/(sinh(n*(beta-alpha))))*(exp(-n*beta)*coth(beta)*sinh(n*(eta-alpha))-exp(-n*alpha)*coth(alpha)*sinh(n*(eta-beta))), (n, 1, m))


yShift1=float((M_*coth(alpha_.subs(M,M_))).subs(F,F_).subs(R2,R2_).subs(R1,R1_).subs(b,abs(shift)))



# ### The actual conformal mapping itself

eta_=log((M**2 + 2*M*(y+gamma) + x**2 + (y+gamma)**2)/(M**2 - 2*M*(y+gamma) + x**2 + (y+gamma)**2))/2

xi_=-atan2(2*M*x,(M**2 - x**2 - (y+gamma)**2))


alphaNum=float(alpha_.subs(M,M_).subs(F,F_).subs(R2,R2_).subs(R1,R1_).subs(b,abs(shift)))
betaNum=float(beta_.subs(M,M_).subs(F,F_).subs(R2,R2_).subs(R1,R1_).subs(b,abs(shift)))


cNum=float(M_.subs(F,F_).subs(R2,R2_).subs(R1,R1_).subs(b,abs(shift)))



velW=u.subs(A,A_).subs(B,B_).subs(Psi,Psi_)
velW=velW.subs(alpha,alpha_).subs(beta,beta_).subs(M,M_).subs(F,F_)
velW=velW.subs(R2,R2_).subs(R1,R1_).subs(b,abs(shift)).subs(gamma,yShift1).subs(u_R,u_R_).subs(mu, mu_).subs(l,l_).subs(dp, 10**5*dp_)


uwNum=lambdify((xi,eta),velW.subs(m,100))

# Geometry creation and plotting
Xi = np.linspace(-np.pi, np.pi,200)
Eta = np.linspace(alphaNum,betaNum,200)
Xi, Eta = np.meshgrid(Xi, Eta)
fig, ax = plt.subplots(figsize=(9,9))
plt.pcolor(Xi,Eta, uwNum(Xi,Eta), cmap='rainbow')

velZ=u.subs(A,A_).subs(B,B_).subs(Psi,Psi_).subs(eta,eta_).subs(xi,xi_)
velZ=velZ.subs(alpha,alpha_).subs(beta,beta_).subs(M,M_).subs(F,F_)
velZ=velZ.subs(R2,R2_).subs(R1,R1_).subs(b,abs(shift)).subs(gamma,yShift1).subs(u_R,u_R_).subs(mu, mu_).subs(l,l_).subs(dp, 10**5*dp_)


uzNum=lambdify((x,y),velZ.subs(m,100))

# Geometry creation and plotting
X = np.linspace(-np.pi, np.pi,200)
Y = np.linspace(alphaNum,betaNum,200)
X, Y = np.meshgrid(X, Y)
zeta = X + 1j*Y
x_y = cNum*np.tan(zeta/2)
X = np.real(x_y)
Y = np.imag(x_y)-yShift1
fig, ax = plt.subplots(figsize=(9,9))
plt.pcolor(X,Y, uzNum(X,Y), cmap='rainbow')

plt.show()
