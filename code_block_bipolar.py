#!/usr/bin/env python
# coding: utf-8


import sympy as sym
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
sns.set()

"""Difference between symbols and expressions

symbols:
symbols in this code are SymPy symbols for algebraic or symbolic
manipulation.

expressions:
Expressions are terms assigned to variable. For example, if R1 is symbol,
R1_ denotes an expression, that may be numeric or not.
Before numerifying SymPy objects, the symbols are substituted by
corresponding expressions.    

"""

# Geometrical constants, pressure drop and prescribed velocity


# Geometrical constants, prescribed pressure drop and prescribed velocity
relativeEccentrcity = 0.5
R2_ = 7.6
R1_ = 5
shift = (R2_-R1_)*relativeEccentrcity
u_R_ = 0.4
dp_ = 10
l_ = 1.55
mu_ = 10.11


R1, R2, gamma, xi, eta, x, y = sym.symbols('R1 R2 gamma xi eta x y', real=True)
b, A, B, C = sym.symbols('b, A B C', real=True)
epsilon, kappa, alpha, beta = sym.symbols('epsilon kappa alpha beta', real=True)
c, M, F, Psi, u_R, mu, l, dp = sym.symbols('c M F Psi u_R mu l dp', real=True)
k, m, n = sym.symbols('k m n', integer=True)


#  The velocity


u = (dp/(mu*l))*M**2*(Psi+A*eta+B-(sym.cosh(eta)-sym.cos(xi))/(4*(sym.cosh(eta)+sym.cos(xi))))
u = u + (u_R/(beta-alpha))*(eta-alpha)


A_ = (sym.coth(alpha)-sym.coth(beta))/(2*(alpha-beta))


B_ = (beta*(1-2*sym.coth(alpha))-alpha*(1-2*sym.coth(beta)))/(4*(alpha-beta))


F_ = (R2**2-R1**2+b**2)/(2*b)


M_ = sym.sqrt(F**2-R2**2)


alpha_ = 0.5*sym.ln((F+M)/(F-M))


beta_ = 0.5*sym.ln((F-b+M)/(F-b-M))



summand = sym.exp(-n*beta)*sym.coth(beta)*sym.sinh(n*(eta-alpha))
summand = summand - sym.exp(-n*alpha)*sym.coth(alpha)*sym.sinh(n*(eta-beta))
summand = (sym.cos(n*xi)/(sym.sinh(n*(beta-alpha))))*summand
Psi_ = sym.Sum((-1)**n*summand, (n, 1, m))


gamma_ = (M_*sym.coth(alpha_.subs(M, M_))).subs(F, F_)
gamma_ = float(gamma_.subs(R2, R2_).subs(R1, R1_).subs(b, abs(shift)))



eta_ = (M**2 - 2*M*(y+gamma) + x**2 + (y+gamma)**2)
eta_ = (M**2 + 2*M*(y+gamma) + x**2 + (y+gamma)**2)/eta_
eta_ = sym.ln(eta_)/2


xi_ = -sym.atan2(2*M*x, (M**2 - x**2 - (y+gamma)**2))


alphaNum = alpha_.subs(M, M_).subs(F, F_).subs(R2, R2_).subs(R1, R1_)
alphaNum = float(alphaNum.subs(b, abs(shift)))

betaNum = beta_.subs(M, M_).subs(F, F_).subs(R2, R2_).subs(R1, R1_)
betaNum = float(betaNum.subs(b, abs(shift)))

cNum = M_.subs(F, F_).subs(R2, R2_).subs(R1, R1_)
cNum = float(cNum.subs(b, abs(shift)))

velW = u.subs(A, A_).subs(B, B_).subs(Psi, Psi_)
velW = velW.subs(alpha, alpha_).subs(beta, beta_)
velW = velW.subs(M, M_).subs(F, F_)
velW = velW.subs(R2, R2_).subs(R1, R1_).subs(b, abs(shift))
velW = velW.subs(gamma, gamma_).subs(u_R, u_R_).subs(mu, mu_)
velW = velW.subs(l, l_).subs(dp, 10**5*dp_)

uwNum = sym.lambdify((xi, eta), velW.subs(m, 100))

# Geometry creation and plotting
Xi = np.linspace(-np.pi, np.pi, 200)
Eta = np.linspace(alphaNum, betaNum, 200)
Xi, Eta = np.meshgrid(Xi, Eta)
fig, ax = plt.subplots(figsize=(9, 9))
plt.pcolor(Xi, Eta, uwNum(Xi, Eta), cmap='rainbow')


# Velocity in z-plane eccentric annulus in $x$ and $y$



velZ = u.subs(A, A_).subs(B, B_).subs(Psi, Psi_)
velZ = velZ.subs(eta, eta_).subs(xi, xi_)
velZ = velZ.subs(alpha, alpha_).subs(beta, beta_)
velZ = velZ.subs(M, M_).subs(F, F_)
velZ = velZ.subs(R2, R2_).subs(R1, R1_).subs(b, abs(shift))
velZ = velZ.subs(gamma, gamma_).subs(u_R, u_R_).subs(mu, mu_)
velZ = velZ.subs(l, l_).subs(dp, 10**5*dp_)


uzNum = sym.lambdify((x, y), velZ.subs(m, 100))

# Geometry creation and plotting
X = np.linspace(-np.pi, np.pi, 200)
Y = np.linspace(alphaNum, betaNum, 200)
X, Y = np.meshgrid(X, Y)
zeta = X + 1j*Y
x_y = cNum*np.tan(zeta/2)
X = np.real(x_y)
Y = np.imag(x_y)-gamma_
fig, ax = plt.subplots(figsize=(9, 9))
plt.pcolor(X, Y, uzNum(X, Y), cmap='rainbow')

plt.show()

