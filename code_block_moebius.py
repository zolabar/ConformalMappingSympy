#!/usr/bin/env python
# coding: utf-8


import sympy as sym
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
sns.set()


# Geometrical constants and prescribed velocity


relativeEccentrcity = 0.5
R2_ = 7.6
R1_ = 5
shift = (R2_-R1_)*relativeEccentrcity
u_R_ = 0.4


# symbols



x, y, z = sym.symbols('x,y,z', real=True)
R1, R2, b = sym.symbols('R1, R2, b', real=True)
R, a, x1, x2, l, delta = sym.symbols('R, a, x1,x2,l, delta', real=True)
u_R, xi, eta, rho, mu = sym.symbols('u_R, xi, eta, rho, mu', real=True)


# Moebius Transform


w = (z+sym.I*a)/(a*z+sym.I)


# scaling, outer radius to 1

w = w.subs(z, x/R2+sym.I*y/R2)


# Separation of real and imaginary part


xi_ = sym.simplify(sym.re(w))


eta_ = sym.simplify(sym.im(w))


# Constants from literature


# a and R from Churchill, Brown

a_ = (1+x1*x2+sym.sqrt((1-x1**2)*(1-x2**2)))/(x1+x2)
a_ = a_.subs(x2, (b-R1))
a_ = a_.subs(x1, (R1+b))
a_ = a_.subs(R1, R1/R2)
a_ = a_.subs(b, b/R2)

R_ = (1-x1*x2+sym.sqrt((1-x1**2)*(1-x2**2)))/(x1-x2)
R_ = R_.subs(x2, (b-R1))
R_ = R_.subs(x1, (R1+b))
R_ = R_.subs(R1, R1/R2)
R_ = R_.subs(b, b/R2)


# Velocity in concentric annulus in w-plane


# Velocity in concentric annulus in w-plane
uw = u_R*sym.log(sym.sqrt(xi**2+eta**2))/sym.log(R)
uwNum = uw.subs(R, R_)
uwNum = uwNum.subs(u_R, u_R_).subs(R1, R1_)
uwNum = uwNum.subs(R2, R2_).subs(b, shift)
uwNum = sym.lambdify((xi, eta), uwNum)


# Plotting in w-Plane
RNum = R_.subs(R1, R1_).subs(R2, R2_)
RNum = float(RNum.subs(b, shift))
rho_ = np.linspace(1, RNum, 200)
theta = np.linspace(0, 2*np.pi, 200)
Rho, Theta = np.meshgrid(rho_, theta)
Xi = Rho * np.cos(Theta)
Eta = Rho * np.sin(Theta)
fig, ax = plt.subplots(figsize=(9, 9))
ax.set_xlim(left=-RNum*1.1, right=RNum*1.1)
ax.set_ylim(bottom=-RNum*1.1, top=RNum*1.1)
plt.pcolor(Xi, Eta, uwNum(Xi, Eta), cmap='rainbow')


# Velocity in eccentric annulus in z-plane

# Velocity in eccentric annulus in z-plane

uzNum = uw.subs(xi, xi_).subs(eta, eta_)
uzNum = uzNum.subs(a, a_)
uzNum = uzNum.subs(R, R_)
uzNum = uzNum.subs(u_R, u_R_).subs(R1, R1_)
uzNum = uzNum.subs(R2, R2_).subs(b, shift)
uzNum = sym.lambdify((x, y), uzNum)


# Geometry creation and plotting

F = (R2_**2-R1_**2+shift**2)/(2*shift)
M = np.sqrt(F**2-R2_**2)
alpha = 0.5*np.log((F+M)/(F-M))
beta = 0.5*np.log((F-shift+M)/(F-shift-M))
yShift1 = float((M*sym.coth(alpha).evalf()))
X = np.linspace(-np.pi, np.pi, 200)
Y = np.linspace(alpha, beta, 200)
X, Y = np.meshgrid(X, Y)
zeta = X + 1j*Y
x_y = M*np.tan(zeta/2)
X = np.real(x_y)
Y = np.imag(x_y)-yShift1
fig, ax = plt.subplots(figsize=(9, 9))
plt.pcolor(X, Y, uzNum(X, Y), cmap='rainbow')


plt.show()
