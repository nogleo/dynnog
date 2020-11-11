# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 15:56:15 2020

@author: LVALEO
"""

import sympy as sm
import sympy.physics.mechanics as me
import matplotlib.pyplot as plt
import numpy as np
from pydy.system import System

q = me.dynamicsymbols('q0:12')
u = me.dynamicsymbols('u0:12')
qd = me.dynamicsymbols('q0:12', 1)
ud = me.dynamicsymbols('u0:12', 1)

m1, m2, g = sm.symbols('m1:3 g')
I = sm.symbols('I1:3x:z')
kg = sm.symbols('kgx:z')
ks = sm.symbols('ksx:z')
la = sm.symbols('la0:4x:z')
lb = sm.symbols('lb0:4x:z')
lc = sm.symbols('lc0:4x:z')
ld = sm.symbols('ld0:4x:z')

N0 = me.ReferenceFrame('N_0')

N1 = N0.orientnew('N_1', 'Body', (q[0], q[1], q[2]), '123')

N2 = N0.orientnew('N_2', 'Body', (q[3], q[4], q[5]), '123')          

P0 = me.Point('P_0')
P0.set_vel(N0, 0)

P1 = P0.locatenew('P_1', q[6]*N0.x + q[7]*N0.y + q[8]*N0.z)
P2 = P0.locatenew('P_2', q[9]*N0.x + q[10]*N0.y + q[11]*N0.z)
P1.set_vel(N0, u[6]*N0.x + u[7]*N0.y + u[8]*N0.z)
P2.set_vel(N0,  u[9]*N0.x + u[10]*N0.y + u[11]*N0.z)

Pa = [P0.locatenew('Pa_0', la[0]*N0.x + la[1]*N0.y + la[2]*N0.z),
      P0.locatenew('Pa_1', la[3]*N0.x + la[4]*N0.y + la[5]*N0.z),
      P0.locatenew('Pa_2', la[6]*N0.x + la[7]*N0.y + la[8]*N0.z),
      P0.locatenew('Pa_3', la[9]*N0.x + la[10]*N0.y + la[11]*N0.z)]

Pb = [P1.locatenew('Pb_0', lb[0]*N1.x + lb[1]*N1.y + lb[2]*N1.z),
      P1.locatenew('Pb_1', lb[3]*N1.x + lb[4]*N1.y + lb[5]*N1.z),
      P1.locatenew('Pb_2', lb[6]*N1.x + lb[7]*N1.y + lb[8]*N1.z),
      P1.locatenew('Pb_3', lb[9]*N1.x + lb[10]*N1.y + lb[11]*N1.z)]

Pc = [P1.locatenew('Pc_0', lc[0]*N1.x + lc[1]*N1.y + lc[2]*N1.z),
      P1.locatenew('Pc_1', lc[3]*N1.x + lc[4]*N1.y + lc[5]*N1.z),
      P1.locatenew('Pc_2', lc[6]*N1.x + lc[7]*N1.y + lc[8]*N1.z),
      P1.locatenew('Pc_3', lc[9]*N1.x + lc[10]*N1.y + lc[11]*N1.z)]

Pd = [P2.locatenew('Pd_0', ld[0]*N2.x + ld[1]*N2.y + ld[2]*N2.z),
      P2.locatenew('Pd_1', ld[3]*N2.x + ld[4]*N2.y + ld[5]*N2.z),
      P2.locatenew('Pd_2', ld[6]*N2.x + ld[7]*N2.y + ld[8]*N2.z),
      P2.locatenew('Pd_3', ld[9]*N2.x + ld[10]*N2.y + ld[11]*N2.z)]

loads = [(P1, -m1 * g * N0.y),
         (P2, -m2 * g * N0.y)]    

for i in range(4):
    Pa[i].v2pt_theory(P0, N0, N0)
    Pb[i].v2pt_theory(P1, N0, N1)
    Pc[i].v2pt_theory(P1, N0, N1)
    Pd[i].v2pt_theory(P2, N0, N2)
    loads.append((Pa[i], -(kg[0] * Pa[i].pos_from(Pb[i]).dot(N0.x) * N0.x + kg[1] * Pa[i].pos_from(Pb[i]).dot(N0.y) * N0.y + kg[2] * Pa[i].pos_from(Pb[i]).dot(N1.z) * N0.z)))
    loads.append((Pb[i], -(kg[0] * Pb[i].pos_from(Pa[i]).dot(N0.x) * N0.x + kg[1] * Pb[i].pos_from(Pa[i]).dot(N0.y) * N0.y + kg[2] * Pb[i].pos_from(Pa[i]).dot(N0.z) * N0.z)))
    loads.append((Pc[i], -(ks[0] * Pc[i].pos_from(Pd[i]).dot(N0.x) * N0.x + ks[1] * Pc[i].pos_from(Pd[i]).dot(N0.y) * N0.y + ks[2] * Pc[i].pos_from(Pd[i]).dot(N0.z) * N0.z)))
    loads.append((Pd[i], -(ks[0] * Pd[i].pos_from(Pc[i]).dot(N0.x) * N0.x + ks[1] * Pd[i].pos_from(Pc[i]).dot(N0.y) * N0.y + ks[2] * Pd[i].pos_from(Pc[i]).dot(N0.z) * N0.z)))

In1 = (me.inertia(N1, I[0], I[1], I[2]), P1)
In2 = (me.inertia(N1, I[3], I[4], I[5]), P2)

body1 = me.RigidBody('b_1', P1, N1, m1, In1)
body2 = me.RigidBody('b_2', P2, N2, m2, In2)

Pd[1].pos_from(P0)


kdes = [q[0].diff() - u[0],
        q[1].diff() - u[1],
        q[2].diff() - u[2],
        q[3].diff() - u[3],
        q[4].diff() - u[4],
        q[5].diff() - u[5],
        q[6].diff() - u[6],
        q[7].diff() - u[7],
        q[8].diff() - u[8],
        q[9].diff() - u[9],
        q[10].diff() - u[10],
        q[11].diff() - u[11]]

KM = me.KanesMethod(N0, q, u, kdes)
fr, frstar = KM.kanes_equations([body1, body2], loads)

consts = {I[0]: 1.0, I[1]: 1.0, I[2]: 1.0, I[3]: 1.0, I[4]: 1.0, I[5]: 1.0,
          kg[0]: 1.0, kg[1]: 1.0, kg[2]: 1.0,
          ks[0]: 1.0, ks[1]: 1.0, ks[2]: 1.0,
          g: 9.81,
          la[0]:  0.01,  la[1]:  0.01,  la[2]:  0.01,
          la[3]: -0.01,  la[4]:  0.01,  la[5]:  0.01,
          la[6]: -0.01,  la[7]:  0.01,  la[8]: -0.01,
          la[9]:  0.01, la[10]:  0.01, la[11]: -0.01,
          lb[0]:  0.01,  lb[1]: -0.01,  lb[2]:  0.01,
          lb[3]: -0.01,  lb[4]: -0.01,  lb[5]:  0.01,
          lb[6]: -0.01,  lb[7]: -0.01,  lb[8]: -0.01,
          lb[9]:  0.01, lb[10]: -0.01, lb[11]: -0.01,
          lc[0]:  0.01,  lc[1]:  0.01,  lc[2]:  0.01,
          lc[3]: -0.01,  lc[4]:  0.01,  lc[5]:  0.01,
          lc[6]: -0.01,  lc[7]:  0.01,  lc[8]: -0.01,
          lc[9]:  0.01, lc[10]:  0.01, lc[11]: -0.01,
          ld[0]:  0.01,  ld[1]: -0.01,  ld[2]:  0.01,
          ld[3]: -0.01,  ld[4]: -0.01,  ld[5]:  0.01,
          ld[6]: -0.01,  ld[7]: -0.01,  ld[8]: -0.01,
          ld[9]:  0.01, ld[10]: -0.01, ld[11]: -0.01,
          m1: 0.1,
          m2: 0.5
          }
ics = [0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0]


sys = System(KM)
sys.constants = consts
sys.times = np.linspace(0, 50.0, 1000000)
yi1 = consts[la[1]]-consts[lb[1]]
yi2 = yi1 + consts[lc[1]]-consts[ld[1]]
sys.initial_conditions = {q[0]: 0.0,  q[1]: 0.0,  q[2]: 0.0,
                          q[3]: 0.0,  q[4]: 0.0,  q[5]: 0.0,
                          q[6]: 0.0,  q[7]: yi1,  q[8]: 0.0,
                          q[9]: 0.0, q[10]: yi2, q[11]: 0.0,
                          u[0]: 0.0,  u[1]: 0.0,  u[2]: 0.0,
                          u[3]: 0.0,  u[4]: 0.0,  u[5]: 0.0,
                          u[6]: 0.0,  u[7]: 0.0,  u[8]: 0.0,
                          u[9]: 0.0, u[10]: 0.0, u[11]: 0.0}
sys.generate_ode_function(generator='cython')
S = sys.integrate()

plt.plot(sys.times, S)