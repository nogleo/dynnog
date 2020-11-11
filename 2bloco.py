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
    loads.append((Pa[i], Pb[i].pos_from(Pa[i]).dot(N0.x) * N0.x + Pb[i].pos_from(Pa[i]).dot(N0.y) * N0.y + Pb[i].pos_from(Pa[i]).dot(N0.z) * N0.z))
    loads.append((Pb[i], Pa[i].pos_from(Pb[i]).dot(N1.x) * N1.x + Pa[i].pos_from(Pb[i]).dot(N1.y) * N1.y + Pa[i].pos_from(Pb[i]).dot(N1.z) * N1.z))

    