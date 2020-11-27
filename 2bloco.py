import sympy as sm
import sympy.physics.mechanics as me
import matplotlib.pyplot as plt
import numpy as np
from pydy.system import System

#me.init_vprinting()

q = me.dynamicsymbols('q0:15')
u = me.dynamicsymbols('u0:15')
qd = me.dynamicsymbols('q0:15', 1)
ud = me.dynamicsymbols('u0:15', 1)
t = me.dynamicsymbols._t

m1, m2, m3, m4, m5, m6, g = sm.symbols('m_1:7 g')
I1 = sm.symbols('I_1_[x:z]')
I2 = sm.symbols('I_2_[x:z]')
I3 = sm.symbols('I_3_[x:z]')
I4 = sm.symbols('I_4_[x:z]')
I5 = sm.symbols('I_5_[x:z]')
I6 = sm.symbols('I_6_[x:z]')
kg = sm.symbols('kg_[x:z]')
ks = sm.symbols('ks_[x:z]')
cg = sm.symbols('cg_[x:z]')
cs = sm.symbols('cs_[x:z]')
l01 = sm.symbols('l_01_0:4_[x:z]')
l10 = sm.symbols('l_10_0:4_[x:z]')
l12 = sm.symbols('l_12_0:4_[x:z]')
l21 = sm.symbols('l_21_0:4_[x:z]')
l23 = sm.symbols('l_23_[x:z]')
l34 = sm.symbols('l_34_[x:z]')
l45 = sm.symbols('l_45_[x]')
l56 = sm.symbols('l_56_[x]')
#l6 = sm.symbols('l_6_[x]')

N0 = me.ReferenceFrame('N_0')
N1 = N0.orientnew('N_1', 'Body', (q[0], q[1], q[2]), '123')
N2 = N0.orientnew('N_2', 'Body', (q[3], q[4], q[5]), '123')
N3 = N2.orientnew('N_3', 'Axis', (q[12], N2.y))
N5 = N2.orientnew('N_5', 'Axis', (sm.asin(l45/l56*me.dot(N3.x, N2.x)), N2.y))
N6 = N2.orientnew('N_6', 'Body', (0, 0, 0), '123')

# Origin Point
P0 = me.Point('P_0')
P0.set_vel(N0, 0*N0.x + 0*N0.y + 0*N0.z)

# Center of Mass and Connection Points
P1 = P0.locatenew('P_1', q[6]*N0.x + q[7]*N0.y + q[8]*N0.z)
P2 = P0.locatenew('P_2', q[9]*N0.x + q[10]*N0.y + q[11]*N0.z)
P3 = P2.locatenew('P_3', l23[0]*N2.x + l23[1]*N2.y + l23[2]*N2.z)
P4 = P3.locatenew('P_4', l34[0]*N3.x + l34[1]*N3.y + l34[2]*N3.z)
P45 = P4.locatenew('P45', l45*N3.x)
P5 = P45.locatenew('P_5', l56/2*N5.x)
P6 = P45.locatenew('P', l56*N5.x)

P6.pos_from(P3)

P01 = [P0.locatenew('P01_0', l01[0]*N0.x + l01[1]*N0.y + l01[2]*N0.z),
       P0.locatenew('P01_1', l01[3]*N0.x + l01[4]*N0.y + l01[5]*N0.z),
       P0.locatenew('P01_2', l01[6]*N0.x + l01[7]*N0.y + l01[8]*N0.z),
       P0.locatenew('P01_3', l01[9]*N0.x + l01[10]*N0.y + l01[11]*N0.z)]

P10 = [P1.locatenew('P10_0', l10[0]*N1.x + l10[1]*N1.y + l10[2]*N1.z),
       P1.locatenew('P10_1', l10[3]*N1.x + l10[4]*N1.y + l10[5]*N1.z),
       P1.locatenew('P10_2', l10[6]*N1.x + l10[7]*N1.y + l10[8]*N1.z),
       P1.locatenew('P10_3', l10[9]*N1.x + l10[10]*N1.y + l10[11]*N1.z)]

P12 = [P1.locatenew('P12_0', l12[0]*N1.x + l12[1]*N1.y + l12[2]*N1.z),
       P1.locatenew('P12_1', l12[3]*N1.x + l12[4]*N1.y + l12[5]*N1.z),
       P1.locatenew('P12_2', l12[6]*N1.x + l12[7]*N1.y + l12[8]*N1.z),
       P1.locatenew('P12_3', l12[9]*N1.x + l12[10]*N1.y + l12[11]*N1.z)]

P21 = [P2.locatenew('P21_0', l21[0]*N2.x + l21[1]*N2.y + l21[2]*N2.z),
       P2.locatenew('P21_1', l21[3]*N2.x + l21[4]*N2.y + l21[5]*N2.z),
       P2.locatenew('P21_2', l21[6]*N2.x + l21[7]*N2.y + l21[8]*N2.z),
       P2.locatenew('P21_3', l21[9]*N2.x + l21[10]*N2.y + l21[11]*N2.z)]


N1.set_ang_vel(N0, u[0]*N0.x + u[1]*N0.y + u[2]*N0.z)
N2.set_ang_vel(N0, u[3]*N0.x + u[4]*N0.y + u[5]*N0.z)
N3.set_ang_vel(N2, u[12]*N2.y)

N5.set_ang_vel(N3, u[13]*N3.y)
N6.set_ang_vel(N2, 0)


P1.set_vel(N0, u[6]*N0.x + u[7]*N0.y + u[8]*N0.z)
P2.set_vel(N0,  u[9]*N0.x + u[10]*N0.y + u[11]*N0.z)
P6.set_vel(N2, u[14]*N2.x)
P3.v2pt_theory(P2, N0, N3)
P4.v2pt_theory(P3, N0, N3)
P5.v2pt_theory(P4, N0, N5)
P6.set_vel(N0, P6.vel(N2).express(N0))
P45.v2pt_theory(P2, N0, N3)

cfgconstr = [P6.pos_from(P4).dot(N2.x) - q[14],
             me.dot(N5.x, N2.x) - sm.cos(q[13])]
velconstr = [P6.vel(N2).dot(N6.z),
             P6.vel(N2).dot(N6.y)]

loads = [(P1, -m1 * g * N0.y),
         (P2, -m2 * g * N0.y),
         (P3, -m3 * g * N0.y),
         (P4, -m4 * g * N0.y),
         (P5, -m5 * g * N0.y),
         (P6, -m6 * g * N0.y)]

for i in range(4):
    P01[i].v2pt_theory(P0, N0, N0)
    P10[i].v2pt_theory(P1, N0, N1)
    P12[i].v2pt_theory(P1, N0, N1)
    P21[i].v2pt_theory(P2, N0, N2)
    loads.append((P01[i], -(kg[0] * P01[i].pos_from(P10[i]).dot(N0.x) * N0.x + kg[1] * P01[i].pos_from(P10[i]).dot(N0.y) * N0.y + kg[2] * P01[i].pos_from(P10[i]).dot(N1.z) * N0.z)))
    loads.append((P10[i], -(kg[0] * P10[i].pos_from(P01[i]).dot(N0.x) * N0.x + kg[1] * P10[i].pos_from(P01[i]).dot(N0.y) * N0.y + kg[2] * P10[i].pos_from(P01[i]).dot(N0.z) * N0.z)))
    loads.append((P12[i], -(ks[0] * P12[i].pos_from(P21[i]).dot(N0.x) * N0.x + ks[1] * P12[i].pos_from(P21[i]).dot(N0.y) * N0.y + ks[2] * P12[i].pos_from(P21[i]).dot(N0.z) * N0.z)))
    loads.append((P21[i], -(ks[0] * P21[i].pos_from(P12[i]).dot(N0.x) * N0.x + ks[1] * P21[i].pos_from(P12[i]).dot(N0.y) * N0.y + ks[2] * P21[i].pos_from(P12[i]).dot(N0.z) * N0.z)))
    loads.append((P10[i], -(cg[0] * P10[i].vel(N0).dot(N0.x) * N0.x + cg[1] * P10[i].vel(N0).dot(N0.y) * N0.y + cg[2] * P10[i].vel(N0).dot(N0.z) * N0.z)))
    loads.append((P12[i], -(cs[0] * P12[i].vel(N0).dot(N0.x) * N0.x + cs[1] * P12[i].vel(N0).dot(N0.y) * N0.y + cs[2] * P12[i].vel(N0).dot(N0.z) * N0.z)))
    loads.append((P21[i], -(cs[0] * P21[i].vel(N0).dot(N0.x) * N0.x + cs[1] * P21[i].vel(N0).dot(N0.y) * N0.y + cs[2] * P21[i].vel(N0).dot(N0.z) * N0.z)))

In1 = (me.inertia(N1, I1[0], I1[1], I1[2]), P1)
In2 = (me.inertia(N2, I2[0], I2[1], I2[2]), P2)
In3 = (me.inertia(N3, I3[0], I3[1], I3[2]), P3)
In4 = (me.inertia(N3, I4[0], I4[1], I4[2]), P4)
In5 = (me.inertia(N5, I5[0], I5[1], I5[2]), P5)
In6 = (me.inertia(N6, I6[0], I6[1], I6[2]), P6)

body1 = me.RigidBody('b_1', P1, N1, m1, In1)
body2 = me.RigidBody('b_2', P2, N2, m2, In2)
body3 = me.RigidBody('b_3', P3, N3, m3, In3)
body4 = me.RigidBody('b_4', P4, N3, m4, In4)
body5 = me.RigidBody('b_5', P5, N5, m5, In5)
body6 = me.RigidBody('b_6', P6, N6, m6, In6)

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
        q[11].diff() - u[11],
        q[12].diff() - u[12],
        q[13].diff() - u[13],
        q[14].diff() - u[14]]

KM = me.KanesMethod(N0, q[0:13], u[0:13], kdes,
                    configuration_constraints=cfgconstr,
                    q_dependent=[q[13], q[14]], u_auxiliary=[u[13], u[14]])
fr, frstar = KM.kanes_equations([body1, body2, body3, body4, body5, body6], loads)

consts = {
          I1[0]: 1.0, I1[1]: 1.0, I1[2]: 1.0,
          I2[0]: 1.0, I2[1]: 1.0, I2[2]: 1.0,
          I3[0]: 1.0, I3[1]: 1.0, I3[2]: 1.0,
          I4[0]: 1.0, I4[1]: 1.0, I4[2]: 1.0,
          I5[0]: 1.0, I5[1]: 1.0, I5[2]: 1.0,
          I6[0]: 1.0, I6[1]: 1.0, I6[2]: 1.0,
          kg[0]: 1000.0, kg[1]: 2000.0, kg[2]: 1000.0,
          ks[0]: 200.0, ks[1]: 100.0, ks[2]: 200.0,
          cg[0]: 0.01, cg[1]: 0.01, cg[2]: 0.01,
          cs[0]: 0.01, cs[1]: 0.01, cs[2]: 0.01,
          g: 9.81,
          l01[0]:  0.02,  l01[1]:  0.01,  l01[2]:  0.01,
          l01[3]: -0.01,  l01[4]:  0.01,  l01[5]:  0.01,
          l01[6]: -0.01,  l01[7]:  0.01,  l01[8]: -0.01,
          l01[9]:  0.01, l01[10]:  0.01, l01[11]: -0.01,
          l10[0]:  0.02,  l10[1]: -0.01,  l10[2]:  0.01,
          l10[3]: -0.01,  l10[4]: -0.02,  l10[5]:  0.01,
          l10[6]: -0.01,  l10[7]: -0.01,  l10[8]: -0.01,
          l10[9]:  0.01, l10[10]: -0.01, l10[11]: -0.01,
          l12[0]:  0.01,  l12[1]:  0.01,  l12[2]:  0.01,
          l12[3]: -0.01,  l12[4]:  0.01,  l12[5]:  0.01,
          l12[6]: -0.01,  l12[7]:  0.01,  l12[8]: -0.01,
          l12[9]:  0.01, l12[10]:  0.01, l12[11]: -0.01,
          l21[0]:  0.01,  l21[1]: -0.01,  l21[2]:  0.01,
          l21[3]: -0.01,  l21[4]: -0.01,  l21[5]:  0.01,
          l21[6]: -0.01,  l21[7]: -0.01,  l21[8]: -0.01,
          l21[9]:  0.01, l21[10]: -0.01, l21[11]: -0.01,
          l23[0]:  0.01,  l23[1]:  0.01,  l23[2]:  0.01,
          l34[0]:  0.01,  l34[1]:  0.01,  l34[2]:  0.01,
          l45:  0.01, l56: 0.02,
          m1: 0.01, m2: 0.05, m3: 0.1,
          m4: 0.01, m5: 0.01 #, m6: 0.02
          }

sys = System(KM)
sys.constants = consts
sys.times = np.linspace(0, 100.0, 1000000)
yi1 = consts[l01[1]]-consts[l10[1]]
yi2 = yi1 + consts[l12[1]]-consts[l21[1]]
x6 = consts[l45] + consts[l56]
sys.initial_conditions = {q[0]: 0.0,  q[1]: 0.0,  q[2]: 0.0,
                          q[3]: 0.0,  q[4]: 0.0,  q[5]: 0.0,
                          q[6]: 0.0,  q[7]: yi1,  q[8]: 0.0,
                          q[9]: 0.0, q[10]: yi2, q[11]: 0.0,
                         q[12]: 0.0, q[13]: np.pi, q[14]: x6,
                          u[0]: 0.0,  u[1]: 0.0,  u[2]: 0.0,
                          u[3]: 0.0,  u[4]: 0.0,  u[5]: 0.0,
                          u[6]: 0.0,  u[7]: 0.0,  u[8]: 0.0,
                          u[9]: 0.0, u[10]: 0.0, u[11]: 0.0,
                         u[12]: 1.0#, u[13]: 0.0, u[14]: 0.0
                          }
sys.generate_ode_function(generator='cython')
S = sys.integrate()


fig, axes = plt.subplots(3, 1, sharex=True)
axes[0].plot(sys.times, S[:,6:12])
axes[0].set_ylabel('Displacement [m]')
axes[1].plot(sys.times, S[:,0:6])
axes[1].set_xlabel('Time [s]')
axes[1].set_ylabel('Angle [deg]')
axes[2].plot(sys.times, S[:,12:15])


fig, axes = plt.subplots(3, 1, sharex=True)
axes[0].plot(sys.times, S[:,21:27])
axes[0].set_ylabel('Linear Velocity [m/s]')
axes[1].plot(sys.times, S[:,15:21])
axes[1].set_xlabel('Time [s]')
axes[1].set_ylabel('Angular Rate [deg/s]')
axes[2].plot(sys.times, S[:, 27:28])
