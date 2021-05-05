import sympy as sm

import sympy.physics.mechanics as me
import autograd as ag
import matplotlib.pyplot as plt
import numpy as np
from pydy.system import System

import pydy


q = me.dynamicsymbols('q:12')
u = me.dynamicsymbols('u0:12')
qd = me.dynamicsymbols('q:12', 1)
ud = me.dynamicsymbols('u0:12', 1)
L = sm.symbols('LA:D0:4x:z')
I = sm.symbols('Ia:bx:z')
g = sm.symbols('g')
ma, mb = sm.symbols('ma:b')
k = sm.symbols('k0:2x:z')
Lgr, Lsp = sm.symbols('Lgr Lsp')


N = me.ReferenceFrame('N')
Orig = me.Point('O')
Orig.set_vel(N, 0*N.x + 0*N.y + 0*N.z)
OA = [Orig.locatenew('OA_0', L[0]*N.x + L[1]*N.y + L[2]*N.z),
      Orig.locatenew('OA_1',-L[3]*N.x + L[4]*N.y + L[5]*N.z),
      Orig.locatenew('OA_2',-L[6]*N.x + L[7]*N.y - L[8]*N.z),
      Orig.locatenew('OA_3', L[9]*N.x + L[10]*N.y - L[11]*N.z)]


A = N.orientnew('A', 'Body', (q[0], q[1], q[2]), '123')
A.set_ang_vel(N, u[0] * N.x + u[1] * N.y + u[2] * N.z)
CMa = Orig.locatenew('CM_s', q[3]*N.x + q[4]*N.y + q[5]*N.z)
CMa.set_vel(N, u[3]*N.x + u[4]*N.y + u[5]*N.z)
AO = [CMa.locatenew('AO_0', L[12]*A.x - L[13]*A.y + L[14]*A.z),
      CMa.locatenew('AO_1',-L[15]*A.x - L[16]*A.y + L[17]*A.z),
      CMa.locatenew('AO_2',-L[18]*A.x - L[19]*A.y - L[20]*A.z),
      CMa.locatenew('AO_3', L[21]*A.x - L[22]*A.y - L[23]*A.z)]
AB = [CMa.locatenew('AB_0', L[24]*A.x + L[25]*A.y + L[26]*A.z),
      CMa.locatenew('AB_1',-L[27]*A.x + L[28]*A.y + L[29]*A.z),
      CMa.locatenew('AB_2',-L[30]*A.x + L[31]*A.y - L[32]*A.z),
      CMa.locatenew('AB_3', L[33]*A.x + L[34]*A.y - L[35]*A.z)]


B = A.orientnew('B', 'Body', (q[6], q[7], q[8]), '123')
B.set_ang_vel(A, u[6]*A.x + u[7]*A.y + u[8]*A.z)
CMb = CMa.locatenew('CM_b', q[9]*A.x + q[10]*A.y + q[11]*A.z)
CMb.set_vel(A, u[9]*A.x + u[10]*A.y + u[11]*A.z)
CMb.set_vel(N, CMb.vel(A).express(N))
BA = [CMb.locatenew('BA_0', L[36]*B.x - L[37]*B.y + L[38]*B.z),
      CMb.locatenew('BA_1',-L[39]*B.x - L[40]*B.y + L[41]*B.z),
      CMb.locatenew('BA_2',-L[42]*B.x - L[43]*B.y - L[44]*B.z),
      CMb.locatenew('BA_3', L[45]*B.x - L[46]*B.y - L[47]*B.z)]

In = [me.inertia(A, I[0], I[1], I[2]),
      me.inertia(B, I[3], I[4], I[5])]
bodies = [me.RigidBody('Shell', CMa, A, ma, (In[0], CMa)),
          me.RigidBody('Block', CMb, B, mb, (In[1], CMb))]

loads = [(CMa, -ma*g*N.y),
         (CMb, -mb*g*N.y)]
# CMa.set_vel(N, CMa.pos_from(Orig).dt(N))
# CMb.set_vel(N, CMb.pos_from(Orig).dt(N))
# CMb.set_vel(A, CMb.pos_from(CMa).dt(A))


for _i in range(4):
    OA[_i].v2pt_theory(Orig, N, N)
    AO[_i].v2pt_theory(CMa, N, A)
    AB[_i].v2pt_theory(CMa, N, A)
    BA[_i].v2pt_theory(CMb, N, B)
    BA[_i].v2pt_theory(CMb, A, B)

    loads.append((OA[_i], -k[0]*OA[_i].pos_from(AO[_i]).dot(N.x)*N.x - k[1]*(Lgr - OA[_i].pos_from(AO[_i]).dot(N.y))*N.y - k[2]*OA[_i].pos_from(AO[_i]).dot(N.z)*N.z))
    loads.append((AO[_i], -k[0]*AO[_i].pos_from(OA[_i]).dot(N.x)*N.x - k[1]*(Lgr - AO[_i].pos_from(OA[_i]).dot(N.y))*N.y - k[2]*AO[_i].pos_from(OA[_i]).dot(N.z)*N.z))
    loads.append((AB[_i], -k[3]*AB[_i].pos_from(BA[_i]).dot(A.x)*A.x - k[4]*(Lsp - AB[_i].pos_from(BA[_i]).dot(A.y))*A.y - k[5]*AB[_i].pos_from(BA[_i]).dot(A.z)*A.z))
    loads.append((BA[_i], -k[3]*BA[_i].pos_from(AB[_i]).dot(A.x)*A.x - k[4]*(Lsp - BA[_i].pos_from(AB[_i]).dot(A.y))*A.y - k[5]*BA[_i].pos_from(AB[_i]).dot(A.z)*A.z))

# eqs = [sm.Eq(u[0], A.ang_vel_in(N).dot(N.x)),
#        sm.Eq(u[1], A.ang_vel_in(N).dot(N.y)),
#        sm.Eq(u[2], A.ang_vel_in(N).dot(N.z)),
#        sm.Eq(u[3], CMa.vel(N).dot(N.x)),
#        sm.Eq(u[4], CMa.vel(N).dot(N.y)),
#        sm.Eq(u[5], CMa.vel(N).dot(N.z)),
#        sm.Eq(u[6], B.ang_vel_in(A).dot(A.x)),
#        sm.Eq(u[7], B.ang_vel_in(A).dot(A.y)),
#        sm.Eq(u[8], B.ang_vel_in(A).dot(A.z)),
#        sm.Eq(u[9],  CMb.vel(A).dot(A.x)),
#        sm.Eq(u[10], CMb.vel(A).dot(A.y)),
#        sm.Eq(u[11], CMb.vel(A).dot(A.z))]

# qdots = [sm.solve(eqs,
#                   q[0].diff(),
#                   q[1].diff(),
#                   q[2].diff(),
#                   q[3].diff(),
#                   q[4].diff(),
#                   q[5].diff(),
#                   q[6].diff(),
#                   q[7].diff(),
#                   q[8].diff(),
#                   q[9].diff(),
#                   q[10].diff(),
#                   q[11].diff())]


# for _i in range(4):
#     AO[_i].set_vel(A, AO[_i].vel(N).express(A))
#     AB[_i].set_vel(A, AB[_i].vel(N).express(A))

# A.set_ang_vel(N, A.ang_vel_in(N).subs(qdots).simplify())
# B.set_ang_vel(A, B.ang_vel_in(A).subs(qdots).simplify())
# CMa.set_vel(N, CMa.vel(N).subs(qdots).simplify())
# CMb.set_vel(A, CMb.vel(A).subs(qdots).simplify())
import sympy as sm
import sympy.physics.mechanics as me
import matplotlib.pyplot as plt
import numpy as np
from pydy.system import System
from scipy.integrate import odeint, solve_ivp
import pydy


q = me.dynamicsymbols('q:12')
u = me.dynamicsymbols('u0:12')
L = sm.symbols('LA:D0:4x:z')
I = sm.symbols('Ia:bx:z')
g = sm.symbols('g')
ma, mb = sm.symbols('ma:b')
k = sm.symbols('k0:2x:z')
Lgr, Lsp = sm.symbols('Lgr Lsp')


N = me.ReferenceFrame('N')
Orig = me.Point('O')
Orig.set_vel(N, 0*N.x + 0*N.y + 0*N.z)
OA = [Orig.locatenew('OA_0', L[0]*N.x + L[1]*N.y + L[2]*N.z),
      Orig.locatenew('OA_1',-L[3]*N.x + L[4]*N.y + L[5]*N.z),
      Orig.locatenew('OA_2',-L[6]*N.x + L[7]*N.y - L[8]*N.z),
      Orig.locatenew('OA_3', L[9]*N.x + L[10]*N.y - L[11]*N.z)]


A = N.orientnew('A', 'Body', (q[0], q[1], q[2]), '123')
A.set_ang_vel(N, u[0] * N.x + u[1] * N.y + u[2] * N.z)
CMa = Orig.locatenew('CM_s', q[3]*N.x + q[4]*N.y + q[5]*N.z)
CMa.set_vel(N, u[3]*N.x + u[4]*N.y + u[5]*N.z)
AO = [CMa.locatenew('AO_0', L[12]*A.x - L[13]*A.y + L[14]*A.z),
      CMa.locatenew('AO_1',-L[15]*A.x - L[16]*A.y + L[17]*A.z),
      CMa.locatenew('AO_2',-L[18]*A.x - L[19]*A.y - L[20]*A.z),
      CMa.locatenew('AO_3', L[21]*A.x - L[22]*A.y - L[23]*A.z)]
AB = [CMa.locatenew('AB_0', L[24]*A.x + L[25]*A.y + L[26]*A.z),
      CMa.locatenew('AB_1',-L[27]*A.x + L[28]*A.y + L[29]*A.z),
      CMa.locatenew('AB_2',-L[30]*A.x + L[31]*A.y - L[32]*A.z),
      CMa.locatenew('AB_3', L[33]*A.x + L[34]*A.y - L[35]*A.z)]


B = A.orientnew('B', 'Body', (q[6], q[7], q[8]), '123')
B.set_ang_vel(A, u[6]*A.x + u[7]*A.y + u[8]*A.z)
CMb = CMa.locatenew('CM_b', q[9]*A.x + q[10]*A.y + q[11]*A.z)
CMb.set_vel(A, u[9]*A.x + u[10]*A.y + u[11]*A.z)
CMb.set_vel(N, CMb.vel(A).express(N))
BA = [CMb.locatenew('BA_0', L[36]*B.x - L[37]*B.y + L[38]*B.z),
      CMb.locatenew('BA_1',-L[39]*B.x - L[40]*B.y + L[41]*B.z),
      CMb.locatenew('BA_2',-L[42]*B.x - L[43]*B.y - L[44]*B.z),
      CMb.locatenew('BA_3', L[45]*B.x - L[46]*B.y - L[47]*B.z)]

In = [me.inertia(A, I[0], I[1], I[2]),
      me.inertia(B, I[3], I[4], I[5])]
bodies = [me.RigidBody('Shell', CMa, A, ma, (In[0], CMa)),
          me.RigidBody('Block', CMb, B, mb, (In[1], CMb))]

loads = [(CMa, -ma*g*N.y),
         (CMb, -mb*g*N.y)]
# CMa.set_vel(N, CMa.pos_from(Orig).dt(N))
# CMb.set_vel(N, CMb.pos_from(Orig).dt(N))
# CMb.set_vel(A, CMb.pos_from(CMa).dt(A))


for _i in range(4):
    OA[_i].v2pt_theory(Orig, N, N)
    AO[_i].v2pt_theory(CMa, N, A)
    AB[_i].v2pt_theory(CMa, N, A)
    BA[_i].v2pt_theory(CMb, N, B)
    BA[_i].v2pt_theory(CMb, A, B)

    loads.append((OA[_i], -k[0]*OA[_i].pos_from(AO[_i]).dot(N.x)*N.x - k[1]*(Lgr - OA[_i].pos_from(AO[_i]).dot(N.y))*N.y - k[2]*OA[_i].pos_from(AO[_i]).dot(N.z)*N.z))
    loads.append((AO[_i], -k[0]*AO[_i].pos_from(OA[_i]).dot(N.x)*N.x - k[1]*(Lgr - AO[_i].pos_from(OA[_i]).dot(N.y))*N.y - k[2]*AO[_i].pos_from(OA[_i]).dot(N.z)*N.z))
    loads.append((AB[_i], -k[3]*AB[_i].pos_from(BA[_i]).dot(A.x)*A.x - k[4]*(Lsp - AB[_i].pos_from(BA[_i]).dot(A.y))*A.y - k[5]*AB[_i].pos_from(BA[_i]).dot(A.z)*A.z))
    loads.append((BA[_i], -k[3]*BA[_i].pos_from(AB[_i]).dot(A.x)*A.x - k[4]*(Lsp - BA[_i].pos_from(AB[_i]).dot(A.y))*A.y - k[5]*BA[_i].pos_from(AB[_i]).dot(A.z)*A.z))

eqs = [sm.Eq(qd[0], A.ang_vel_in(N).dot(N.x)),
       sm.Eq(qd[1], A.ang_vel_in(N).dot(N.y)),
       sm.Eq(qd[2], A.ang_vel_in(N).dot(N.z)),
       sm.Eq(qd[3], CMa.vel(N).dot(N.x)),
       sm.Eq(qd[4], CMa.vel(N).dot(N.y)),
       sm.Eq(qd[5], CMa.vel(N).dot(N.z)),
       sm.Eq(qd[6], B.ang_vel_in(A).dot(A.x)),
       sm.Eq(qd[7], B.ang_vel_in(A).dot(A.y)),
       sm.Eq(qd[8], B.ang_vel_in(A).dot(A.z)),
       sm.Eq(qd[9],  CMb.vel(A).dot(A.x)),
       sm.Eq(qd[10], CMb.vel(A).dot(A.y)),
       sm.Eq(qd[11], CMb.vel(A).dot(A.z))]

qdots = [sm.solve(eqs,
                  q[0].diff(),
                  q[1].diff(),
                  q[2].diff(),
                  q[3].diff(),
                  q[4].diff(),
                  q[5].diff(),
                  q[6].diff(),
                  q[7].diff(),
                  q[8].diff(),
                  q[9].diff(),
                  q[10].diff(),
                  q[11].diff())]


# for _i in range(4):
#     AO[_i].set_vel(A, AO[_i].vel(N).express(A))
#     AB[_i].set_vel(A, AB[_i].vel(N).express(A))

# A.set_ang_vel(N, A.ang_vel_in(N).subs(qdots).simplify())
# B.set_ang_vel(A, B.ang_vel_in(A).subs(qdots).simplify())
# CMa.set_vel(N, CMa.vel(N).subs(qdots).simplify())
# CMb.set_vel(A, CMb.vel(A).subs(qdots).simplify())

# for _i in range(4):
#    OA[_i].set_vel(N, OA[_i].vel(N).subs(qdots).simplify())
#    AO[_i].set_vel(N, AO[_i].vel(N).subs(qdots).simplify())
#    AB[_i].set_vel(N, AB[_i].vel(N).subs(qdots).simplify())
#    BA[_i].set_vel(N, BA[_i].vel(N).subs(qdots).simplify())


# kdes = [eqs[0].rhs - eqs[0].lhs,
#         eqs[1].rhs - eqs[1].lhs,
#         eqs[2].rhs - eqs[2].lhs,
#         eqs[3].rhs - eqs[3].lhs,
#         eqs[4].rhs - eqs[4].lhs,
#         eqs[5].rhs - eqs[5].lhs,
#         eqs[6].rhs - eqs[6].lhs,
#         eqs[7].rhs - eqs[7].lhs,
#         eqs[8].rhs - eqs[8].lhs,
#         eqs[9].rhs - eqs[9].lhs,
#         eqs[10].rhs - eqs[10].lhs,
#         eqs[11].rhs - eqs[11].lhs]

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


# for _i in range(4):
#    OA[_i].set_vel(N, OA[_i].vel(N).subs(qdots).simplify())
#    AO[_i].set_vel(N, AO[_i].vel(N).subs(qdots).simplify())
#    AB[_i].set_vel(N, AB[_i].vel(N).subs(qdots).simplify())
#    BA[_i].set_vel(N, BA[_i].vel(N).subs(qdots).simplify())


rhs =  [eqs[0].rhs,
        eqs[1].rhs,
        eqs[2].rhs,
        eqs[3].rhs,
        eqs[4].rhs,
        eqs[5].rhs,
        eqs[6].rhs,
        eqs[7].rhs,
        eqs[8].rhs,
        eqs[9].rhs,
        eqs[10].rhs,
        eqs[11].rhs]

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

kane = me.KanesMethod(N, q, u, kd_eqs=kdes)

Fr, Frst = kane.kanes_equations(bodies, loads=loads)

sys = System(kane)

sys.constants = {I[0]: 1.0, I[1]: 1.0, I[2]: 1.0,
                 I[3]: 1.0, I[4]: 1.0, I[5]: 1.0,

                 L[0]: 0.3, L[1]: 0.0, L[2]: 0.5,
                 L[3]: 0.3, L[4]: 0.0, L[5]: 0.5,
                 L[6]: 0.3, L[7]: 0.0, L[8]: 0.6,
                 L[9]:  0.3, L[10]: 0.0, L[11]: 0.6,

                 L[12]: 0.3, L[13]: 0.02, L[14]: 0.5,
                 L[15]: 0.3, L[16]: 0.02, L[17]: 0.5,
                 L[18]: 0.3, L[19]: 0.02, L[20]: 0.6,
                 L[21]: 0.3, L[22]: 0.02, L[23]: 0.6,

                 L[24]: 0.01, L[25]: 0.03, L[26]: 0.02,
                 L[27]: 0.01, L[28]: 0.03, L[29]: 0.02,
                 L[30]: 0.01, L[31]: 0.03, L[32]: 0.02,
                 L[33]: 0.01, L[34]: 0.03, L[35]: 0.02,

                 L[36]: 0.01, L[37]: 0.04, L[38]: 0.02,
                 L[39]: 0.01, L[40]: 0.04, L[41]: 0.02,
                 L[42]: 0.01, L[43]: 0.04, L[44]: 0.02,
                 L[45]: 0.01, L[46]: 0.04, L[47]: 0.02,
                 Lgr: 0.01,
                 Lsp: 0.04,
                 g: 9.81,
                 k[0]: 100.0, k[1]: 100.0, k[2]: 100.0,
                 k[3]: 50.0, k[4]: 10.0, k[5]: 50.0,
                 ma: 0.2,
                 mb: 1.0}

sys.initial_conditions = {q[0]: 0.0,
                          q[1]: 0.0,
                          q[2]: 0.0,
                          q[3]: 0.0,
                          q[4]: 0.03,
                          q[5]: 0.0,
                          q[6]: 0.0,
                          q[7]: 0.0,
                          q[8]: 0.0,
                          q[9]: 0.0,
                          q[10]: 0.14,
                          q[11]: 0.0,
                          u[0]: 0.0,
                          u[1]: 0.0,
                          u[2]: 0.0,
                          u[3]: 0.0,
                          u[4]: 0.0,
                          u[5]: 0.0,
                          u[6]: 0.0,
                          u[7]: 0.0,
                          u[8]: 0.0,
                          u[9]: 0.0,
                          u[10]: 0.0,
                          u[11]: 0.0}

sys.times = np.linspace(0, 10.0, 1000000)


sys.generate_ode_function(generator='cython')

x = sys.integrate()
plt.plot(sys.times, x)


plt.plot(sys.times, x[:, 0:3])
plt.plot(sys.times, x[:, 3:6])
plt.plot(sys.times, x[:, 6:9])
plt.plot(sys.times, x[:, 9:12])
