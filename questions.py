import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from scipy.integrate import odeint
if 1:
    matplotlib.use("pgf")
    matplotlib.rcParams['axes.unicode_minus'] = False
    matplotlib.rcParams.update({
        "pgf.texsystem": "pdflatex",
        'font.family': 'serif',
        'text.usetex': True,
        'pgf.rcfonts': False,
    })

def fd2(u, t, U0, w, nu, N, H):
    h = H/N
    new_u = [-U0*w*np.sin(w*t)]
    for i in range(1, N-1):
        new_u.append(nu*(u[i+1]-2*u[i]+u[i-1])/h**2)
    new_u.append(-U0*w*np.sin(w*t))
    return np.array(new_u)

def a_sol(y, t, U0, w, nu, H):
    E = np.exp((1+1j)*H*np.sqrt(w/(2*nu)))
    E_inv = np.exp(-(1+1j)*H*np.sqrt(w/(2*nu)))
    def f():
        return ((1-E)*np.exp(-(1+1j)*y*np.sqrt(w/(2*nu))) + (E_inv - 1)*np.exp((1+1j)*y*np.sqrt(w/(2*nu))))/(E_inv - E)
    return U0*np.real(np.exp((w*t)*1j)*f())

H=1
U0 = 1
w = 1
#nu = 1
t_array = np.linspace(0, 5)


#N_array = np.linspace(100, 1000, num=10)
#error_array = np.empty((4,10))
#nu_array = np.array([1.,1.5,2.,2.5])
#for i in range(4):
#    nu = nu_array[i]
#    t_error_array = []
#    for N2 in N_array:
#        y_array = np.linspace(0, H, num=int(N2))
#        init_sol = a_sol(y_array, 0, U0, w, nu, H)
#        num_sol = odeint(fd2, init_sol, t_array, args=(U0, w, nu, int(N2), H))
#        ac_sol = np.array([a_sol(y_array, t, U0, w, nu, H) for t in t_array])
#        t_error_array.append(np.mean(abs(num_sol - ac_sol)))
#    error_array[i,:] = np.array(t_error_array)
#fig = plt.figure()
#ax = fig.add_subplot()
#ax.set_ylabel("Residual")
#ax.set_xlabel("N")
#ax.set_title("Residual of odeint vs the analytical solution")
#ax.plot(N_array, error_array[0], label='nu = 1')
#ax.plot(N_array, error_array[1], label='nu = 1.5')
#ax.plot(N_array, error_array[2], label='nu = 2')
#ax.plot(N_array, error_array[3], label='nu = 2.5')
#plt.legend()
#plt.savefig('residuals.pgf')



def fd3(u, t, U0, w, nu, N, H, n, dt):
    h = H/N
    def v():
        new_v = []
        for i in range(N-1):
            new_v.append(abs((u[i+1]-u[i-1])/2*h)**(n-1))
        new_v.append(new_v[0])
        return np.array(new_v)
    new_v = v()
    new_u = [-U0*w*np.sin(w*(t-dt)) if t>0 else 0]
    for i in range(1, N-1):
        new_u.append(nu*(new_v[i]*(u[i+1]-2*u[i]+u[i-1])/h**2 + \
             ((new_v[i+1]-new_v[i-1])/2*h)*(u[i+1]-u[i-1])/2*h))
    new_u.append(-U0*w*np.sin(w*(t-dt)) if t>0 else 0)
    return np.array(new_u)


nu = 1
H = 10
N = 200
w = 2
t_array = np.linspace(0, 3, num=300)
dt = t_array[1] - t_array[0]
y_array = np.linspace(0, H, num=N)
init_s = a_sol(y_array, 0, U0, w, nu, H)
#pow_sol = odeint(fd3, init_s, t_array, args=(U0, w, nu, N, H, 1.5, dt)) #works fine with n=1. Any other value of n, though...
#print(output)
#newt_sol = np.array([a_sol(y_array, t, U0, w, nu, H) for t in t_array])

#fig, (newt_ax, pow_ax) = plt.subplots(1, 2)
#newt_ax.set_xlabel('u, velocity of fluid')
#newt_ax.set_ylabel('y')
#newt_ax.plot(newt_sol[0], y_array, label='t=0')
#newt_ax.plot(newt_sol[99], y_array, label='t=1')
#newt_ax.plot(newt_sol[199], y_array, label='t=2')
#newt_ax.plot(newt_sol[299], y_array, label='t=3')

#pow_ax.set_xlabel('u, velocity of fluid')
#pow_ax.set_ylabel('y')
#pow_ax.plot(pow_sol[0], y_array, label='t=0')
#pow_ax.plot(pow_sol[99], y_array, label='t=1')
#pow_ax.plot(pow_sol[199], y_array, label='t=2')
#pow_ax.plot(pow_sol[299], y_array, label='t=3')
#plt.show()

#------------now see the effect of changing n---------------------
N=400
H=5
w=2
nu=2
y_array = np.linspace(0, H, num=N)
n_array = np.array([0.8, 0.9, 1])
t_array = np.arange(0, 1.01, 0.01)
init_s = a_sol(y_array, 0, U0, w, nu, H)
ns = np.array([odeint(fd3, init_s, t_array, args=(U0, w, nu, N, H, n, dt)) for n in n_array])

fig, (n1, n2, n3) = plt.subplots(1, 3)
n1.set(xlabel='u', ylabel='y')
n2.set(xlabel='u', ylabel='y')
n3.set(xlabel='u', ylabel='y')

n1.plot(ns[0, 0], y_array, label="t=0")
n1.plot(ns[0, 19], y_array, label=f"t={t_array[19]}")
n1.plot(ns[0, 39], y_array, label=f"t={t_array[39]}")
n1.plot(ns[0, 59], y_array, label=f"t={t_array[59]}")
n1.plot(ns[0, 79], y_array, label=f"t={t_array[79]}")
n1.plot(ns[0, 99], y_array, label=f"t={t_array[99]}")

n2.plot(ns[1, 0], y_array, label="t=0")
n2.plot(ns[1, 19], y_array, label=f"t={t_array[19]}")
n2.plot(ns[1, 39], y_array, label=f"t={t_array[39]}")
n2.plot(ns[1, 59], y_array, label=f"t={t_array[59]}")
n2.plot(ns[1, 79], y_array, label=f"t={t_array[79]}")
n2.plot(ns[1, 99], y_array, label=f"t={t_array[99]}")

n3.plot(ns[2, 0], y_array, label="t=0")
n3.plot(ns[2, 19], y_array, label=f"t={t_array[19]}")
n3.plot(ns[2, 39], y_array, label=f"t={t_array[39]}")
n3.plot(ns[2, 59], y_array, label=f"t={t_array[59]}")
n3.plot(ns[2, 79], y_array, label=f"t={t_array[79]}")
n3.plot(ns[2, 99], y_array, label=f"t={t_array[99]}")

plt.legend()
plt.show()
#plt.savefig('different_n_lessthan1.pgf')

def sigma(u, h, mu, n):
    N = len(u)
    dudy = np.array([(u[1]-u[0])/h]+ [(u[i+1] - u[i-1])/(2*h) for i in range(1, N-1)] + [(u[-1]-u[-2])/h])
    return mu*pow(abs(dudy),n)
H=5
w=2
U0=1
N=200
nu=3
n = 1.3
dt=0.1
t_array = np.arange(0, 6+dt, dt)
y_array = np.linspace(0, H, num=N)
init_s = a_sol(y_array, 0, U0, w, nu, H)
ns = odeint(fd3, init_s, t_array, args=(U0, w, nu, N, H, n, dt))
h = y_array[1] - y_array[0]
newt_sol = np.array([a_sol(y_array, t, U0, w, nu, H) for t in t_array])
sigma_t = np.array([sigma(ns[i], h, 1, n) for i in range(len(t_array))])
sigma_t_newt = np.array([sigma(newt_sol[i], h, 1, 1) for i in range(len(t_array))])

#fig, ax = plt.subplots(1,1)
#ax.set(xlabel=r'$\sigma_{12}$', ylabel='y', title='The stress on a power law fluid at various times')
#for ti in [0, 99, 199, 299]:
#    ax.plot(sigma_t[ti], y_array, label=f't={t_array[ti]}')
#plt.legend()
#plt.show()

T, Y = np.meshgrid(t_array, y_array)
fig, (newt_ax, pow_ax) = plt.subplots(1, 2)
pow_ax.set(xlabel='t', ylabel='y', title='Power law', adjustable='box', aspect='equal')
im = pow_ax.pcolor(T, Y, abs(sigma_t.T), shading='auto', cmap='hot', rasterized=True)
cbar = fig.colorbar(im, ax=pow_ax, fraction=0.046, pad=0.04)
cbar.set_label(r'$\sigma_{12}$')

newt_ax.set(xlabel='t', ylabel='y', title='Newtonian', adjustable='box', aspect='equal')
im = newt_ax.pcolor(T, Y, abs(sigma_t_newt.T), shading='auto', cmap='hot', rasterized=True)
cbar = fig.colorbar(im, ax=newt_ax, fraction=0.046, pad=0.04)
cbar.set_label(r'$\sigma_{12}$')
fig.tight_layout()
plt.savefig('heatmap_sigma_powervsnewt.pgf')

def eff(u, h, mu, n):
    N = len(u)
    dudy = np.array([(u[1]-u[0])/h]+ [(u[i+1] - u[i-1])/(2*h) for i in range(1, N-1)] + [(u[-1]-u[-2])/h])
    return mu*pow(abs(dudy), n-1)

mu_eff = np.array([eff(ns[i], h, 1, n) for i in range(len(t_array))])
fig, pow_ax = plt.subplots()
pow_ax.set(xlabel='t', ylabel='y', title='Effective viscosity of power law fluid')
im = pow_ax.pcolor(T, Y, mu_eff.T, shading='auto', cmap='hot', rasterized=True)
cbar = fig.colorbar(im, ax=pow_ax)
cbar.set_label(r'$\eta_{eff}$')


plt.savefig('heatmap_nu_power.pgf')
#plt.show()