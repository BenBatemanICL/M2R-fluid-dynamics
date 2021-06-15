import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

from scipy.integrate import odeint

def a_sol(y, t, U0, w, nu, H):
    E = np.exp((1+1j)*H*np.sqrt(w/(2*nu)))
    E_inv = np.exp(-(1+1j)*H*np.sqrt(w/(2*nu)))
    def f():
        return ((1-E)*np.exp(-(1+1j)*y*np.sqrt(w/(2*nu))) + (E_inv - 1)*np.exp((1+1j)*y*np.sqrt(w/(2*nu))))/(E_inv - E)
    return U0*np.real(np.exp((w*t)*1j)*f())

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
        new_u.append(nu*(new_v[i]*(u[i+1]-2*u[i]+u[i-1])/h**2 + ((new_v[i+1]-new_v[i-1])/2*h)*(u[i+1]-u[i-1])/2*h))
    new_u.append(-U0*w*np.sin(w*(t-dt)) if t>0 else 0)
    return np.array(new_u)

def old_b_sol(y, t, U0, w, eta, rho, tau, G, H):
    nu = eta/rho
    a1 = tau/(1 + (w**2)*(tau**2))
    a2 = w*(tau**2)/(1+(w**2)*(tau**2))
    b =1j*rho*w/(rho*nu + G*(a1 - 1j*a2))
    E = np.exp(H*np.sqrt(b))
    return U0*np.real(((1-E**-1)*np.exp(w*t*1j + y*np.sqrt(b)) + (E -1)*np.exp(w*t*1j -y*np.sqrt(b)))/(E - E**-1))

pi = np.pi
H=5
N=200
U0=1
w=2
n=1.2
dt=0.01
nu=1.2
t_array = np.arange(0, 8*pi/w + dt, dt)
y_array = np.linspace(0, H, num=N)
init_s = a_sol(y_array, 0, U0, w, nu, H)
newt_sol = np.array([a_sol(y_array, t, U0, w, nu, H) for t in t_array])
pow_sol = odeint(fd3, init_s, t_array, args=(U0, w, nu, N, H, n, dt))

fig, ax = plt.subplots()
ax.set(xlabel=r'$u$', ylabel=r'$y$', title='Velocity field for a newtonian and power law fluid.',
        xlim=(-U0, U0), ylim=(0,H))
newt_line, = ax.plot([], [], label='Newtonian', color='blue')
pow_line, = ax.plot([], [], label='Power law', color='orange')
plt.legend()

def animate(i):
    ax.set(xlim=(-U0, U0), ylim=(0, H))
    newt_line.set_data(newt_sol[i], y_array)
    pow_line.set_data(pow_sol[i], y_array)

an = ani.FuncAnimation(fig, animate, len(t_array), interval=10)
an.save('powvsnewt.mp4')

nu = 1.4
H = 10
N = 100
w = 2
U0=1
mu = rho = 1
tau = 0.5
G = 1
t_array = np.arange(0, 8*pi/w, dt)
y_array = np.linspace(0, H, num=N)
init_s = old_b_sol(y_array, 0, U0, w, mu, rho, tau, G, H)
pow_s = odeint(fd3, init_s, t_array, args=(U0, w, nu, N, H, 1.2, dt))
old = [old_b_sol(y_array, t, U0, w, mu, rho, tau, G, H) for t in t_array]

fig, ax = plt.subplots()
ax.set(xlabel=r'$u$', ylabel=r'$y$', title='Velocity field for an Oldroyd-B and power law fluid.',
        xlim=(-U0, U0), ylim=(0,H))
old_line, = ax.plot([], [], label='Oldroyd-B', color='blue')
pow_line, = ax.plot([], [], label='Power law', color='orange')
plt.legend()

def animate2(i):
    ax.set(xlim=(-U0, U0), ylim=(0, H))
    old_line.set_data(old[i], y_array)
    pow_line.set_data(pow_s[i], y_array)

an = ani.FuncAnimation(fig, animate2, len(t_array), interval=10)
plt.show()
an.save('powvsold.mp4')