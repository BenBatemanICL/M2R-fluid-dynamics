import numpy as np
import matplotlib.pyplot as plt
import matplotlib
pgf = True
if pgf:
    matplotlib.use("pgf")
    matplotlib.rcParams['axes.unicode_minus'] = False
    matplotlib.rcParams.update({
        "pgf.texsystem": "pdflatex",
        'font.family': 'serif',
        'text.usetex': True,
        'pgf.rcfonts': False,
    })

def a_sol(y, t, U0, w, nu, H):
    E = np.exp((1+1j)*H*np.sqrt(w/(2*nu)))
    def f():
        return (U0/2)*((1-E)*np.exp(-(1+1j)*y*np.sqrt(w/(2*nu))) \
            + (E**-1 - 1)*np.exp((1+1j)*y*np.sqrt(w/(2*nu))))/(E**-1 - E)
    return 2*np.real(np.exp((w*t)*1j)*f())

t_array = np.arange(1,6,1)
y_array = np.linspace(0, 10)
sol = [a_sol(y_array, t, 1, 1, 1.2, 10) for t in t_array]

fig, ax = plt.subplots()
ax.set(xlabel=r'$u$', ylabel=r'$y$', title=r'Velocity field at time $t$')
for t in [0,1,2,3,4]:
    ax.plot(sol[t], y_array, label=f'time={t}')
plt.legend()
plt.savefig('newt_visualisation.pgf')
