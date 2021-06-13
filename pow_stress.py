import numpy as np
import matplotlib.pyplot as plt
import matplotlib
if True:
    matplotlib.use("pgf")
    matplotlib.rcParams['axes.unicode_minus'] = False
    matplotlib.rcParams.update({
        "pgf.texsystem": "pdflatex",
        'font.family': 'serif',
        'text.usetex': True,
        'pgf.rcfonts': False,
    })

def sigma(u, h, mu, n):
    N = len(u)
    dudy = np.array([0]+ [(u[i+1] - u[i-1])/(2*h) for i in range(1, N-1)] + [0])
    return mu*abs(dudy)**(n-1)*dudy

