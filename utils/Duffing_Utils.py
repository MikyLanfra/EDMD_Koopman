import numpy as np

def sqrt_int(n):
    """
    Function to obtain two closest integers multiplying to n
    """
    a = np.floor(np.sqrt(n))
    while bool(n%a):
        a -= 1
    b = n//a
    return a, b


def generate_coord_grid(nx:int, ny:int,xlim=2,ylim=2):
    """
    Function Used to create a numpy.ndarray.\n
    Output has shape [nx,ny] representing the state space in the interval [-xlim, xlim], [-ylim, ylim].
    """
    x = np.linspace(-xlim,xlim,nx)
    y = np.linspace(-ylim,ylim,ny)
    coords = [[x[i],y[j]] for i in range(nx) for j in range(ny)]
    return np.array(coords)

def duffing_ode(state, t, de=0.5,be=-1,al=1):
    """
    Function describing the Duffing ODE governing the system.\n
    Takes as input 3 parameters:\n
    - de, defined to be the damping δ\n
    - be, defined to be the restoring force β\n
    - al, defined to be the stifness α
    """
    x,y = state
    dx =  y
    dy = -de*y-be*x-al*x**3
    return [dx, dy]

