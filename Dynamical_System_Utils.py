from Duffing_Utils import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint

class DynamicalSystem2D():
    """
    Parent Class with predefined collect_data Method.\n
    All Children Classes to be passed to the EDMD class must contain:\n
    - self.mode: either 'Discrete' or 'Continuous', defining the type of Dynamical System in Exam\n
    - collect_data: a function to define two 2-dimensional np.ndarray of shape (2,N) containing the snapshots of the system
    """
    def __init__(self):
        self.mode = "Discrete"
    
    def collect_data(self):
        X = np.random.randn(2,100)
        self.X = X
        self.Y = self.X


class Linear_System(DynamicalSystem2D):
    """
    Children Class of DynamicalSystem2D for Discrete-time Linear Systems.\n
    Takes as input a matrix J defined to be the Transfer Matrix of the Linear System.
    """
    def __init__(self, J:np.ndarray = np.array([[0.9,-.1],[0,0.8]])):
        self.J = J
        self.mode = "Discrete"
    
    def collect_data(self):
        """
        collect_data copies the random initialization of the parent class for matrix X and 
        obtains matrix Y by multiplication by the transfer matrix: Y = J @ X
        """
        super().collect_data()
        self.Y = self.J @ self.X
 


class Duffing_System(DynamicalSystem2D):
    """
    Children Class of DynamicalSystem2D for Linear System.
    Takes as input 3 parameters: 
    - de, defined to be the damping δ
    - be, defined to be the restoring force β
    - al, defined to be the stifness α
    """
    def __init__(self, de:float=0.5,be:float=-1.0,al:float=1.0):
        self.de = de
        self.be = be
        self.al = al
        self.mode = "Continuous"


    def collect_data(self, n:int=10**3, T:float=10.0, timesteps:int=101, mode:str="rand"):
        """
        The function collect_data takes (possibly) additional inputs:
        - n: defining the number of initialization points.
        - T: defining the end-time of our simulation.
        - timesteps: defining the number of timesteps of our simulation.
        - mode: either 'rand' or 'grid', defining the type of initialization, either randomized points or grid-like points.\n
        Note that T/timesteps defines the resolution of the simulation: the smaller, the more accurate the simulation.
        """

        if mode == "rand":
            np.random.seed(0)
            inits = (np.random.rand(n,2)-.5)*6
        elif mode == "grid":
            inits = generate_coord_grid(sqrt_int(n))
        else:
            raise ValueError("Invalid Mode Input, Accepted Values are 'rand' for random inits or 'grid' for grid-space")
        
        t = np.linspace(0, T, timesteps)
        X = np.zeros((2,0))
        Y = np.zeros((2,0))

        for i in range(len(inits)):
            initial_conditions = inits[i]
            solution = odeint(duffing_ode, initial_conditions, t, args=(self.de,self.be,self.al))
            x_data, y_data = solution.T
            X_new = np.vstack((x_data[:100], y_data[:100]))
            Y_new = np.vstack((x_data[1:], y_data[1:]))
            X = np.hstack((X,X_new))
            Y = np.hstack((Y,Y_new))

        self.inits = inits
        self.t = t
        self.X = X
        self.Y = Y
        
    def visualize_system(self, n_traj=100):
        """
        Function for visualization.
        After collecting data from the system with collect_data(), creates plots regarding the simulation:
        - A sample evolution of both coordinates in time
        - The plot of n_traj trajectories in the state space
        """
        try:
            _ = self.X
        except:
            raise Exception("Data Not Generated, Run .collect_data() first!")
        
        fig,axs = plt.subplots(1, 2, figsize=(15,4))
        axs[0].plot(self.X[0,:len(self.t)-1], label="X")
        axs[0].plot(self.X[1,:len(self.t)-1], label="Y")
        axs[0].legend()

        for i in range(n_traj):
            axs[1].plot(self.X[0,i*(len(self.t)-1):(i+1)*(len(self.t)-1)], self.X[1,i*(len(self.t)-1):(i+1)*(len(self.t)-1)])

        axs[0].grid(True)
        axs[1].grid(True)
        plt.suptitle("Sample Trajectories Visualization")
        plt.show()