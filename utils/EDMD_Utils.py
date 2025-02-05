import numpy.linalg as la
from numpy.polynomial.hermite import hermval as he

from utils.Dynamical_System_Utils import *
from utils.RBF_Utils import *


class EDMD():
    """
    Parent Class of all implementation of EDMD in the code.\n
    Initialization takes as input an object of type DynamicalSystem2D, and collects data if it hasn't been collected yet.\n
    All Children Class derived from this should contain the following methods:
    - psivec(self,x): a function (with possible additional parameters) evaluating the single column of X on the function dictionary
    - get_G_A_K(self): not needed to be changed by children classes, defines the evaluation of matrices G, A and K by the algorithm
    - get_eigendecomp(self): not needed to be changed by children classes, computes the eigendecomposition of K
    """
    def __init__(self, sys:DynamicalSystem2D):
        try:
            _=sys.X
        except:
            sys.collect_data()
        self.X, self.Y = sys.X, sys.Y
        self.mode = sys.mode
        if self.mode == "Continuous":
            self.t = sys.t
            self.dt = self.t[1]-self.t[0]
    
    def psivec(self, x):
        return x

    def get_G_A_K(self):
        """
        Taken the function defined in the class psivec, taking as input a point in the state space and evaluating all functions
        of the dictionary at that point, computes the matrices G,A and K following the EDMD algorithm as follows:
        - G = 1/M sum(psivec(x_i)*psivec(x_i))
        - A = 1/M sum(psivec(x_i)*psivec(y_i))
        - K = pseudoinverse(G) @ A
        """
        M = len(self.X.T)
        G,A = 0, 0
        for i in range(M):
            x1 = self.psivec(self.X[:,i])
            G += np.outer(x1,x1)/M
            A += np.outer(x1, self.psivec(self.Y[:,i]))/M
        K = la.pinv(G) @ A
        self.G = G
        self.A = A
        self.K = K
        return G,A,K
    
    def get_eigendecomp(self):
        """
        Computes the eigendecomposition of matrix K.\n
        If the system is a Continous-time Dynamical System changes the notation of the eigenvalues to the Continuous-Time notation.
        """
        mu, w = la.eig(self.K)
        sort_idx = np.argsort(mu)[::-1]
        mu = mu[sort_idx]
        w = w[:,sort_idx]
        if self.mode == "Continuous":
            mu = np.log(mu)/self.dt
        self.eigvals = mu
        self.w = w
        return mu,w
    


class Hermite_EDMD(EDMD):
    """
    Children Class of EDMD, built with Hermite Polynomials of degree up to d as dictionary.\n
    Hermite Polynomials are particularly suited for EDMD in the case of space R^n and Normally-distributed data.
    """

    def __init__(self, sys: DynamicalSystem2D, degree:int):
        super().__init__(sys)
        self.degree = degree

    def psivec(self, x):
        """
        Creates a function dictionary with all the products of Hermite polynomials of degree from 0 to d in x and polynomials of degree from 0 to n in y (ex. H_0(x)*H_1(y)).
        """
        tuples_list = []
        for i in range(self.degree + 1):
            tuple_i = tuple(1 if j == i else 0 for j in range(i+1))
            tuples_list.append(tuple_i)
        result = np.array([[he(x[0], tuples_list[i]) * he(x[1], tuples_list[j]) for i in range(self.degree+1) for j in range(self.degree+1)]])
        return result
    
    def plot_eigenfunctions(self, nx:int=50, lim=5):
        """
        Plots the evaluation of the eigenfunctions on a square grid [-lim,lim]^2 evaluated at nx datapoints in each direction.
        """
        phi = np.zeros([nx,nx,len(self.w)])
        x = np.linspace(-lim,lim,nx)
        for k in range(len(self.w)):
            for i in range(nx):
                p = self.psivec([x[i],x[:]]).transpose((0,2,1))
                phi[:,i,k] = p @ self.w.T[k]

        fig,axs = plt.subplots(2,4)
        for i in range(2):
            for j in range(4):
                im = axs[i,j].imshow(phi[:,::-1,1+4*i+j],extent=([-5,5,-5,5]))
                axs[i,j].set_title('μ={:.3f}'.format(np.real(self.eigvals[1+4*i+j])))
        plt.tight_layout()



class RBF_EDMD(EDMD):
    """
    Children Class of EDMD, built with RBF functions as dictionary.\n
    It uses KMeans to identify n_centers cluster points from the data (proportionally to the point density in the state space), 
    and defines them to be the interpolation nodes for the RBFs of the algorithm.\n
    It is particularly suited for systems with complex geometries and mesh-free methods.
    """
    def __init__(self, sys: DynamicalSystem2D, n_centers=200):
        super().__init__(sys)
        phix, c = phiMat(self.X, n_centers,return_center=True)
        phix = np.concatenate([np.ones((1,len(phix.T))),phix],axis=0)
        self.phix = phix
        self.centers = c

    def psivec(self,x):
        """
        Creates a Function Dictionary through an RBF interpolation with all the centers obtained during the initialization.\n
        The Kernel used by the RBF is the so-called Thin Plate Spline, that is ϕ(r) = r^2 * log(r)
        """
        r = np.sqrt((x[0]-self.centers[:,0])**2 + (x[1]-self.centers[:,1])**2)
        outlist = np.array([1])
        outlist = np.append(outlist,TPS(r))
        return outlist
    
    def plot_fdict(self):
        """
        Utility Function for plotting the Centers used for interpolation, as well as an example function from the RBF Dictionary.
        """
        fig,axs = plt.subplots(1, 2, figsize=(15,4))
        axs[0].scatter(self.centers[:,0],self.centers[:,1],s=10)
        axs[0].set_title("Cluster Centroids")
        im = axs[1].scatter(self.X[0,:],self.X[1,:], c=self.phix[10,:])
        axs[1].set_title("Plot of random dictionary function")
        plt.colorbar(im,ax=axs)
        axs[0].grid(True)
        axs[1].grid(True)
        plt.show()

    def plot_eigenfunctions(self, n:int=4, lim=2):
        """
        Plots the first n eigenfunctions, evaluated on a square grid [-lim,lim]^2.
        """
        phi_list = []
        for i in range(n):
            phi_list.append(self.phix.T @ self.w[:,i])
        a = n//5 +1
        fig,axs = plt.subplots(a,n//5 + n%5,gridspec_kw={'wspace':.3,'hspace':0},figsize=(18,3))
        for i in range(n):
            im = axs[i].scatter(self.X[0,:],self.X[1,:],c=phi_list[i])
            axs[i].set_xlim([-lim, lim])  # Set x-axis limits
            axs[i].set_ylim([-lim, lim]) 
            axs[i].set_title(f"λ = {np.real(self.eigvals[i]):.5f}+{np.imag(self.eigvals[i]):.2f}i")
        plt.colorbar(im, ax=axs[i])
        plt.show()

    def linearity_check(self, n=4):
        """
        Confronts the evolution of n eigenfunctions approximated by the algorithm and computed on the interpolation nodes, 
        and confronts the obtained result with the theoretical linear evolution φ(t) = φ(0)e^λt.\n
        The returned list contained the approximation errors, computed as the infinity norm between the linear evolution
        and the empirical result of the evolution of the eigenfunction.
        """
        if self.mode != "Continuous":
            return None
        t_max = len(self.t)
        phi_list = np.array([self.phix.T @ self.w[:,i] for i in range(n)])
        phi_th = phi_list[:,:t_max]
        phi_pr = np.array([[phi_th[j,0]*(np.e**(self.eigvals[j]*n*self.dt)) for n in range(t_max)]for j in range(phi_th.shape[0])])
        errs = [max((phi_th[i,:]-phi_pr[i,:])**2) for i in range(phi_th.shape[0])]
        errs = ["{:.5f}".format(np.real(error)) for error in errs]
        print(f"The Plotted Eigenfunctions are approximately linear, with the following Linearity Errors: {errs}")
        # return errs