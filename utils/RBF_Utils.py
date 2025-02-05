import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import euclidean_distances
from scipy.interpolate import RBFInterpolator

def centers(data, n_clusters:int):
    """
    Function Used to determine a number of clusters from an array-like of data.\n
    Returns a list containing the cluster centers.
    """
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(data)
    return kmeans.cluster_centers_

def TPS(r):
    """
    Function used as Kernel the evaluation of the Radial Basis Function.\n
    The choice r^2*log(r) corresponds to the so-called Thin Plate Spline (TPS).
    """
    return np.power(r, 2)*np.log(r+10e-4)

def phiMat(X, n_clusters=100, return_center=False):
    """
    Function to evaluate the RBF associated to the data given as input.\n
    return_center is a parameter defining what to return:\n
    - False: only the array ϕ of RBFs computed on the data\n
    - True: both the array ϕ of RBFs and the cluster centers list
    """
    c = centers(X.T, n_clusters)
    eudist = euclidean_distances(X.T, c)
    phi =  TPS(eudist)
    return (phi.T,c) if return_center else phi.T