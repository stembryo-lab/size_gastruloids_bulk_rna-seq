import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
import statsmodels.api as sm


def qc_threshold(data: np.ndarray = None, n_sigma: int = 3, min_rel: int = 1):
    """
    Calculates thresholds for quality control based on a Kernel Density Estimate (KDE) of the input data.
    
    Parameters:
    -----------
    data : np.ndarray
        The input data, a 1D or 2D array (usually values for QC analysis).
    n_sigma : int, optional
        The number of standard deviations to define the upper threshold (default is 3).
    min_rel : int, optional
        Minimum relative difference threshold for bimodal detection (default is 1).
    
    Returns:
    --------
    min_ret : float
        The lower threshold value.
    max_ret : float
        The upper threshold value.
    """

    # Reshape the data to a 2D array with one column (necessary for fitting)
    X=data.reshape(-1, 1)
    m = np.asarray(sm.robust.scale.huber(X))
    
    # Compute the mean and standard deviation of the data
    mean = m[0]
    std = m[1]
    # loop over bandwith until just one (or none) minima is found below mean
    for bandwidth in np.linspace(std/10, std, 100):
        # Fit a Kernel Density Estimate (KDE) using a Gaussian kernel and specified bandwidth
        kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(X)

        # Generate a range of values from the minimum to maximum of the data for KDE evaluation
        x = np.linspace(X.min(), X.max(), 1000).reshape(-1, 1)

        # Score the sample (log-density estimation) at points in x
        y = kde.score_samples(x)

        # Find local minima in the KDE curve (points where the curve goes down)
        mins = argrelextrema(y, np.less)
        local_min = y[mins]

        # If there's exactly one minimum below the mean
        if sum(local_min < mean) == 1:
            # Use that minimum as the threshold
            x_min = x[mins[0][0]]
            idx_min = mins[0][0]
            break
        # If no minima are above the mean, use the first value in x as the minimum
        elif sum(local_min > mean) == 0:
            x_min = x[0]
            idx_min = 0
            break
    # Find local maxima in the KDE curve (points where the curve goes up)
    maxs = argrelextrema(y, np.greater)
    local_max = x[maxs]

    # If only one local maximum is found
    if len(local_max) == 1:
        # Set x_max1 to the first value in x, and x_max2 to the found local maximum
        x_max1 = x[0]
        x_max2 = local_max[0]
        idy_max2 = maxs[0][0]
    else:
        # Otherwise, use the first two maxima
        x_max1 = local_max[0]
        x_max2 = local_max[1]
        idy_max2 = maxs[0][1]
    # If the distribution is bimodal and the relative difference at x_min is less than min_rel
    if x_min > x_max1 and x_min < x_max2 and (np.exp(y[idx_min] - y[idy_max2]) < min_rel):
        # Set new minimum based on data above x_min
        m = np.asarray(sm.robust.scale.huber([X > x_min]))
        min_ret = x_min

    else:
        min_ret = min(data)

    # Define the upper threshold based on the mean and n_sigma times the standard deviation
    m = np.asarray(sm.robust.scale.huber(X[X>min_ret]))
    mean = m[0]
    std = min(m[1],np.std(X[X>min_ret]))
    max_ret = mean + n_sigma * std

    return min_ret, max_ret