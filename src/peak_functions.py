import numpy as np
from scipy.signal import find_peaks, peak_widths

from dataclasses import dataclass

@dataclass
class PeakPars:
    peaks:     np.array
    widths:    np.array
    proms :    np.array
    left_ips:  np.array
    right_ips: np.array
    lcuts:     np.array
    rcuts:     np.array

def find_peak_params(cwf, itmn, itmx, prominence=1000, distance=20, plateau_size=0, nsigma=2):
    """
    Find peaks using Scipy function

    """
    peaks, props = find_peaks(cwf[itmn:itmx], prominence=prominence, distance=distance, plateau_size=plateau_size)
    widths, _, left_ips, right_ips = peak_widths(cwf, peaks, rel_height=0.5)
    lcuts = np.array([int(left_ips[i] -  nsigma * widths[i]) for i in range(len(widths))])
    rcuts = np.array([int(right_ips[i] +  nsigma * widths[i]) for i in range(len(widths))])
    return PeakPars(peaks, widths, props["prominences"], left_ips, right_ips,lcuts, rcuts)


def rebin_2d(xx, r):
    """
    Rebins a 2D NumPy array x of shape (n1, n2) into y of shape (n1, N2),
    where N2 = n2 // r, by averaging blocks of size r along the second axis.

    Parameters:
    -----------
    x : np.ndarray
        Input array of shape (n1, n2).
    r : int
        Rebin factor.

    Returns:
    --------
    y : np.ndarray
        Rebinned array of shape (n1, N2), obtained by averaging groups
        of 'r' elements along the second dimension.
    """
    x = np.copy(xx)
    n1, n2 = x.shape
    # Calculate the new number of columns
    N2 = n2 // r

    # Trim x so its second dimension is divisible by r
    trimmed = x[:, :N2*r]

    # Reshape so each row is split into blocks of size r
    trimmed_reshaped = trimmed.reshape(n1, N2, r)

    # Average across the last axis (the blocks of size r)
    y = trimmed_reshaped.mean(axis=2)

    return y


def rebin_sum(arr, r):
    """
    Rebin a 1D NumPy array `arr` by a factor of `r`.
    The result is the sum of each chunk of size `r`.
    If `len(arr)` is not divisible by `r`, the extra elements are discarded.
    """
    xx = np.copy(arr)
    n = xx.size
    m = n // r
    trimmed = xx[:m*r]
    return trimmed.reshape(m, r).sum(axis=1)


def s12_energy(cwf, ps):
    """
    Sum the energy of waveform cwf in the interval defined by ps

    """
    return [np.sum(cwf[ps.lcuts[i]:ps.rcuts[i]]) for i in range(len(ps.lcuts))]



def print_peak_pars(pp, tspmt, mode="wvfm"):
    if mode == "wvfm":
        print(f"found peaks -->{pp.peaks}, time position (mus): {tspmt*np.array(pp.peaks)}")
        print(f"prominences = {np.array(pp.proms)}")
        print(f"widths (mus) = {tspmt*np.array(pp.widths)}")
        print(f"left ips (mus) = {tspmt*np.array(pp.left_ips)}")
        print(f"right ips (mus) = {tspmt*np.array(pp.right_ips)}")
        print(f"left cuts = {pp.lcuts}, right cuts = {pp.rcuts}")
        print(f"left cuts (mus) = {tspmt * pp.lcuts}, right cuts = {tspmt * pp.rcuts}")
        
    else:
        print(f"found peaks -->{pp.peaks}, time position (samples): {np.array(pp.peaks)}")
        print(f"prominences = {np.array(pp.proms)}")
        print(f"widths (samples) = {np.array(pp.widths)}")
        print(f"left ips (samples) = {np.array(pp.left_ips)}")
        print(f"right ips (samples) = {np.array(pp.right_ips)}")
        print(f"left cuts = {pp.lcuts}, right cuts = {pp.rcuts}")


def compute_moving_average(w, window_size=100):
    """
    Compute the moving average of array w with a given window_size.
    The result has the same shape as w using 'same' mode in convolution.
    """
    kernel = np.ones(window_size) / window_size
    maw = np.convolve(w, kernel, mode='same')
    return maw


def apply_threshold(w, maw, th):
    """
    Keep elements of w that are above (maw + th); set others to 0.
    """
    return np.where(w > maw + th, w, 0)


def get_s2_windows(sipmwf, ps):
    SIPMW = []
    
    ipl = ps.lcuts
    ipr = ps.rcuts
    print(ipl, ipr)
    
    for i in range(len(ipl)):
        sipmw = sipmwf[:,ipl[i]:ipr[i]]
        
        print(f"peak number = {i}")
        print(f"ipl = {ipl[i]}, ipr={ipr[i]}")
        print(f"sipm window shape {sipmw.shape}")
        SIPMW.append(sipmw)
    return SIPMW


def s2_windows_sum(SIPMW, thr=3):
    QSIPM = []
    for i, sipmw in enumerate(SIPMW):
        qsipm =apply_threshold(sipmw, thr, 0)
        QSIPM.append(np.sum(qsipm, axis=1))
    return QSIPM


def sipm_xg(XSi, YSi, QSIPM):
    XG = []
    YG = []
    for i in range(len(QSIPM)):
        xsi = np.sum(XSi * QSIPM[i])/np.sum(QSIPM[i])
        ysi = np.sum(YSi * QSIPM[i])/np.sum(QSIPM[i])
        print(f"xsi = {xsi}, ysi = {ysi}")
        XG.append(xsi)
        YG.append(ysi)
    return XG,YG

