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

def find_peak_params(cwf, prominence=1000, distance=20, plateau_size=0, nsigma=2):
    """
    Find peaks using Scipy function

    """
    peaks, props = find_peaks(cwf, prominence=prominence, distance=distance, plateau_size=plateau_size)
    widths, _, left_ips, right_ips = peak_widths(cwf, peaks, rel_height=0.5)
    lcuts = [int(left_ips[i] -  nsigma * widths[i]) for i in range(len(widths))]
    rcuts = [int(right_ips[i] +  nsigma * widths[i]) for i in range(len(widths))]
    return PeakPars(peaks, widths, props["prominences"], left_ips, right_ips,lcuts, rcuts)


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


def print_peak_pars(pp, tspmt):
    print(f"found peaks -->{pp.peaks}, time position (mus): {tspmt*np.array(pp.peaks)}")
    print(f"prominences = {np.array(pp.proms)}")
    print(f"widths (mus) = {tspmt*np.array(pp.widths)}")
    print(f"left ips (mus) = {tspmt*np.array(pp.left_ips)}")
    print(f"right ips (mus) = {tspmt*np.array(pp.right_ips)}")
    print(f"left cuts = {pp.lcuts}, right cuts = {pp.rcuts}")
