import numpy as np
from invisible_cities.cities.components import deconv_pmt 
from invisible_cities.database.load_db import DataSiPM
from invisible_cities.database.load_db import DataPMT
from invisible_cities.database import load_db


def get_adc_to_pes():
    DataPMTX    = load_db.DataPMT("next100", run_number = 14637)
    adc_to_pes = np.abs(DataPMTX.adc_to_pes.values)
    adc_to_pes = adc_to_pes[adc_to_pes > 0]
    adc_to_pes[(adc_to_pes > 1e+3)] = 0
    return adc_to_pes


def deconvolve(run_number, pmtrwf, deconv_par):
    """
    Calls IC deconvolve
    
    """
    deconv = deconv_pmt("next100", int(run_number), deconv_par)
    return deconv(pmtrwf)



def suppress_glow(cwf, peaks, lwindow, rwindow):
    """
    Return a wvf where the glow peaks are suppressed.
    
    """
    cwf2 = np.copy(cwf)
    for i, pk in enumerate(peaks):
        
        xmin = lwindow[i]
        xmax = rwindow[i]
        cwf2[xmin:xmax]=0
    
    return cwf2


def find_dead_pmts(qpmt, stds,  num_pmts, nstd=7):
    """
    Return indices of dead pmts
    """

    dpmts = []
    for i in range(num_pmts):
        if qpmt[i] < nstd * stds[i] :
            print(f" DEAD PMT number = {i}, Q = {qpmt[i]}, std = {stds[i]}")
            dpmts.append(i)
    return dpmts


def corrected_sum_without_dead_pmts(cwf_corr, dpmts):
    """
    Return the sum of corrected waveforms not in dead pmt list
    
    """
    cwf_corr_sum = np.zeros(cwf_corr[0].shape)
    for i in range(len(cwf_corr)):
        if i not in dpmts:
            cwf_corr_sum += cwf_corr[i]
    return cwf_corr_sum


def active_pmts(qpmt, dpmts):
    apmts =[]
    for i in range(len(qpmt)):
        if i not in dpmts:
            apmts.append(qpmt[i])
    return apmts


def sgn_sipms(qsipm, nsigma=3):
    sipmMax   = np.max(qsipm) 
    sipmMean  = np.mean(qsipm)
    sipmStd   = np.std(qsipm)
    print(f" sipm charge: max = {sipmMax} mean = {sipmMean} std = {sipmStd}")
    condition = qsipm > sipmMean + nsigma * sipmStd
    qsipm2 = np.where(condition, qsipm, 0)
    return qsipm2
