from invisible_cities.core.system_of_units import adc, pes, mus, ns
from dataclasses import dataclass

@dataclass
class StMinMax:
    stmx:  float 
    stmn:  float 
    istmx: int 
    istmn: int

irene_pars = {
    "detector_db": "next100",
    "compression": "ZLIB4",
    "n_baseline": 90000,
    "n_baseline_units": "nounits",
    "n_maw": 100,
    "n_maw_units": "nounits",
    "thr_maw": 10,
    "thr_maw_unit": "adc",
    "thr_csum_s1": 0.1,
    "thr_csum_s1_unit": "pes",
    "thr_csum_s2": 0.5,
    "thr_csum_s2_unit": "pes",
    "thr_sipm": 1.0,
    "thr_sipm_unit": "pes",
    "thr_sipm_typ": "common",
    "pmt_samp_wid": 25,
    "pmt_samp_wid_unit": "ns",
    "sipm_samp_wid": 1,
    "sipm_samp_wid_unit": "mus",
    "s1_lmax": 99,
    "s1_lmax_unit": "nounit",
    "s1_lmin": 4,
    "s1_lmin_unit": "nounit",
    "s1_tmax": 1450,
    "s1_tmax_unit": "mus",
    "s1_tmin": 0,
    "s1_tmin_unit": "mus",
    "s1_stride": 4,
    "s1_stride_unit": "nounit",
    "s1_rebin_stride": 1,
    "s1_rebin_stride_unit": "nounit",
    "s2_lmax": 10000,
    "s2_lmax_unit": "nounit",
    "s2_lmin": 400,
    "s2_lmin_unit": "nounit",
    "s2_tmax": 2000,
    "s2_tmax_unit": "mus",
    "s2_tmin": 1450,
    "s2_tmin_unit": "mus",
    "s2_stride": 10,
    "s2_stride_unit": "nounit",
    "s2_rebin_stride": 80,
    "s2_rebin_stride_unit": "nounit",
    "thr_sipm_s2": 5.0,
    "thr_sipm_s2_units": "pes"
}

fkr_globs = {"detector_db": "next100",
             "compression": "ZLIB4",
             "pmt_samp_wid_mus": 25E-3,
             "tbin_pmt_ns" : 25,
             "n_maw_samples": 100,
             "thr_maw_adcs": 10,
             "thr_sipm_calib_pes": 1.0,
             "thr_sipm_pes": 5.0}

fkr_pars = {

    "glow_peak_prominence": 1000, 
    "glow_peak_distance": 20,
    "s2_peak_prominence": 100, 
    "s2_peak_distance" : 10,
    "s1_peak_prominence": 5, 
    "s1_peak_distance" : 10,

    "s1_rebin_stride": 4,
    "s2_rebin_stride": 80,
    
    "thr_sipm_s2_mus": 5.0,

    "s1_tmin_mus": 0,
    "s1_tmax_mus": 1450,
    "s2_tmin_mus": 1450,
    "s2_tmax_mus": 2000,
    "masked_sipm" : []
}

fkr_pars_run_1463={    
    "glow_peak_prominence": 1000, 
    "glow_peak_distance": 20,
    "s2_peak_prominence": 100, 
    "s2_peak_distance" : 10,
    "s1_peak_prominence": 5, 
    "s1_peak_distance" : 10,

    "s1_rebin_stride": 4,
    "s2_rebin_stride": 80,
    
    "thr_sipm_s2_mus": 5.0,

    "s1_tmin_mus": 0,
    "s1_tmax_mus": 250,
    "s2_tmin_mus": 0,
    "s2_tmax_mus": 800,
    "masked_sipm" : [2425]
    }



def get_units(sunit):
    if sunit == "adc":
        return adc
    elif sunit == "mus":
        return mus
    elif sunit == "ns":
        return ns
    elif sunit == "pes":
        return pes
    else:
        return 1.0


def get_tmin_tmax(parameters, tbin):

    s1tmx = parameters['s1_tmax_mus'] * mus
    s1tmn = parameters['s1_tmin_mus'] * mus
    s2tmx = parameters['s2_tmax_mus'] * mus
    s2tmn = parameters['s2_tmin_mus'] * mus
    
    is1tmx = int(s1tmx/tbin)
    is1tmn = int(s1tmn/tbin)
    is2tmx = int(s2tmx/tbin)
    is2tmn = int(s2tmn/tbin)

    return StMinMax(s1tmx, s1tmn, is1tmx, is1tmn), StMinMax(s2tmx, s2tmn, is2tmx, is2tmn)


def get_maw(pars):
    n_maw = pars['n_maw'] * get_units(pars['n_maw'])
    thr_maw = pars['thr_maw'] * get_units(pars['thr_maw'])
    return int(n_maw), thr_maw

def get_thr_sipmw(pars):
    thr_sipm = pars['thr_sipm'] * get_units(pars['thr_sipm'])
    return thr_sipm

def get_var_with_units(var, pars):
    return pars[var] * get_units(pars[var])
    



    