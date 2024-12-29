import os
import tables as tb
import pandas as pd
import numpy as np 

from dataclasses import dataclass

@dataclass
class SensorPars:
    pmt_time_bins: float
    num_pmts: int
    sipm_time_bins: float
    num_sipms: int


@dataclass
class Waveforms:
    run_number:   float
    event_number: float
    n_baseline:   int 
    pmtrwf:       np.array
    sipmrwf:      np.array
    
@dataclass
class PmtCWF:
    ccwfs:      np.array 
    ccwfs_maw:   np.array 
    cwf_sum:     np.array 
    cwf_sum_maw: np.array

def check_folder_exists(run_number, base_path, path_to_wf):
    """
    Check if folder exists
    
    """
    folder_path=base_path+"/"+str(run_number)+"/"+ path_to_wf

    if os.path.exists(folder_path):
        filelist = sorted([f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))])
        file_count = len(filelist)


        if file_count>0:
            outcomand = os.popen(f'ptdump {folder_path}{filelist[0]} | grep "pmtrwf" | grep -m 7 -o "[0-9][0-9]" | grep -m 1 "[0-9][0-9]"').read()
        else:
            return False, 0, 0
            
        nevents_per_file=int(outcomand)
        return folder_path, filelist, nevents_per_file
            
    else:
        print(f"Run does not exist")
        return False, 0, 0


def get_filename_event(filelist, folder_path, nevents_per_file, event_number):
    """
    Return the filename and local event number
    
    """
    
    filenum_inlist = event_number // nevents_per_file
    local_ev_number = event_number % nevents_per_file
    filename = os.path.join(folder_path, filelist[filenum_inlist])
    return filename, local_ev_number 


def read_waveforms(filename, local_ev_number, xbaseline=0.9):
    """
    Return the raw waveforms from PMTs and SiPMs
    
    """
    
    try:
        f = tb.open_file(filename, 'r')  # Initialize `f` and open the file  
        pmtrwf_full= f.get_node('/RD/pmtrwf')[:]
        sipmrwf_full=f.get_node('/RD/sipmrwf')[:]
    except NameError:
        print(f"filename ={filename} does not exist")
        f = tb.open_file(filename, 'r')  # Initialize `f` and open the file
        
    
    pmtrwf = pmtrwf_full[local_ev_number, :, :]
    sipmrwf = sipmrwf_full[local_ev_number, :, :]

    time_bins = pmtrwf.shape[1]
    sipm_time_bins = sipmrwf.shape[1]
    n_baseline = int(xbaseline * time_bins)
    num_pmts=pmtrwf.shape[0]
    num_sipms = sipmrwf.shape[0]
    
    return pmtrwf, sipmrwf, n_baseline, SensorPars(time_bins, num_pmts, sipm_time_bins, num_sipms)
    

def load_waveforms(wfdct):
    """
    Given a dictionary of waveforms load the waveforms for PMTs and SiPMs

    """

# read the raw waveforms
    pmtrwf = wfdct["pmt"]
    sipmrwf = wfdct["sipm"]
    run_number = wfdct["run_number"]
    event_number = wfdct["event_number"]
    
    print(f"run_number = {run_number}, event_number = {event_number}")

    num_pmts=        pmtrwf.shape[0]
    num_sipms =      sipmrwf.shape[0]
    time_bins =      pmtrwf.shape[1]
    sipm_time_bins = sipmrwf.shape[1]
    n_baseline =     int(0.9 * time_bins)
                                
    print(f"number of PMTs = {num_pmts}, number of PMT time bins = {time_bins}")
    print(f"number of SiPM = {num_sipms}, number of SiPMs time bins = {sipm_time_bins}")
    print(f"n_baseline = {n_baseline}")

    wf = Waveforms(run_number, event_number, n_baseline, pmtrwf, sipmrwf)
    sp = SensorPars(time_bins, num_pmts, sipm_time_bins, num_sipms)
    return wf, sp 
    


def make_temp_file(config_file_path, unshown_var, variables, units):
    
    with open(config_file_path, "w") as config_file:
        # Write unshown_var contents
        for key, value in unshown_var.items():
            # Add quotes around string values
            if isinstance(value, str):
                config_file.write(f'{key} = "{value}"\n')
            else:
                config_file.write(f"{key} = {value}\n")

        for key, value in variables.items():
            if key in units:
                config_file.write(f'{key} = {value} {units[key]}\n')
            else:
                config_file.write(f'{key} = {value}\n')
                    

def launch_reco_irene(config_file_path, filename, fout, unshown_var, variables, units):
    
    make_temp_file(config_file_path, unshown_var, variables, units)
    if os.path.exists(fout):
        os.system(f"rm {fout}")
        
    cmd = f"city irene {config_file_path} -i {filename} -o {fout}"    

    os.system(cmd+"> /dev/null 2>&1")



def read_irene(irene_file):
    
    f_Irene = tb.open_file(irene_file, 'r')
    
    df_S1_all = pd.DataFrame(f_Irene.get_node('/PMAPS/S1')[:])
    df_S2_all = pd.DataFrame(f_Irene.get_node('/PMAPS/S2')[:])
    df_S1_sing_PMT = pd.DataFrame(f_Irene.get_node('/PMAPS/S1Pmt')[:])
    df_S2_sing_PMT = pd.DataFrame(f_Irene.get_node('/PMAPS/S2Pmt')[:])
    df_S2_sing_SiPM = pd.DataFrame(f_Irene.get_node('/PMAPS/S2Si')[:])
    unique_SiPM_list=np.sort(df_S2_sing_SiPM['nsipm'].unique())
    return df_S1_all,df_S2_all, df_S1_sing_PMT, df_S2_sing_PMT, df_S2_sing_SiPM, unique_SiPM_list


def launch_reco_sophronia(config_file_path, filename, fout, unshown_var, variables, units):
    
    make_temp_file(config_file_path, unshown_var, variables, units)
    if os.path.exists(fout):
        os.system(f"rm {fout}")
        
    cmd = f"city irene {config_file_path} -i {filename} -o {fout}"    

    os.system(cmd+"> /dev/null 2>&1")
