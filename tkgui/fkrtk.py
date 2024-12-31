"""
Fast Krypton TKinter based monitor.

- Read waveforms, performs deconvolution.
- Search and suppressed glow peaks
- Finds S2 and S1 peaks
- Finds position of S2 in SiPMs.

"""
import tkinter as tk
from tkinter import messagebox, filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np


import os
import sys
#import time
#import configparser

module_dir = os.path.abspath('../src')  
sys.path.append(module_dir)

#import tables as tb
#import pandas as pd
import matplotlib.pyplot as plt
#from   matplotlib.colors import LogNorm
#import matplotlib.gridspec as gridspec
import numpy as np


from irene_params import fkr_globs, fkr_pars, fkr_pars_run_1463, masked_sipm
#import irene_params as irp

from irene_functions import  wf_from_files
from aux_functions import  load_waveforms
from peak_functions import  rebin_2d, rebin_sum, get_sipm_max

from ic_functions import suppress_glow 

from peak_functions import find_peak_params, s12_energy, print_peak_pars, apply_threshold
from peak_functions import get_s2_windows, s2_windows_sum, sipm_xg

from invisible_cities.cities.components import deconv_pmt, calibrate_pmts, calibrate_sipms 

from invisible_cities.types.symbols import WfType
from invisible_cities.database.load_db import DataSiPM
from invisible_cities.database.load_db import DataPMT
from invisible_cities.core.system_of_units import adc, pes, mus, ns

from irene_params import get_s1_tmin_tmax, get_s2_tmin_tmax

from plot_functions1 import plot_adc_to_pes
from plot_functions import plot_sum_waveform_tk, plot_waveform_zoom_peaks_tk, plot_sipmw_tk
from plot_functions import plot_waveform_tk,  plot_sipm_tk

# Globals

# Positions of PMTs and SiPMs from Data Base
PMTPos = DataPMT('next100', 0).filter('XY')
SiPMpos = DataSiPM('next100', 0).filter('XY')
X, Y = PMTPos["X"], PMTPos["Y"]
XSi, YSi = SiPMpos["X"], SiPMpos["Y"]

class FkTkApp:
    def __init__(self, master):
        self.master = master
        self.master.title("Fast Krypton Configuration")

        self.parameters = fkr_pars
        self.entries = {}

        # Create a frame for the scrollable area
        self.scroll_frame = tk.Frame(self.master)
        self.scroll_frame.pack(fill=tk.BOTH, expand=True)

        # Create a canvas and a scrollbar for parameters
        self.canvas = tk.Canvas(self.scroll_frame)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.scrollbar = tk.Scrollbar(self.scroll_frame, orient=tk.VERTICAL, command=self.canvas.yview)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        # Create a frame inside the canvas for the parameters
        self.param_frame = tk.Frame(self.canvas)
        self.canvas.create_window((0,0), window=self.param_frame, anchor="nw")

        # Bind event to update scroll region
        self.param_frame.bind("<Configure>", self._on_frame_configure)

        # Populate param_frame with entries

        self.PARS = []
        self.VALS = []
        self.nc =3       # number of columns to arrange the parameters 

        for (param, value) in self.parameters.items():
            self.PARS.append(param)
            self.VALS.append(value)
        
        self.nr =int(len(self.PARS)/self.nc)+1

        ip = 0
        for i in range(self.nr):
            k=0
            for j in range(self.nc):
                if ip < len(self.PARS):
                    label = tk.Label(self.param_frame, text= self.PARS[ip])
                    label.grid(row=i, column=j+k, sticky="e", padx=5, pady=2)
                
                    entry = tk.Entry(self.param_frame, width=10)
                    entry.insert(0, str( self.VALS[ip]))
                    entry.grid(row=i, column=j+k+ 1, sticky="w", padx=5, pady=2)
                    self.entries[self.PARS[ip]] = entry
                    k = k+1
                ip+=1
      
        # for (param, value) in enumerate(self.parameters.items()):
        #     label = tk.Label(self.param_frame, text=param)
        #     label.grid(row=i, column=j%2, sticky="e", padx=5, pady=2)
        #     entry = tk.Entry(self.param_frame, width=10)
        #     entry.insert(0, str(value))
        #     entry.grid(row=i, column=1, sticky="w", padx=5, pady=2)
        #     self.entries[param] = entry


        # Label for event number
        self.label_event_num = tk.Label(self.param_frame, text="Event Number:")
        self.label_event_num.grid(row=i+1, column=0, sticky="e", padx=5, pady=2)

        # Entry widget to accept the event number
        self.event_number_entry = tk.Entry(self.param_frame, width=5)
        self.event_number_entry.insert(0, -1)
        self.event_number_entry.grid(row=i+1, column=1, sticky="w", padx=5, pady=2)
        
        # Button frame
        self.button_frame = tk.Frame(self.master)
        self.button_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)

        # "Plot DataADC/PES" button calls external plotting logic
        #self.plot_data_button = tk.Button(self.button_frame, text="Plot ADC/PES", 
         #                                 command=self.plot_data)
        #self.plot_data_button.pack(side=tk.LEFT, padx=5)

        # "Load File" button
        self.load_file_button = tk.Button(self.button_frame, text="Read Data File", 
                                          command=self.load_file)
        self.load_file_button.grid(row=0, column=0, sticky="e", padx=5, pady=2)

        # "Load event" button
        self.load_events_button = tk.Button(self.button_frame, text="Load Event", 
                                            command=self.load_events, state="disabled")
        self.load_events_button.grid(row=0, column=1, sticky="e", padx=5, pady=2) 

        # "Print parameters" button 
        self.print_params_button = tk.Button(self.button_frame, text="Print Parameters", 
                                             command=self.print_parameters)
        self.print_params_button.grid(row=0, column=2, sticky="e", padx=5, pady=2)      
    
        # "Plot waveforms" button 
        self.plot_cwf_button = tk.Button(self.button_frame, text="Plot WFS", 
                                         command=self.plot_cwf, state="disabled")
        self.plot_cwf_button.grid(row=1, column=0, sticky="e", padx=5, pady=2)

        # "Find glow" button 
        self.find_glow_button = tk.Button(self.button_frame, text="Find Glow", 
                                         command=self.find_glow, state="disabled")
        self.find_glow_button.grid(row=1, column=1, sticky="e", padx=5, pady=2)

        #  Plot glow, glow zooms and corrected waveform
        self.glow_plot_wvf_button   = tk.Button(self.button_frame, text="Plot Glow", 
                                                command=self.plot_glow_wvfm, state="disabled")
        self.glow_plot_wvf_button.grid(row=1, column=2, sticky="e", padx=5, pady=2)

        self.glow_zoom_peaks_button = tk.Button(self.button_frame, text="Zoom Glow", 
                                                command=self.plot_glow_zoom_peaks, state="disabled")
        self.glow_zoom_peaks_button.grid(row=1, column=3, sticky="e", padx=5, pady=2)
        
        self.plot_corr_wvfm_button = tk.Button(self.button_frame, text="Plot CWF", 
                                                command=self.plot_corr_wvfm, state="disabled")
        self.plot_corr_wvfm_button.grid(row=1, column=4, sticky="e", padx=5, pady=2)
       
       # S2 search
        self.s2_search_button   = tk.Button(self.button_frame, text="S2 Search", 
                                                command=self.s2_search, state="disabled")
        self.s2_search_button.grid(row=2, column=0, sticky="e", padx=5, pady=2)

        # S1 search
        self.s1_search_button   = tk.Button(self.button_frame, text="S1 Search", 
                                                command=self.s1_search, state="disabled")
        self.s1_search_button.grid(row=2, column=1, sticky="e", padx=5, pady=2)

        # Plot SiPMs
        self.plot_sipm_button   = tk.Button(self.button_frame, text="Plot SiPMs", 
                                                command=self.plot_sipm, state="disabled")
        self.plot_sipm_button.grid(row=3, column=0, sticky="e", padx=5, pady=2)

        # Get SiPMs in S2 window(s)
        self.get_sipm_s2_window_button   = tk.Button(self.button_frame, text="Get SiPMs in S2 window", 
                                                command=self.get_sipm_s2_window, state="disabled")
        self.get_sipm_s2_window_button.grid(row=3, column=1, sticky="e", padx=5, pady=2)

        # SiPMs in S2 sum in window(s)
        self.get_sipm_s2_window_sum_button   = tk.Button(self.button_frame, text="Sum SiPMs in S2 window", 
                                                command=self.get_sipm_s2_window_sum, state="disabled")
        self.get_sipm_s2_window_sum_button.grid(row=3, column=2, sticky="e", padx=5, pady=2)

        # SiPMs Baricenter
        self.get_sipm_baricenter_button   = tk.Button(self.button_frame, text="Compute Baricenter", 
                                                command=self.compute_baricenter, state="disabled")
        self.get_sipm_baricenter_button.grid(row=3, column=3, sticky="e", padx=5, pady=2)

        
        # Frame to hold the matplotlib plot
        self.plot_frame = tk.Frame(self.master)
        self.plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Initialize the figure and canvas
        self.fig, self.ax = plt.subplots(figsize=(12,6))
        self.fig.clear()
        self.canvas_plot = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas_plot.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    

    def load_file(self):
        """
        Open a file dialog for the user to select a file

        """ 
        file_path = filedialog.askopenfilename(
            title="Select a File",
            filetypes=[("All Files", "*.*")]
        )

        if file_path:
            print(f"File selected: {file_path}")
            self.current_event=0
            self.file_path = file_path
            self.wfg = wf_from_files([file_path], WfType.rwf)
            messagebox.showinfo("File Loaded", f"File selected:\n{file_path}")

            # Activate load events futtom
            self.load_events_button.config(state="normal") # activate button

            # Deactivate all other buttons:
            self.plot_cwf_button.config(state="disabled") # activate button to plot calibrated waveforms
            self.find_glow_button.config(state="disabled") # activate button to find glow.
            self.glow_plot_wvf_button.config(state="disabled")
            self.glow_zoom_peaks_button.config(state="disabled")
            self.plot_corr_wvfm_button.config(state="disabled")
            self.s2_search_button.config(state="disabled")
            self.s1_search_button.config(state="disabled")
            self.plot_sipm_button.config(state="disabled")
            self.get_sipm_s2_window_button.config(state="disabled")
            self.get_sipm_s2_window_sum_button.config(state="disabled")
            self.get_sipm_baricenter_button.config(state="disabled")
        
           
    def load_events(self):
        """
        Load the a dictionary with the waveforms from file.
        Read Waveforms and Sensor Parameters, packed in two data classes

            wf = Waveforms(run_number, event_number, n_baseline, pmtrwf, sipmrwf)
            sp = SensorPars(time_bins, num_pmts, sipm_time_bins, num_sipms)
        
        Redefine tmin, tmax for s1, s2 in the case of run 14643 (800 mus buffer)
        Get parameters for tmin, tmax (s1, s2) packed in the data class. 

        StMinMax(s1tmx, s1tmn, is1tmx, is1tmn)
        StMinMax(s2tmx, s2tmn, is2tmx, is2tmn)

        Get Waveforms (corrected) for PMTs and SiPMs.
    
        """
        event_number = int(self.event_number_entry.get())

        print(f"Reading local event number = {event_number}")

        if event_number < 0:
            events_to_skip = 0
            self.current_event = self.current_event+1

        elif event_number >=  self.current_event: # skip the difference between actual and current
            events_to_skip = event_number - self.current_event
            self.current_event = event_number
        else: # we need to rewind
            self.wfg = wf_from_files([self.file_path], WfType.rwf)
            events_to_skip = event_number
            self.current_event = event_number
        
        
        for i in range(events_to_skip): # read but do not process event up to event_number
            wfdct = next(self.wfg)  
        
        wfdct = next(self.wfg) # read waveform dict for event number
        
        self.wf, self.sp = load_waveforms(wfdct)  # get waveform and sensor parameters
        self.tspmt =fkr_globs["pmt_samp_wid_mus"]  # time conversion for PMTs

        print(f" Event loaded: run number = {self.wf.run_number}, event number = {self.wf.event_number}")    
        print(f"number of PMTs = {self.sp.num_pmts}, number of PMT time bins = {self.sp.pmt_time_bins}")
        print(f"number of SiPM = {self.sp.num_sipms}, number of SiPMs time bins = {self.sp.sipm_time_bins}")
        print(f"n_baseline = {self.wf.n_baseline}")
        print(f"time sample pmt  = {self.tspmt} mus")

        if self.wf.run_number == 14643:
            self.parameters = fkr_pars_run_1463
            self.PARS = []
            self.VALS = []

            for (param, value) in self.parameters.items():
                self.PARS.append(param)
                self.VALS.append(value)
           
            ip = 0
            for i in range(self.nr):
                k=0
                for j in range(self.nc):
                    if ip < len(self.PARS):
                        entry = tk.Entry(self.param_frame, width=10)
                        entry.insert(0, str(self.VALS[ip]))
                        entry.grid(row=i, column=j+k+ 1, sticky="w", padx=5, pady=2)
                        self.entries[self.PARS[ip]] = entry
                        k = k+1
                    ip+=1

           

            # Label for event number
            self.label_event_num = tk.Label(self.param_frame, text="Event Number:")
            self.label_event_num.grid(row=i+1, column=0, sticky="e", padx=5, pady=2)

            # Entry widget to accept the event number
            self.event_number_entry = tk.Entry(self.param_frame, width=5)
            self.event_number_entry.insert(0, -1)
            self.event_number_entry.grid(row=i+1, column=1, sticky="w", padx=5, pady=2)

        tbin=   fkr_globs["tbin_pmt_ns"] 
        self.st1 = get_s1_tmin_tmax(self.parameters, tbin)
        self.st2 = get_s2_tmin_tmax(self.parameters, tbin)

        print(f"s1tmx = {self.st1.stmx/mus} mus, s1tmn = {self.st1.stmn/mus} mus")
        print(f"is1tmx = {self.st1.istmx}, s1tmn = {self.st1.istmn}")
        print(f"s2tmx = {self.st2.stmx/mus} mus, s2tmn = {self.st2.stmn/mus} mus")
        print(f"is2tmx = {self.st2.istmx}, is2tmn = {self.st2.istmn}")

        self.get_cwfs()     # Calibrated (deconvoluted) waveform for PMTs
        self.get_csipm()    # Calibrated waveforms for SiPMs

        self.plot_cwf_button.config(state="normal") # activate button to plot calibrated waveforms
        self.find_glow_button.config(state="normal") # activate button to find glow.
        self.glow_plot_wvf_button.config(state="disabled")
        self.glow_zoom_peaks_button.config(state="disabled")
        
        messagebox.showinfo("Event Loaded:", 
                            f"run number = {self.wf.run_number}, event number = {self.wf.event_number}")


    def get_cwfs(self):
        """
        Get deconvolved and calibrated waveforms for PMTs
        
        """

        n_maw =   fkr_globs["n_maw_samples"]
        thr_maw = fkr_globs["thr_maw_adcs"]
        print(f"n_maw = {n_maw} samples, thr_max = {thr_maw} adcs")

        # deconvolution
        deconv = deconv_pmt(fkr_globs["detector_db"], int(self.wf.run_number), self.wf.n_baseline)
        cwf    = deconv(self.wf.pmtrwf)

        # calibrated waveforms
        calibpmts = calibrate_pmts(fkr_globs["detector_db"], int(self.wf.run_number), n_maw, thr_maw)
        
        # Output: calibrated waveforms, cwf with mau, calibrated sum, cal sum with mau                           
        self.ccwfs, self.ccwfs_maw, self.cwf_sum, self.cwf_sum_maw   = calibpmts(cwf)
        

    def get_csipm(self):
        """
        Get calibrated waveforms for SiPMs
        
        """
      
        thr_sipm = fkr_globs["thr_sipm_calib_pes"]  
        print(f"thr_sipm = {thr_sipm} pes")   

        sipm_rwf_to_cal  = calibrate_sipms(fkr_globs["detector_db"], 
                                           self.wf.run_number, thr_sipm)
        sipmcwf = sipm_rwf_to_cal(self.wf.sipmrwf)

        self.sipm_sum = np.sum(sipmcwf, axis=0)  # sum
        self.sipmrbwf = rebin_2d(sipmcwf, 2)     # rebinned waveform at 2 mus


        # msipm = []
        # for run_number, value in masked_sipm.items():
        #  if run_number == self.wf.run_number:
        #     msipm = value
        #     break         
            
        if len(masked_sipm) > 0:
            self.sipmrbwf[masked_sipm] = 0
    

    def find_glow(self):
        """
        Find glow peaks, if any 
        Returns an instance of Peak Parameters data class

        PeakPars(peaks, widths, props["prominences"], left_ips, right_ips,lcuts, rcuts)
        """

        prominence= self.parameters["glow_peak_prominence"]
        distance  = self.parameters["glow_peak_distance"]

        self.pp = find_peak_params(self.cwf_sum_maw, 0, len(self.cwf_sum_maw), 
                                   prominence, distance, nsigma=3)
        print_peak_pars(self.pp, self.tspmt)
        
        if len(self.pp.proms) > 0:

            self.glow_plot_wvf_button.config(state="normal")
            self.glow_zoom_peaks_button.config(state="normal")
            self.plot_corr_wvfm_button.config(state="normal")

            self.sum_cwf_corr = suppress_glow(self.cwf_sum_maw, 
                                              self.pp.peaks, self.pp.lcuts, self.pp.rcuts)
            
        else:
            self.sum_cwf_corr = self.cwf_sum_maw
            self.glow_plot_wvf_button.config(state="disabled")
            self.glow_zoom_peaks_button.config(state="disabled")
            self.plot_corr_wvfm_button.config(state="normal")
        
        self.s2_search_button.config(state="normal")
        self.s1_search_button.config(state="normal")
        self.plot_sipm_button.config(state="normal")
        

    def s2_search(self):
        """
        Rebin s2
        Search for S2 peaks

        """

        s2_rebin = int(self.parameters["s2_rebin_stride"])
        self.cwf_s2 = rebin_sum(self.sum_cwf_corr, s2_rebin)
        self.tbins2 = self.tspmt * s2_rebin

        print(f"rebined s2 bins (mus) = {self.tbins2}")
        print(f"cwf_2 lengths = {self.cwf_s2.shape[0]}")

        prominence= self.parameters["s2_peak_prominence"]
        distance  = self.parameters["s2_peak_distance"]

        tbin=   int(fkr_globs["tbin_pmt_ns"] *  s2_rebin)
        self.st2 = get_s2_tmin_tmax(self.parameters, tbin)

        print(f"Searching S2 with s2 tmin = {self.st2.stmn/mus}, s2 tmax = {self.st2.stmx/mus}")
        print(f"index s2 tmin = {self.st2.istmn}, index s1 tmax = {self.st2.istmx}")
        print(f"prominence = {prominence}, distance to next peak = {distance}")

        self.ps2 = find_peak_params(self.cwf_s2, self.st2.istmn, self.st2.istmx, prominence, distance)
        print_peak_pars(self.ps2, self.tbins2)

        print(f"Energy of S2 = {s12_energy(self.cwf_s2, self.ps2)} pes")

        self.plot_s2_wvfm()
        # activate buttom to search SiPMs on window S2
        self.get_sipm_s2_window_button.config(state="normal")

    
    def s1_search(self):
        """
        Rebin s1
        Search for S1 peaks

        """
        self.s1_rebin = int(self.parameters["s1_rebin_stride"])
        self.cwf_s1 = rebin_sum(self.sum_cwf_corr, self.s1_rebin)
        self.tbins1 = self.tspmt * self.s1_rebin

        print(f"rebined s1 bins (mus) = {self.tbins1}")
        print(f"cwf_s1 length = {self.cwf_s1.shape[0]}")

        tbin=   int(fkr_globs["tbin_pmt_ns"] *  self.s1_rebin)
        self.st1 = get_s1_tmin_tmax(self.parameters, tbin)

        prominence= self.parameters["s1_peak_prominence"]
        distance  = self.parameters["s1_peak_distance"]
        

        print(f"Searching S1 with s1 tmin = {self.st1.stmn/mus}, s1 tmax = {self.st1.stmx/mus}")
        print(f"index s1 tmin = {self.st1.istmn}, index s1 tmax = {self.st1.istmx}")
        print(f"prominence = {prominence}, distance to next peak = {distance}")

        self.ps1 = find_peak_params(self.cwf_s1, self.st1.istmn, self.st1.istmx, prominence, distance)
        print_peak_pars(self.ps1, self.tbins1)

        print(f"Energy of S1 = {s12_energy(self.cwf_s1, self.ps1)} pes")

        self.plot_s1_wvfm()


    def get_sipm_s2_window(self):
        """
        Get SiPMs in S2 window(s)
        """
        self.SIPMW = get_s2_windows(self.sipmrbwf, self.ps2)

        self.fig.clear()
        axs = self.fig.subplots(1, len(self.SIPMW))

        plot_sipmw_tk(axs,self.SIPMW)
        
        self.fig.tight_layout() 
        self.canvas_plot.draw()

        self.get_sipm_s2_window_sum_button.config(state="normal")


    def get_sipm_s2_window_sum(self):
        """
        Sum the signals of SiPMs in the windows (with a threshold) 
        """

        thr = self.parameters["thr_sipm_s2_pes"]
        self.QSIPM = s2_windows_sum(self.SIPMW, thr)

        _ = get_sipm_max(self.QSIPM)

        self.fig.clear()
        axs = self.fig.subplots(1, len(self.QSIPM))
        plot_sipmw_tk(axs,self.QSIPM)
        
        self.fig.tight_layout() 
        self.canvas_plot.draw()

        self.get_sipm_baricenter_button.config(state="normal")

        
    def compute_baricenter(self):
        """
        Compute Baricenter

        """
        XG, YG = sipm_xg(XSi, YSi, self.QSIPM)

        self.fig.clear()
        axs = self.fig.subplots(1, len(self.QSIPM))

        plot_sipm_tk(axs,XSi,YSi,self.QSIPM, XG, YG)
        self.canvas_plot.draw()
        

    def plot_sipm(self):
        """
        Plot SiPMs

        """
        self.fig.clear()
        axs = self.fig.subplots(1, 1)

        plot_sipmw_tk(axs,[self.sipmrbwf])
        
        self.fig.tight_layout() 
        self.canvas_plot.draw()


    def plot_s2_wvfm(self):
        """
        Plot s2 rebinned with peak search

        """
        self.fig.clear()
        axs = self.fig.subplots(1, 1)

        plot_waveform_tk(axs, self.cwf_s2, self.cwf_s2.shape[0], 
                             self.wf.run_number, self.wf.event_number, 
                             self.ps2.peaks, self.ps2.widths, self.ps2.left_ips, self.ps2.right_ips,
                             tbin=self.tbins2)
        
        self.fig.tight_layout() 
        self.canvas_plot.draw()
    

    def plot_s1_wvfm(self):
        """
        Plot s1 rebinned with peak search

        """
        self.fig.clear()
        axs = self.fig.subplots(1, 1)
        imn = int(self.st1.istmn)
        imx = int(self.st1.istmx)
        id = imx - imn

        print(f"imx = {imx}, imn = {imn}, id = {id}, tbins s1 = {self.tbins1}")

        plot_waveform_tk(axs, self.cwf_s1[imn:imx], id,  
                         self.wf.run_number, self.wf.event_number, 
                         self.ps1.peaks, self.ps1.widths, self.ps1.left_ips, self.ps1.right_ips,
                         tbin=self.tbins1)
        
        self.fig.tight_layout() 
        self.canvas_plot.draw()


    def plot_cwf(self):
        """
        Plot the calibrated sums for PMTs and SiPMs
        
        """
        self.fig.clear()
        axs = self.fig.subplots(1, 2)
        
        plot_sum_waveform_tk(axs, self.sipm_sum, self.cwf_sum_maw, 
                             self.sp.sipm_time_bins, self.sp.pmt_time_bins, 
                             self.wf.run_number,
                             self.wf.event_number,
                             self.st1.stmx/mus, self.st1.stmn/mus, 
                             self.st2.stmx/mus, self.st2.stmn/mus, self.tspmt)
        
        self.fig.tight_layout() 
        self.canvas_plot.draw()


    def plot_glow_wvfm(self):
        """
        Plot the uncorrected waveform showing all the peaks identified as glow

        """
        self.fig.clear()
        axs = self.fig.subplots(1, 1)

        plot_waveform_tk(axs, self.cwf_sum_maw, self.sp.pmt_time_bins, 
                             self.wf.run_number, self.wf.event_number, 
                             self.pp.peaks, self.pp.widths, self.pp.left_ips, self.pp.right_ips,
                             tbin=self.tspmt)
        
        self.fig.tight_layout() 
        self.canvas_plot.draw()


    def plot_glow_zoom_peaks(self):
        """
        Plot a zoom of the peaks identified as glow.
        """
        self.fig.clear()
        axs = self.fig.subplots(1, len(self.pp.peaks))
        plot_waveform_zoom_peaks_tk(axs, self.cwf_sum_maw, self.wf.run_number, self.wf.event_number, 
                                     self.pp.peaks, self.pp.lcuts, self.pp.rcuts, twindows =self.pp.widths, 
                                     tbin=self.tspmt,
                                     tscale=True)
        
        self.fig.tight_layout() 
        self.canvas_plot.draw()
        
    
    def plot_corr_wvfm(self):
        """
        Plot the waveform corrected after suppressing glow 

        """
        self.fig.clear()
        axs = self.fig.subplots(1, 1)

        plot_waveform_tk(axs, self.sum_cwf_corr, self.sp.pmt_time_bins, 
                             self.wf.run_number, self.wf.event_number, 
                             peaks=[], widths=[], left_ips=[], right_ips=[],
                             tbin=self.tspmt)
        
        self.fig.tight_layout() 
        self.canvas_plot.draw()
            
    
    def print_parameters(self):
        """
        Print current parameters
        """

        self._update_parameters_from_entries()
        msg = "\n".join([f"{k}: {v}" for k, v in self.parameters.items()])
        print("Current Parameters:\n" + msg)
        #messagebox.showinfo("Current Parameters", msg)


    def plot_data(self):
        """
        Plot histogram of data to adc
        """
        self.fig.clear()
        plot_adc_to_pes(self.ax)
        self.canvas_plot.draw()
 
    
    def _on_frame_configure(self, event):
        """
        Update the scroll region of the canvas
        
        """
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))


    def show_parameters(self):
        """
        Fetch the current values from each Entry widget
        and display them (converted to float or int if possible).

        """
        updated_parameters = {}
        for param_name, entry_widget in self.param_entries.items():
            raw_value = entry_widget.get()
            # Attempt to convert to float or int; fallback to string if it fails
            try:
                if "." in raw_value:
                    value = float(raw_value)
                else:
                    value = int(raw_value)
            except ValueError:
                value = raw_value
            updated_parameters[param_name] = value
            print(updated_parameters)


    def _update_parameters_from_entries(self):
        for param, entry in self.entries.items():
            val = entry.get()
            # Attempt numeric conversion
            try:
                if '.' in val:
                    val = float(val)
                else:
                    val = int(val)
            except ValueError:
                # Keep as string if conversion fails
                pass
            self.parameters[param] = val


if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("1200x1200")  # Increased window size for clarity
    app = FkTkApp(root)
    root.mainloop()
