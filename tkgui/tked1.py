import tkinter as tk
from tkinter import messagebox, filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np


import os
import sys
import time
import configparser

module_dir = os.path.abspath('../src')  # Adjust path as needed
sys.path.append(module_dir)

import tables as tb
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.signal import find_peaks, peak_widths

import irene_params as irp

from irene_functions import  wf_from_files
from aux_functions import  load_waveforms
from peak_functions import  rebin_2d, rebin_sum

from ic_functions import suppress_glow 

from peak_functions import find_peak_params, s12_energy, print_peak_pars, apply_threshold
from peak_functions import get_s2_windows, s2_windows_sum, sipm_xg

from invisible_cities.cities.components import deconv_pmt, calibrate_pmts, calibrate_sipms 

from invisible_cities.types.symbols import WfType
from invisible_cities.database.load_db import DataSiPM
from invisible_cities.database.load_db import DataPMT
from invisible_cities.core.system_of_units import adc, pes, mus, ns

from irene_params import get_s1_tmin_tmax, get_maw, get_s2_tmin_tmax, get_var_with_units

from plot_functions1 import plot_adc_to_pes
from plot_functions import plot_sum_waveform_tk, plot_waveform_zoom_peaks_tk
from plot_functions import plot_waveform_tk


# Positions of PMTs and SiPMs from Data Base
PMTPos = DataPMT('next100', 0).filter('XY')
SiPMpos = DataSiPM('next100', 0).filter('XY')
X, Y = PMTPos["X"], PMTPos["Y"]
XSi, YSi = SiPMpos["X"], SiPMpos["Y"]

class IreneParamApp:
    def __init__(self, master):
        self.master = master
        self.master.title("Irene Parameters Configuration")

        self.parameters = irp.irene_pars
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
        for i, (param, value) in enumerate(self.parameters.items()):
            label = tk.Label(self.param_frame, text=param)
            label.grid(row=i, column=0, sticky="e", padx=5, pady=2)
            entry = tk.Entry(self.param_frame, width=30)
            entry.insert(0, str(value))
            entry.grid(row=i, column=1, sticky="w", padx=5, pady=2)
            self.entries[param] = entry

        # Button frame
        self.button_frame = tk.Frame(self.master)
        self.button_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)

        # "Plot DataADC/PES" button calls external plotting logic
        #self.plot_data_button = tk.Button(self.button_frame, text="Plot ADC/PES", 
         #                                 command=self.plot_data)
        #self.plot_data_button.pack(side=tk.LEFT, padx=5)

        # "Load File" button
        self.load_file_button = tk.Button(self.button_frame, text="Load Raw Data File", 
                                          command=self.load_file)
        self.load_file_button.pack(side=tk.LEFT, padx=5)

        # "Load event" button
        self.load_events_button = tk.Button(self.button_frame, text="Load Event", 
                                            command=self.load_events, state="disabled")
        self.load_events_button.pack(side=tk.LEFT, padx=5)

        # "Plot waveforms" button 
        self.plot_cwf_button = tk.Button(self.button_frame, text="Plot WF", 
                                         command=self.plot_cwf, state="disabled")
        self.plot_cwf_button.pack(side=tk.LEFT, padx=5)

        # "Find glow" button 
        self.find_glow_button = tk.Button(self.button_frame, text="Find Glow", 
                                         command=self.find_glow, state="disabled")
        self.find_glow_button.pack(side=tk.LEFT, padx=5)

        # "Hidden" buttons: Plot glow zooms
        self.glow_plot_wvf_button   = tk.Button(self.master, text="Plot Glow", 
                                       command=self.plot_glow_wvfm)
        self.glow_zoom_peaks_button = tk.Button(self.master, text="Zoom Glow", 
                                         command=self.plot_glow_zoom_peaks)
        
        self.plot_corr_wvfm_button = tk.Button(self.master, text="Plot CWF", 
                                         command=self.plot_corr_wvfm)
       
        
        # "Print parameters" button 
        self.print_params_button = tk.Button(self.button_frame, text="Print Parameters", 
                                             command=self.print_parameters)
        self.print_params_button.pack(side=tk.LEFT, padx=5)

        # Frame to hold the matplotlib plot
        self.plot_frame = tk.Frame(self.master)
        self.plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Initialize the figure and canvas
        self.fig, self.ax = plt.subplots(figsize=(10,4))
        self.fig.clear()
        self.canvas_plot = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas_plot.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _on_frame_configure(self, event):
        # Update the scroll region of the canvas
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


    def plot_data(self):
        self.fig.clear()
        plot_adc_to_pes(self.ax)
        self.canvas_plot.draw()
    

    def plot_cwf(self):
        self.fig.clear()
        axs = self.fig.subplots(1, 2)
        #time = np.linspace(0, self.sp.pmt_time_bins * tbin,self.sp.pmt_time_bins)
        #axs.plot(time, self.cwf_sum_maw, label=f"{self.wf.run_number}", color="red")
        
        plot_sum_waveform_tk(axs, self.sipm_sum, self.cwf_sum_maw, 
                             self.sp.sipm_time_bins, self.sp.pmt_time_bins, 
                             self.wf.run_number,
                             self.wf.event_number,
                             self.st1.stmx/mus, self.st1.stmn/mus, 
                             self.st2.stmx/mus, self.st2.stmn/mus, self.tspmt)
        #plot_sum_waveform(axs, self.sipm_sum, self.cwf_sum_maw, 
        #                  self.sipm_time_bins, self.time_bins, 
        #                  self.run_number,
        #                  self.event_number)
        self.fig.tight_layout() 
        self.canvas_plot.draw()


    def plot_glow_wvfm(self):
        """
        Plot the uncorrected waveform showing all the peaks identified as glow

        """
        self.fig.clear()
        axs = self.fig.subplots(1, len(self.pp.peaks))

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
        

    def find_glow(self, prominence=1000, distance=20):
        """
        Find glow peaks, if any 
        """
        self.pp = find_peak_params(self.cwf_sum_maw, 0, len(self.cwf_sum_maw), prominence, distance)
        print_peak_pars(self.pp, self.tspmt)
        
        if len(self.pp.proms) > 0:
            self.glow_plot_wvf_button.pack(side=tk.LEFT, padx=5)
            self.glow_zoom_peaks_button.pack(side=tk.LEFT, padx=5)
            self.plot_corr_wvfm_button.pack(side=tk.LEFT, padx=5)

            self.sum_cwf_corr = suppress_glow(self.cwf_sum_maw, 
                                              self.pp.peaks, self.pp.lcuts, self.pp.rcuts)
            
        else:
            self.sum_cwf_corr = self.cwf_sum_maw
            self.plot_corr_wvfm_button.pack(side=tk.LEFT, padx=5)

    
    def print_parameters(self):

        self._update_parameters_from_entries()
        msg = "\n".join([f"{k}: {v}" for k, v in self.parameters.items()])
        print("Current Parameters:\n" + msg)
        #messagebox.showinfo("Current Parameters", msg)


    def load_file(self):
        # Open a file dialog for the user to select a file
        file_path = filedialog.askopenfilename(
            title="Select a File",
            filetypes=[("All Files", "*.*")]
        )

        if file_path:
           
            print(f"File selected: {file_path}")

            self.wfg = wf_from_files([file_path], WfType.rwf)
            self.load_events_button.config(state="normal") # activate button
            messagebox.showinfo("File Loaded", f"File selected:\n{file_path}")
            

    def get_cwfs(self):
        """
        Get deconvolved and calibrated waveforms for PMTs
        
        """

        n_maw, thr_maw = get_maw(self.parameters)
        print(f"n_maw = {n_maw}, thr_max = {thr_maw} adc")

        deconv = deconv_pmt(self.parameters["detector_db"], int(self.wf.run_number), self.wf.n_baseline)
        cwf    = deconv(self.wf.pmtrwf)

        print(f"run number = {int(self.wf.run_number)}")

        calibpmts = calibrate_pmts(self.parameters["detector_db"], int(self.wf.run_number), n_maw, thr_maw)
        #ccwfs, ccwfs_maw, cwf_sum, cwf_sum_maw   = calibpmts(cwf)

        #deconv = deconv_pmt(self.parameters["detector_db"], int(self.run_number), self.n_baseline)
        #cwf    = deconv(self.pmtrwf)
        
        #calibpmts = calibrate_pmts(self.parameters["detector_db"], int(self.run_number),
        #                           self.parameters["n_maw"], self.parameters["thr_maw"])
                                   
        self.ccwfs, self.ccwfs_maw, self.cwf_sum, self.cwf_sum_maw   = calibpmts(cwf)
        

    def get_csipm(self):
        ### Get calibrated waveforms for SiPMs

        thr_sipm = get_var_with_units("thr_sipm", self.parameters) 
        print(f"thr_sipm = {thr_sipm} pes")   

        sipm_rwf_to_cal  = calibrate_sipms(self.parameters["detector_db"], 
                                           self.wf.run_number, thr_sipm)
        sipmcwf = sipm_rwf_to_cal(self.wf.sipmrwf)


        #sipm_sum = np.sum(sipmcwf, axis=0)
        #stdsipm = np.array([np.std(sipmcwf[i][0:300]) for i in range(num_sipms)])
        #qsipmx   = np.array([np.max(sipmcwf[i]) for i in range(num_sipms)])

        #sipm_thr = self.parameters["thr_sipm"] * irp.get_units(self.parameters["thr_sipm_unit"])
        #sipm_rwf_to_cal  = calibrate_sipms(self.parameters["detector_db"], 
        #                                   int(self.run_number), sipm_thr)
        #sipmcwf = sipm_rwf_to_cal(self.sipmrwf)
        #self.num_sipms = self.sipmrwf.shape[0]
        #self.stdsipm = np.array([np.std(sipmcwf[i][0:300]) for i in range(self.num_sipms)])
        #self.qsipm   = np.array([np.max(sipmcwf[i]) for i in range(self.num_sipms)])
        self.sipm_sum = np.sum(sipmcwf, axis=0)
        self.sipmrbwf = rebin_2d(sipmcwf, 2)
    
    
    def load_events(self):
        """
        Load the a dictionary with the waveforms.
        """
        wfdct = next(self.wfg)

        #wf = Waveforms(run_number, event_number, n_baseline, pmtrwf, sipmrwf)
        #sp = SensorPars(time_bins, num_pmts, sipm_time_bins, num_sipms)
    
        self.wf, self.sp = load_waveforms(wfdct)

        # read the raw waveforms
        #self.pmtrwf = wfdct["pmt"]
        #self.sipmrwf = wfdct["sipm"]
        #self.run_number = wfdct["run_number"]
        #self.event_number = wfdct["event_number"]
        #print(f"run_number = {self.run_number}, event_number = {self.event_number}")

        #messagebox.showinfo("Event Loaded:", 
        #                    f"run number = {self.wf.run_number}, event number = {self.wf.event_number}")
        

        #self.num_pmts=self.pmtrwf.shape[0]
        #self.num_sipms = self.sipmrwf.shape[0]
        #self.time_bins = self.pmtrwf.shape[1]
        #self.sipm_time_bins = self.sipmrwf.shape[1]
        #self.n_baseline = int(0.9 * self.time_bins)
                                   
        print(f"number of PMTs = {self.sp.num_pmts}, number of PMT time bins = {self.sp.pmt_time_bins}")
        print(f"number of SiPM = {self.sp.num_sipms}, number of SiPMs time bins = {self.sp.sipm_time_bins}")
        print(f"n_baseline = {self.wf.n_baseline}")

        self.tpmt = get_var_with_units("pmt_samp_wid", self.parameters)
        self.tspmt = self.tpmt/mus
        print(f"time pmt = {self.tpmt} ns")
        print(f"time sample pmt  = {self.tspmt} mus")

        if self.wf.run_number == 14643:
            self.parameters['s1_tmax']=250
            self.parameters['s1_tmin']=0
            self.parameters['s2_tmax']=800
            self.parameters['s2_tmin']=0

            self.st1 = get_s1_tmin_tmax(self.parameters, tbin_pmt_ns = 25)

            #self.s1tmx, self.s1tmn, self.is1tmx, self.is1tmn = get_s1_tmin_tmax(ip.irene_pars, tbin_pmt_ns = 25)
            print(f"s1tmx = {self.st1.stmx/mus} mus, s1tmn = {self.st1.stmn/mus} mus")
            print(f"is1tmx = {self.st1.istmx} mus, s1tmn = {self.st1.istmn} mus")
            #self.s2tmx, self.s2tmn, self.is2tmx, self.is2tmn = get_s2_tmin_tmax(ip.irene_pars, tbin_pmt_ns = 25)

            self.st2 = get_s2_tmin_tmax(self.parameters, tbin_pmt_ns = 25)
            print(f"s2tmx = {self.st2.stmx/mus} mus, s2tmn = {self.st2.stmn/mus} mus")
            print(f"is2tmx = {self.st2.istmx}, is2tmn = {self.st2.istmn}")
        
        for i, (param, value) in enumerate(self.parameters.items()):
            entry = tk.Entry(self.param_frame, width=30)
            entry.insert(0, str(value))
            entry.grid(row=i, column=1, sticky="w", padx=5, pady=2)
            self.entries[param] = entry


        self.get_cwfs()
        self.get_csipm()

        self.plot_cwf_button.config(state="normal") # activate button
        self.find_glow_button.config(state="normal")
        
        #handle_load_events()
        messagebox.showinfo("Event Loaded:", 
                            f"run number = {self.wf.run_number}, event number = {self.wf.event_number}")


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
    root.geometry("800x800")  # Increased window size for clarity
    app = IreneParamApp(root)
    root.mainloop()
