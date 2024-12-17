import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import numpy as np

def plot_sum_waveform(wvf, cwf, n_time_bins, label="Waveform",
                      twvf="Sum Raw Waveform", tcwf="Sum Deconv Waveform",
                      figsize=(18, 6), tbin=25e-3): 
    """
    Plots the raw and deconv sum waveforms. Time in ns
    """
    time = np.linspace(0, n_time_bins * tbin, n_time_bins)  
    
    fig, axs = plt.subplots(1, 2, figsize=figsize)
    axs[0].plot(time, wvf, label=label, color="blue")
    axs[0].set_title(twvf)
    axs[0].set_xlabel("Time [μs]")
    axs[0].set_ylabel("Amplitude")
    #axs[0].set_yscale('log')
    axs[0].grid(True)
    axs[0].legend()
    axs[1].plot(time, cwf, label=label, color="red")
    axs[1].set_title(tcwf)
    axs[1].set_xlabel("Time [μs]")
    axs[1].set_ylabel("Amplitude")
    #axs[0].set_yscale('log')
    axs[1].grid(True)
    axs[1].legend()


def plot_waveform(wvf, n_time_bins, peaks=None, window= 20, label="Waveform", 
                  log=False, figsize=(18, 6), tbin=25e-3): 
    """
    Plot a waveform. Optionaly mark the position of peaks
    
    """
    
    fig, axs = plt.subplots(1, 1, figsize=figsize)
    time = np.linspace(0, n_time_bins * tbin, n_time_bins)  
    axs.plot(time, wvf, label=label, color="blue")

    if len(peaks)>0:
        tmin = [(pk - window)*tbin for pk in peaks]
        tmax = [(pk + window)*tbin for pk in peaks]
        for i in range(len(tmin)):
            axs.axvline(tmin[i], color='red', linestyle='--', linewidth=2)
        for i in range(len(tmax)):
            axs.axvline(tmax[i], color='red', linestyle='--', linewidth=2)
    axs.set_title("Waveform")
    axs.set_xlabel("Time [μs]")
    axs.set_ylabel("Amplitude")
    if log:
        axs.set_yscale('log')
    axs.grid(True)
    axs.legend()
    
    
def plot_waveform_pmts(wvflist, n_time_bins, tbin=25e-3): 
    """
    Plot PMT waveforms
    
    """
    
    time = np.linspace(0, n_time_bins * tbin, n_time_bins)

    rows, cols = 10, 6  # Adjust rows and columns for clarity
    fig, axes = plt.subplots(rows, cols, figsize=(20, 15))  # Increase figure size for clarity

    #fig, axes = plt.subplots(20, 3, figsize=figsize)  # 10 rows, 6 columns
    axs = axes.flatten()

    # Plot each waveform
    for i, wvf in enumerate(wvflist):
        axs[i].plot(time, wvf, label=f"pmt = {i}", color="blue")
        axs[i].tick_params(labelsize=8)
        axs[i].set_xlabel("Time [μs]")
        axs[i].set_ylabel("Amplitude")
        axs[i].grid(True)
        axs[i].legend()

    plt.subplots_adjust(wspace=0.4, hspace=0.6)  # Increase spacing between subplots

    ###plt.tight_layout()
    


def plot_waveform_zoom_peak(wvf, tpeak, twindow, label="Waveform",
                            log=False, figsize=(18, 6), tbin=25e-3, tscale=False): 
    """
    Plot a waveform between tmin and tamx.  Time in ns
    """
    tmin = tpeak-twindow
    tmax = tpeak+twindow
    time = np.linspace(tmin, tmax, num=2*twindow)
    #print(tmin, tmax, time.shape)
    #print(wvf.shape, wvf[tmin:tmax].shape)

    fig, axs = plt.subplots(1, 1, figsize=figsize)
    if tscale:
        axs.plot(time * tbin, wvf[tmin:tmax], label=label, color="blue")
    else:
        axs.plot(time, wvf[tmin:tmax], label=label, color="blue")
    axs.set_title("Waveform")
    axs.set_xlabel("Time [μs]")
    axs.set_ylabel("Amplitude")
    if log:
        axs.set_yscale('log')
    axs.grid(True)
    axs.legend()
    plt.show()


def plot_waveform_zoom_peaks(wvf, tpeak, twindow, label="Waveform",
                            log=False, figsize=(18, 6), tbin=25e-3, tscale=False): 
    """
    Plot a waveform between tmin and tamx.  Time in ns
    """
    
    fig, axs = plt.subplots(1, len(tpeak), figsize=figsize)
    
    for i in range(len(tpeak)):
        tmin = tpeak[i]-twindow[i]
        tmax = tpeak[i] + twindow[i]
        time = np.linspace(tmin, tmax, num=2*twindow[i])
        #print(tmin, tmax, time.shape)
        #print(wvf.shape, wvf[tmin:tmax].shape)

        if tscale:
            axs[i].plot(time * tbin, wvf[tmin:tmax], label=label, color="blue")
        else:
            axs[i].plot(time, wvf[tmin:tmax], label=label, color="blue")
        axs[i].set_title("Waveform")
        axs[i].set_xlabel("Time [μs]")
        axs[i].set_ylabel("Amplitude")
    
        axs[i].grid(True)
        axs[i].legend()
    
    plt.tight_layout()
    plt.show()

    
    
def plot_waveform_right_peak(wvf, tpeak, twindow, label="Waveform",
                       log=False, figsize=(18, 6), tbin=25e-3): 
    """
    Plot a waveform between tmin and tamx.  Time in ns
    """
    tmin = tpeak
    tmax = tpeak+twindow
    time = np.linspace(tmin, tmax, num=twindow)
    #print(tmin, tmax, time.shape)
    #print(wvf.shape, wvf[tmin:tmax].shape)

    fig, axs = plt.subplots(1, 1, figsize=figsize)
    axs.plot(time, wvf[tmin:tmax], label=label, color="blue")
    axs.set_title("Waveform")
    axs.set_xlabel("Time [μs]")
    axs.set_ylabel("Amplitude")
    if log:
        axs.set_yscale('log')
    axs.grid(True)
    axs.legend()
    plt.show()


def plot_waveform_right_peaks(wvf, tpeak, twindow, label="Waveform",
                              log=False, figsize=(18, 6), tbin=25e-3, tscale=False): 
    """
    Plot a waveform between tmin and tamx.  Time in ns
    """

    fig, axs = plt.subplots(1, len(tpeak), figsize=figsize)
    
    for i in range(len(tpeak)):
        tmin = tpeak[i]
        tmax = tpeak[i] + twindow[i]
        time = np.linspace(tmin, tmax, num=twindow[i])
        
        if tscale:
            axs[i].plot(time * tbin, wvf[tmin:tmax], label=label, color="blue")
        else:
            axs[i].plot(time, wvf[tmin:tmax], label=label, color="blue")


        axs[i].set_title("Waveform")
        axs[i].set_xlabel("Time [μs]")
        axs[i].set_ylabel("Amplitude")
    
        axs[i].grid(True)
        axs[i].legend()
    
    plt.tight_layout()
    plt.show()

            
    
def plot_pmt_max_rms(X,Y,qpmt, std, dpmts):
    """
    rms and max signal for PMTs
    """
    
    stdmax = np.max(std)
    qmax = np.max(qpmt)

    for i in range(len(qpmt)):
        if i in dpmts:
            std[i] = 0.
            qpmt[i] = 0.
            
    fig, axs = plt.subplots(1, 2, figsize=(14, 6),dpi=160)

    # Plot 1: STD values
    scatter1 = axs[0].scatter(X, Y, c=std, s=100, cmap='plasma', alpha=0.8, edgecolors='k')
    axs[0].set_title("Baseline RMS vs PMT Positions")
    axs[0].set_xlabel("X [mm]")
    axs[0].set_ylabel("Y [mm]")
    plt.colorbar(scatter1, ax=axs[0], label="stds")

    # Plot 2: Max values
    scatter2 = axs[1].scatter(X, Y, c=qpmt, s=100, cmap='plasma', alpha=0.8, edgecolors='k')
    axs[1].set_title("Max Amp vs PMT Positions")
    axs[1].set_xlabel("X [mm]")
    axs[1].set_ylabel("Y [mm]")
    plt.colorbar(scatter2, ax=axs[1], label="Max Value")

    # Adjust layout and display
    plt.tight_layout()
    plt.show()


def histo_q_sigma(qpmt, stds):
    fig, axs = plt.subplots(1, 2, figsize=(14, 6),dpi=160)
    axs[0].hist(qpmt, bins=20, color='blue', edgecolor='black', alpha=0.7)

    # Add labels and title
    axs[0].set_xlabel('aDC counts')
    axs[0].set_ylabel('Frequency')
    axs[0].set_title('Sensor Q')

    axs[1].hist(stds, bins=20, color='blue', edgecolor='black', alpha=0.7)

    # Add labels and title
    axs[1].set_xlabel('aDC counts')
    axs[1].set_ylabel('Frequency')
    axs[1].set_title('Sensor STD')

    plt.show()



def histo_q(qpmt):
    fig, axs = plt.subplots(1, 1, figsize=(14, 6),dpi=160)
    axs.hist(qpmt, bins=20, color='blue', edgecolor='black', alpha=0.7)

    # Add labels and title
    axs.set_xlabel('aDC counts')
    axs.set_ylabel('Frequency')
    axs.set_title('Sensor Q')


    plt.show()


def plot_sipm_max_rms(X,Y,qmax, STD, units="adc"):
    fig, axs = plt.subplots(1, 2, figsize=(14, 6),dpi=100)

    # Plot 1: STD values
    scatter1 = axs[0].scatter(
            X, Y, 
            c=STD, 
            s=10, 
            cmap='plasma', 
            alpha=0.8, 
            edgecolors='k'
        )
    axs[0].set_title("Baseline RMS vs SiPM Positions")
    axs[0].set_xlabel("X [mm]")
    axs[0].set_ylabel("Y [mm]")
    plt.colorbar(scatter1, ax=axs[0], label="STD")

    # Plot 2: Max values
    if units == "adc":
        scale = 10
    else:
        scale = 1
    scatter2 = axs[1].scatter(X, Y, c=qmax, s=np.array(qmax)/scale,
                              cmap='plasma', alpha=0.8, edgecolors='k')
    axs[1].set_title("Max Amp vs SiPM Positions")
    axs[1].set_xlabel("X [mm]")
    axs[1].set_ylabel("Y [mm]")
    plt.colorbar(scatter2, ax=axs[1], label="Max Value")

    # Adjust layout and display
    plt.tight_layout()
    plt.show()




def plot_sum_PMT(df, log=False, offset=500):
    
    average_times = df.groupby('peak')['time'].mean()

    fig, axs = plt.subplots(1, 1, figsize=(14, 6),dpi=160)
    
  #  axs.plot(df['time'].values / 1000, df['ene'].values,
  #           marker='.', linestyle='', color='b', label='Energy vs Time')
    axs.plot(df['time'].values / 1000, df['ene'].values, label='Energy vs Time')

    # Add average time markers
    ###axs.scatter(average_times / 1000, [df['ene'].max() + offset] * len(average_times), 
    #                color='red', marker='1', label='Avg Time per Peak', zorder=5)

    # Add labels, title, and grid
    axs.set_xlabel('Time ($\mu$s)')
    axs.set_ylabel('Energy')
    axs.set_title('Energy vs Time ')
    axs.grid(True)
    if(log):
        axs.yscale('log')
    axs.legend()


    plt.show()



def plot_single_PMT(dfall, dfPMT, X, Y, npmt, log=False):


    # Extract time and energy data
    time = dfall['time']
    ene = dfPMT[(dfPMT['npmt'] == npmt)]['ene']

    # Create a figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(16, 4),
                             dpi=160, gridspec_kw={'width_ratios': [1, 3]})

    # Left subplot: X vs Y with a yellow dot
    ax1 = axes[0]
    ax1.scatter(X, Y, color='gray', alpha=0.7, label='PMT positions')
    ax1.scatter(X[npmt], Y[npmt], color='yellow', s=100, label=f'PMT {npmt}')
    for i in range(len(X)):
        ax1.annotate(str(i), (X[i], Y[i]), textcoords="offset points", xytext=(8, 0),
                     ha='center', fontsize=8, color='red')
        ax1.set_xlabel('X Position')
        ax1.set_ylabel('Y Position')
        ax1.set_title('PMT Map')

    # Right subplot: Energy vs Time
    ax2 = axes[1]
    ax2.plot(time, ene, marker='.', linestyle='', color='b', label='Energy vs Time')
    ax2.set_xlabel('Time ($\mu$s)')
    ax2.set_ylabel('Energy')
    ax2.set_title('Energy vs Time')
    if log:
        ax2.set_yscale('log')

    # Adjust layout and show the plots
    plt.tight_layout()
    plt.show()


def plot_peaks(df, tpeak, twindow, label="Waveform",
               log=False, figsize=(18, 6), tbin=1000): 
    """
    Plot a waveform between tmin and tamx.  Time in mus (divide by 1000)
    """
    tmin = tpeak-twindow
    tmax = tpeak+ 2*twindow

    print(tmin, tmax)
    ene =  df['ene'].values
    time = df['time'].values
    ptime = time[tmin:tmax]/tbin
    pene = ene[tmin:tmax]

    fig, axs = plt.subplots(1, 1, figsize=figsize)
    axs.plot(ptime, pene, label=label, color="blue")
    axs.set_title("S1/S2")
    axs.set_xlabel("Time [μs]")
    axs.set_ylabel("Amplitude")
    if log:
        axs.set_yscale('log')
    axs.grid(True)
    axs.legend()
