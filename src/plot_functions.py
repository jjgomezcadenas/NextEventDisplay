import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import numpy as np

def plot_adc_to_pes(axs):
    adc_to_pes = get_adc_to_pes()
    print(adc_to_pes)
    axs.hist(adc_to_pes, bins=20, color='blue', edgecolor='black', alpha=0.7)
    axs.set_xlabel('ADC to PES')
    axs.set_ylabel('Frequency')
    axs.set_title('ADC to PES for NEXT-100')


def plot_sum_waveform(siwf, cwf, n_timesi_bins, n_time_bins, run_number, event_number,
                      s1tmx = 2600, s1tmn=0, s2tmx = 2600, s2tmn=0,
                      figsize=(18, 6), tbin=25e-3): 
    """
    Plots the sum of the waveforms (calibrated) for SiPMs and PMTs. 
    Time in mus (thus multiplity the time of PMTs by 25 10^3, or 25 ns)
    Interface for notebook

    """

    fig, axs = plt.subplots(1, 2, figsize=figsize)

    plot_sum_waveform_tk(axs, siwf, cwf, n_timesi_bins, n_time_bins, run_number, event_number,
                      s1tmx, s1tmn, s2tmx, s2tmn, tbin)
    
    fig.tight_layout()


def plot_sum_waveform_tk(axs, siwf, cwf, n_timesi_bins, n_time_bins, run_number, event_number,
                         s1tmx, s1tmn, s2tmx, s2tmn, tbin):
    """
    Plots the sum of the waveforms (calibrated) for SiPMs and PMTs. 
    Time in mus (thus multiplity the time of PMTs by 25 10^3, or 25 ns)

    """
    
    
    label=f"evt={event_number} run={run_number}"
    timesi = np.linspace(0, n_timesi_bins, n_timesi_bins)
    time = np.linspace(0, n_time_bins * tbin, n_time_bins) 

    axs[0].plot(timesi, siwf, label=label, color="blue")
    axs[0].set_title("Calibrated Sum SiPM")
    axs[0].set_xlabel("Time [μs]")
    axs[0].set_ylabel("Amplitude")
    axs[0].axvline(s1tmn, color='yellow', linestyle='--', linewidth=2)
    axs[0].axvline(s1tmx, color='yellow', linestyle='--', linewidth=2)
    axs[0].axvline(s2tmn, color='green', linestyle='--', linewidth=2)
    axs[0].axvline(s2tmx, color='green', linestyle='--', linewidth=2)
    #axs[0].set_yscale('log')
    axs[0].grid(True)
    axs[0].legend()

    axs[1].plot(time, cwf, label=label, color="red")
    axs[1].axvline(s1tmn, color='yellow', linestyle='--', linewidth=2)
    axs[1].axvline(s1tmx, color='yellow', linestyle='--', linewidth=2)
    axs[1].axvline(s2tmn, color='green', linestyle='--', linewidth=2)
    axs[1].axvline(s2tmx, color='green', linestyle='--', linewidth=2)
    axs[1].set_title("Calibrated Sum PMT")
    axs[1].set_xlabel("Time [μs]")
    axs[1].set_ylabel("Amplitude")
    #axs[0].set_yscale('log')
    axs[1].grid(True)
    axs[1].legend()
    
    
    
def plot_waveform(wvf, n_time_bins, run_number, event_number, peaks=[], widths=[], left_ips=[], right_ips=[], 
                  figsize=(18, 6), tbin=25e-3): 
    
    """
    Plot a waveform. Optionaly mark the position of peaks
    
    """
    fig, axs = plt.subplots(1, 1, figsize=figsize)

    plot_waveform_tk(axs, wvf, n_time_bins, run_number, event_number, 
                     peaks, widths, left_ips, right_ips, tbin)
    fig.tight_layout()


def plot_waveform_tk(axs, wvf, n_time_bins, run_number, event_number, 
                     peaks, widths, left_ips, right_ips, tbin): 
    
    label=f"evt={event_number} run={run_number}"
    
    time = np.linspace(0, n_time_bins * tbin, n_time_bins)  
    axs.plot(time, wvf, label=label, color="blue")
    
    for i, pk in enumerate(peaks):
        lcut = tbin * (left_ips[i] -  1 * widths[i])
        rcut = tbin * (right_ips[i] +  2 * widths[i])
        print(f"left cut = {lcut}, right cut = {rcut}")
        axs.axvline(lcut, color='red', linestyle='--', linewidth=2)
        axs.axvline(rcut, color='red', linestyle='--', linewidth=2)

    axs.set_title("Waveform")
    axs.set_xlabel("Time [μs]")
    axs.set_ylabel("Amplitude")
   
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

    fig.tight_layout()



def plot_sipmw(SIPMW,figsize=(18, 6)): 
    """
    Plot the SIPMW
    
    """
    
    fig, axs = plt.subplots(1, len(SIPMW), figsize=figsize)  # Increase figure size for clarity
    
    if len(SIPMW) == 1:
        axs.plot(SIPMW[0])
        axs.set_xlabel(f"SiPM number")
        axs.set_ylabel("Energy in PES")
        axs.grid(True)
    else:

        for i, sipmw in enumerate(SIPMW):
            axs[i].plot(sipmw)
            axs[i].set_xlabel(f"SiPM number in window {i}")
            axs[i].set_ylabel("Energy in PES")
            axs[i].grid(True)
    
    fig.tight_layout()


def plot_signal_sipms(esi, run_number, event_number,figsize=(18, 6)): 
    """
    Plot the energy in the SiPMs
    
    """
    label=f"evt={event_number} run={run_number}"
    fig, axs = plt.subplots(1, 1, figsize=figsize)  # Increase figure size for clarity
    axs.plot(esi,label=label, color="red")
    axs.set_xlabel("SiPM number")
    axs.set_ylabel("Energy in PES")
    axs.grid(True)
    axs.legend()
    fig.tight_layout()

    
def plot_waveform_zoom_peak(wvf, tpeak, tleft, tright, twindow, run_number, event_number,
                            figsize=(18, 6), tbin=25e-3, tscale=False): 
    """
    Plot a waveform between tmin and tamx.  Time in ns
    """
    tmin = tpeak - twindow
    tmax = tpeak + twindow
    time = np.linspace(tmin, tmax, num=2*twindow)
    
    label=f"evt={event_number} run={run_number}"

    fig, axs = plt.subplots(1, 1, figsize=figsize)
    if tscale:
        axs.plot(time * tbin, wvf[tmin:tmax], label=label, color="blue")
        axs.axvline(tleft * tbin, color='red', linestyle='--', linewidth=2)
        axs.axvline(tright * tbin, color='red', linestyle='--', linewidth=2)
        axs.set_xlabel("Time [μs]")
    else:
        axs.plot(time, wvf[tmin:tmax], label=label, color="blue")
        axs.axvline(tleft, color='red', linestyle='--', linewidth=2)
        axs.axvline(tright, color='red', linestyle='--', linewidth=2)
        axs.set_xlabel("Time [samples]")
        
    axs.set_title("Waveform")
    
    axs.set_ylabel("Amplitude")
    
    axs.grid(True)
    axs.legend()
    plt.show()


def plot_waveform_zoom_peaks(wvf, run_number, event_number, tpeaks, tlefts, trights, twindows, 
                            figsize=(18, 6), tbin=25e-3, tscale=False): 
    """
    Plot a waveform between tmin and tamx.  Time in ns
    """
    
    fig, axs = plt.subplots(1, len(tpeaks), figsize=figsize)
    plot_waveform_zoom_peaks_tk(axs, wvf, run_number, event_number, tpeaks, tlefts, trights, twindows, 
                               tbin, tscale)
    fig.tight_layout()
   

def plot_waveform_zoom_peaks_tk(axs, wvf, run_number, event_number, tpeaks, tlefts, trights, twindows, 
                                tbin, tscale): 
    """
    Plot a waveform between tmin and tamx.  Time in ns
    """
    
    label=f"evt={event_number} run={run_number}"

    if len(tpeaks) == 1:
        tmin = tlefts[0] - int(twindows[0])
        tmax = trights[0] + int(twindows[0])
        time = np.linspace(tmin, tmax, num=tmax-tmin)
    
        if tscale:
            axs.plot(time * tbin, wvf[tmin:tmax], label=label, color="blue")
            axs.axvline(tlefts[0] * tbin, color='red', linestyle='--', linewidth=2)
            axs.axvline(trights[0] * tbin, color='red', linestyle='--', linewidth=2)
            axs.set_xlabel("Time [μs]")
            axs.set_ylabel("Amplitude")
            axs.grid(True)
            axs.legend()
        else:
            axs.plot(time, wvf[tmin:tmax], label=label, color="blue")
            axs.axvline(tlefts[0], color='red', linestyle='--', linewidth=2)
            axs.axvline(trights[0], color='red', linestyle='--', linewidth=2)
            axs.set_xlabel("Time [samples]")
            axs.set_ylabel("Amplitude")
            axs.grid(True)
            axs.legend()
    else:

        for i, tpeak in enumerate(tpeaks):
            #tmin = tpeak - twindows[i]
            #tmax = tpeak + twindows[i]
            tmin = tlefts[i] - int(twindows[i])
            tmax = trights[i] + int(twindows[i])
            time = np.linspace(tmin, tmax, num=tmax-tmin)
        
            if tscale:
                axs[i].plot(time * tbin, wvf[tmin:tmax], label=label, color="blue")
                axs[i].axvline(tlefts[i] * tbin, color='red', linestyle='--', linewidth=2)
                axs[i].axvline(trights[i] * tbin, color='red', linestyle='--', linewidth=2)
                axs[i].set_xlabel("Time [μs]")
                axs[i].set_ylabel("Amplitude")
                axs[i].grid(True)
                axs[i].legend()
            else:
                axs[i].plot(time, wvf[tmin:tmax], label=label, color="blue")
                axs[i].axvline(tlefts[i], color='red', linestyle='--', linewidth=2)
                axs[i].axvline(trights[i], color='red', linestyle='--', linewidth=2)
                axs[i].set_xlabel("Time [samples]")
                axs[i].set_ylabel("Amplitude")
                axs[i].grid(True)
                axs[i].legend()

    
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


def plot_sipm(X,Y,QSIPM, XG, YG, scale=1,figsize=(14, 6)):
    
    fig, axs = plt.subplots(1, len(QSIPM), figsize=figsize)

    if len(QSIPM) == 1:
        qmax = np.max(QSIPM[0])
        scatter = axs.scatter(X, Y, c=QSIPM[0], s=np.array(QSIPM[0])/scale,
                                cmap='plasma', alpha=0.8, edgecolors='k')
        scatter = axs.scatter(XG[0], YG[0], c=10*qmax, s=np.array(qmax)/(scale/3),
                                cmap='plasma', alpha=0.8, edgecolors='k')
        axs.set_title("Amp vs SiPM Positions")
        axs.set_xlabel("X [mm]")
        axs.set_ylabel("Y [mm]")
        plt.colorbar(scatter, ax=axs)
    else:
        for i, qsipm in enumerate(QSIPM):
            qmax = np.max(QSIPM[i])
            scatter = axs[i].scatter(X, Y, c=qsipm, s=np.array(qsipm)/scale,
                                cmap='plasma', alpha=0.8, edgecolors='k')
            scatter = axs[i].scatter(XG[i], YG[i], c=10*qmax, s=np.array(qmax)/(scale/3),
                                cmap='plasma', alpha=0.8, edgecolors='k')
            axs[i].set_title("Amp vs SiPM Positions")
            axs[i].set_xlabel("X [mm]")
            axs[i].set_ylabel("Y [mm]")
            plt.colorbar(scatter, ax=axs)


    # Adjust layout and display
    #fig.tight_layout()
    


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
