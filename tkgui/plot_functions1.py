import numpy as np
from ic_functions import get_adc_to_pes


def plot_adc_to_pes(axs):
    adc_to_pes = get_adc_to_pes()
    print(adc_to_pes)
    axs.hist(adc_to_pes, bins=20, color='blue', edgecolor='black', alpha=0.7)
    axs.set_xlabel('ADC to PES')
    axs.set_ylabel('Frequency')
    axs.set_title('ADC to PES for NEXT-100')


def plot_sum_waveform(axs, siwf, cwf, n_timesi_bins, n_time_bins, run_number, event_number,
                      tbin=25e-3): 
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
    #axs[0].set_yscale('log')
    axs[0].grid(True)
    axs[0].legend()
    axs[1].plot(time, cwf, label=label, color="red")
    axs[1].set_title("Calibrated Sum PMT")
    axs[1].set_xlabel("Time [μs]")
    axs[1].set_ylabel("Amplitude")
    #axs[0].set_yscale('log')
    axs[1].grid(True)
    axs[1].legend()
    

def plot_sipm_max_rms(axs, X,Y,qmax, STD, units="adc"):
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