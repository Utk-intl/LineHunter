import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
from scipy.signal import savgol_filter, find_peaks, peak_prominences, peak_widths, fftconvolve
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d, CubicSpline, UnivariateSpline

def load_and_process_spectra(obs_filename, syn_filename='m05t4000g3.0_0.spc.r00'):
    """
    Load observed and synthetic spectra files and perform continuum normalization
    
    Parameters:
    obs_filename: str - filename of observed spectra
    syn_filename: str - filename of synthetic spectra (default provided)
    
    Returns:
    tuple: (wl_obs, flux_obs, normalized_flux, continuum, spline_fit, wl_syn, flux_syn)
    """
    # 1. Load the observed spectra file
    data_obs = np.loadtxt(obs_filename)
    wl_obs = np.array(data_obs[:,0])
    flux_obs = np.array(data_obs[:,1])

    #2. load synthetic spectra file 
    data_syn = np.loadtxt(syn_filename)
    flux_syn = np.array(data_syn[:,1])
    wl_syn = np.array(data_syn[:,0])

    # Maximum and minimum value of wavelength
    wl_max = np.max(wl_obs)
    wl_min = np.min(wl_obs)

    # Bin width
    bin_width = 5
    bin_edges = np.arange(wl_min, wl_max+bin_width, bin_width)

    # Defining a dataframe and binning
    df = pd.DataFrame({'wl': wl_obs, 'flux': flux_obs})
    df['wl_bins'] = pd.cut(df['wl'], bins=bin_edges)

    anchor_wl = []
    anchor_flux = []

    # Group by wavelength bins
    grouped = df.groupby('wl_bins',observed='bool')

    for bin_range, group in grouped:
        bin_flux = group['flux'].values
        bin_wl = group['wl'].values  
        
        if len(bin_flux) == 0:
            continue
        
        min_noise = np.inf
        best_wl, best_flux = None, None
        for p in np.arange(96, 99.9, 1.0):
            
            clipped_data = sigmaclip(bin_flux, low=2.0, high=3.0)[0]
            
            
            if len(clipped_data) == 0:
                continue
            
            target_flux = np.percentile(clipped_data, p)
            idx = np.argmin(np.abs(bin_flux - target_flux))
            try_wl = bin_wl[idx]  
            try_flux = bin_flux[idx]  

            if try_flux < 0.99 :
                continue 
            # Calculate local normalization
            local_norm = bin_flux/try_flux
            smoothed = savgol_filter(local_norm, min(11, len(local_norm) - (len(local_norm) % 2) - 1), 3)
            noise = np.std(local_norm - smoothed)
            
            if noise < min_noise:
                min_noise = noise  
                best_wl, best_flux = try_wl, try_flux
        
        if best_wl is not None:
            anchor_wl.append(best_wl)
            anchor_flux.append(best_flux)

    if len(anchor_wl) >= 3:  
        sorted_indices = np.argsort(anchor_wl)
        anchor_wl = np.array(anchor_wl)[sorted_indices]
        anchor_flux = np.array(anchor_flux)[sorted_indices]
        s = len(anchor_wl)*1e-6
        spline_fit = UnivariateSpline( anchor_wl,anchor_flux ,s=s)
        continuum = spline_fit(wl_obs)
        normalized_flux = flux_obs/continuum
    else:
        raise ValueError("Not enough anchor points found for continuum fitting")
    
    return wl_obs, flux_obs, normalized_flux, continuum, spline_fit, wl_syn, flux_syn

def plot_continuum_normalization(wl_obs, flux_obs, continuum, normalized_flux):
    """Plot the continuum normalization results"""
    plt.figure(figsize=(12, 6))
    plt.plot(wl_obs, flux_obs, label='Original Flux', alpha=0.6)
    plt.plot(wl_obs, continuum, label='Fitted Continuum', linestyle='--')
    plt.plot(wl_obs, normalized_flux, label='Normalized Flux', linewidth=1.2)
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Flux")
    plt.legend()
    plt.title("Continuum Normalization")
    plt.grid(True)
    plt.show()

# define an inverted gaussian function
def gaussian(x, a, c, x0, sigma):
    '''Here a=amplitude, x0=centre of gaussian, c=continuum level, sigma=standard deviation'''
    return c - a*np.exp(-(x-x0)**2 / (2*sigma**2))  # Fixed missing parentheses

#----------------------------------------------------------------------------------------------------------------
#define a function to analyze the section of the spectrum with quality check 
def analyze_spectrum(wl_section, flux_section, section_num, spline_fit, wl_obs):
    #data smoothing 
    valid_peaks = []
    try:
        smoothed_flux = savgol_filter(flux_section, window_length=11, polyorder=3)
    except:
        smoothing_window = 5
        smoothed_flux = fftconvolve(flux_section, np.ones(smoothing_window)/smoothing_window, mode='same')
    
    #inverting the flux to assume the absorption/emission lines as peaks
    flux_range = np.max(smoothed_flux) - np.min(smoothed_flux)
    if flux_range < 0.01:  # Added a small threshold instead of just < 0
        print(f'Section {section_num} is flat, skipped!!')
        return None  # Fixed lowercase 'none' to 'None'
    
    # Use flux_section instead of normalized_flux (which isn't defined in this scope)
    inv_flux = 1 - flux_section
    
    #rough calculation of parameters for peak detection
    #initial rough peaks
    noise_level = np.std(inv_flux)
    clip_noise_frac = 3
    prom_rough = noise_level * clip_noise_frac
    
    # Fixed variable name consistency
    peaks_r, _ = find_peaks(inv_flux, prominence=prom_rough)
    if len(peaks_r) < 2:
        prom_rough = noise_level * (clip_noise_frac/2)  
        peaks_r, _ = find_peaks(inv_flux, prominence=prom_rough)
    
    #distance parameter   
    if len(peaks_r) > 1:
        mean_sp = np.mean(np.diff(wl_section[peaks_r]))  
        dx = np.median(np.diff(wl_section))  
        dist = max(1, int(0.5 * (mean_sp/dx)))
    else:
        dist = 5  # Default value if not enough peaks
        
    #width parameter
    if len(peaks_r) > 0:
        width_r = peak_widths(inv_flux, peaks_r, rel_height=0.5)[0]
        min_w, max_w = np.percentile(width_r, [10, 90])
    else:
        min_w, max_w = 3, 15
    
    #prominence parameter
    candidate_peaks, _ = find_peaks(inv_flux, distance=dist, width=(min_w, max_w))
    
    # Use the peak_prominences function to calculate prominences
    if len(candidate_peaks) > 0:
        prominences = peak_prominences(inv_flux, candidate_peaks)[0]
        
        # Get peaks with sufficient prominence
        final_peaks = []
        for i, peak in enumerate(candidate_peaks):
            if prominences[i] >= noise_level*clip_noise_frac:
                final_peaks.append(peak)
    else:
        final_peaks = []
    
    peak_indices = final_peaks
    
    if len(peak_indices) == 0:
        print(f'No significant absorption line found in section {section_num}')
    
    #verify that each candidate peak is a "real" absorption line of sufficient depth
    min_depth = 0.05 
    for i in peak_indices:
        rel_depth = inv_flux[i]
        if rel_depth >= min_depth: 
            valid_peaks.append((i, rel_depth))
    
    if not valid_peaks:
        
        return None
        
    # Sort valid_peaks by depth (descending)
    valid_peaks.sort(key=lambda x: x[1], reverse=True)  
    
    results = []
    
    avg_steps = 0.01  # Wavelength step size
    continuum = np.mean(spline_fit(wl_obs))  # Assuming normalized flux has continuum level at 1.0
    
    fig, axs = plt.subplots(1, min(3, len(valid_peaks)), figsize=(15, 5), constrained_layout=True)
    if len(valid_peaks) == 1:
        axs = [axs]
    elif len(valid_peaks) == 0:
        return None
        
    for i in range(min(3, len(valid_peaks))):
        idx, depth = valid_peaks[i]
        window = 50
        left_idx = max(0, idx-window)
        right_idx = min(len(wl_section), idx+window)
        line_wl = wl_section[left_idx:right_idx]
        line_flux = flux_section[left_idx:right_idx]
        
        if len(line_wl) < 10:  # Fixed condition (was: line_wl < 10)
            print(f'Not enough points for fitting peak at {wl_section[idx]}')
            continue
            
        centre_wave = wl_section[idx]
        amplitude = continuum - np.min(line_flux)
        half_depth = line_flux[idx-left_idx] + amplitude/2
        
        left_half_idx = right_half_idx = idx-left_idx
        
        for j in range(idx-left_idx, 0, -1):
            if line_flux[j] >= half_depth:
                left_half_idx = j
                break
                
        for k in range(idx-left_idx, len(line_flux)-1):
            if line_flux[k] >= half_depth:  
                right_half_idx = k
                break
                
        measured_fwhm = max(line_wl[right_half_idx] - line_wl[left_half_idx], avg_steps*3)
        sigma_guess = measured_fwhm / 2.3548
        
        p0 = [amplitude, continuum, centre_wave, sigma_guess]
        
        try:
            lower_bounds = [0, 0, line_wl[0], avg_steps]
            upper_bounds = [amplitude*2, continuum*2, line_wl[-1], (line_wl[-1] - line_wl[0])/2]
            
            popt, pcov = curve_fit(gaussian, line_wl, line_flux, p0=p0, bounds=(lower_bounds, upper_bounds))
            
            # Unpack fitted parameters properly
            amp_fit, cont_fit, centre_fit, sigma_fit = popt
            
            fwhm = sigma_fit * 2.3548
            resolving_power = centre_fit / fwhm  
            
            residuals = line_flux - gaussian(line_wl, *popt)
            ss_residuals = np.sum(residuals**2)
            ss_total = np.sum((line_flux - np.mean(line_flux))**2) 
            r_squared = 1 - (ss_residuals / ss_total if ss_total > 0 else 0)
            
            results.append({
                'center': centre_fit, 
                'sigma': sigma_fit, 
                'fwhm': fwhm,
                'resolving_power': resolving_power,  
                'amplitude': amp_fit,  
                'continuum': cont_fit,  
                'r_squared': r_squared, 
                'line_wl': line_wl, 
                'line_flux': line_flux
            })
            
            # Handle plotting
            if i < len(axs):
                ax = axs[i]
                ax.plot(line_wl, line_flux, 'o', markersize=3, label='data')
                
                x_fit = np.linspace(line_wl[0], line_wl[-1], 1000)
                y_fit = gaussian(x_fit, *popt)
                ax.plot(x_fit, y_fit, 'r-', label='fit')
                
                half_depth = cont_fit - amp_fit/2
                half_width = sigma_fit * np.sqrt(2 * np.log(2))
                left_half = centre_fit - half_width
                right_half = centre_fit + half_width
                
                # Fixed string interpolation
                ax.axvline(x=centre_fit, color='green', linestyle='--')
                ax.axhline(y=half_depth, color='blue', linestyle='--')  # Fixed: x= -> y=
                ax.axvline(x=left_half, color='purple', linestyle=':')
                ax.axvline(x=right_half, color='red', linestyle=':')
               
                # Add resolving power to title
                ax.set_title(f"λ={centre_fit:.2f}Å, R={resolving_power:.0f}")
                ax.set_xlabel('Wavelength (Å)')
                ax.set_ylabel('Flux')
                ax.grid(alpha=0.3)
                ax.legend(fontsize='small')
                
        except Exception as e:
            print(f"Fitting failed for peak at {wl_section[idx]:.2f} Å: {e}")
            
            # Plot the raw data if fitting failed
            if i < len(axs):
                ax = axs[i]
                ax.plot(line_wl, line_flux, 'o', markersize=3)
                ax.set_title(f"Fitting failed at λ={wl_section[idx]:.2f}Å")
                ax.set_xlabel('Wavelength (Å)')
                ax.set_ylabel('Flux')
                ax.grid(alpha=0.3)
    
    plt.suptitle(f"Section {section_num}: {wl_section[0]:.2f} - {wl_section[-1]:.2f} Å", fontsize=16)
    plt.tight_layout()
    plt.savefig(f'section_{section_num}_peaks.png', dpi=200)
    plt.close()
    
    return results

#--------------------------------------------------------------------------------------------------------------------------------

# Main code to run the analysis (assuming wl_obs and normalized_flux are defined from previous code)
def run_resolution_analysis(wl_obs, normalized_flux, spline_fit):
    num_sections = 10
    section_size = len(wl_obs) // num_sections
    
    all_results = []
    best_results_per_section = []
    
    for i in range(num_sections):
        start = i*section_size 
        end = (i+1)*section_size if i < num_sections-1 else len(wl_obs)  # Fixed: num_section -> num_sections
        wl_section = wl_obs[start:end]
        flux_section = normalized_flux[start:end]
        
        
        
        # Fixed function name: analyze_section -> analyze_spectrum
        results = analyze_spectrum(wl_section, flux_section, i+1, spline_fit, wl_obs) 
        if results:
            all_results.extend(results)
            # Fixed key name: res_power -> resolving_power
            best_result = max(results, key=lambda x: x['resolving_power'])
            best_results_per_section.append(best_result)
    
    # Create a summary plot
    if all_results:
        plt.figure(figsize=(15,8))
        plt.plot(wl_obs, normalized_flux, 'b-', alpha=0.5)  # Changed flux_obs to normalized_flux for consistency
        
        # Mark all the identified absorption lines
        for result in all_results:
            plt.axvline(x=result['center'], color='r', alpha=0.3, linestyle=':')
        
        # Mark the best absorption lines more prominently
        for result in best_results_per_section:
            plt.axvline(x=result['center'], color='g', alpha=0.7, linestyle='-', linewidth=1.5)
            plt.text(result['center'], 0.9, f"R={result['resolving_power']:.0f}", 
                     rotation=90, verticalalignment='top', fontsize=8)
        
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Flux')
        plt.title('Full Spectrum with Identified Absorption Lines')
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig('full_spectrum_analysis.png', dpi=300)  # Reduced dpi for faster processing
        
        # Calculate average resolving power
        resolving_powers = [r['resolving_power'] for r in all_results]
        avg_resolving_power = np.mean(resolving_powers)
        median_resolving_power = np.median(resolving_powers)
        print(f"\nResults:")
        print(f"Average Resolving Power: {avg_resolving_power:.0f}")
        print(f"Median Resolving Power: {median_resolving_power:.0f}")
        
        return avg_resolving_power, all_results, best_results_per_section
    else:
        print("No valid absorption lines found for resolving power calculation.")
        return None, [], []

def plot_best_absorption_lines(all_results, gaussian):
    """Plot the best absorption lines with Gaussian fits"""
    if all_results:
        # sort and get best two
        all_results.sort(key=lambda x: x['resolving_power'], reverse=True)
        best_result = all_results[0]
        second_best_result = all_results[1] if len(all_results) > 1 else None

        fig, axs = plt.subplots(1, 2 if second_best_result else 1, figsize=(15, 6))
        if second_best_result:
            ax1, ax2 = axs
        else:
            ax1 = axs

        # plot best
        wl1 = best_result['line_wl']
        fl1 = best_result['line_flux']
        ax1.plot(wl1, fl1, 'o', markersize=3, label='Data')
        x1 = np.linspace(wl1[0], wl1[-1], 1000)
        y1 = gaussian(x1, best_result['amplitude'], best_result['continuum'],
                      best_result['center'], best_result['sigma'])
        ax1.plot(x1, y1, 'r-', linewidth=2, label='Gaussian Fit')
        hd1 = best_result['continuum'] - best_result['amplitude']/2
        hw1 = best_result['sigma'] * np.sqrt(2 * np.log(2))
        lh1 = best_result['center'] - hw1
        rh1 = best_result['center'] + hw1
        ax1.axvline(x=best_result['center'], linestyle='--')
        ax1.axhline(y=hd1, linestyle='--')
        ax1.axvline(x=lh1, linestyle=':')
        ax1.axvline(x=rh1, linestyle=':')
        ax1.text(0.05, 0.05,
                 f"R: {best_result['resolving_power']:.0f}\nR²: {best_result['r_squared']:.3f}",
                 transform=ax1.transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.8))
        ax1.set_xlabel('Wavelength (Å)')
        ax1.set_ylabel('Flux')
        ax1.set_title('Best Absorption Line')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # plot second best if exists
        if second_best_result:
            wl2 = second_best_result['line_wl']
            fl2 = second_best_result['line_flux']
            ax2.plot(wl2, fl2, 'o', markersize=3, label='Data')
            x2 = np.linspace(wl2[0], wl2[-1], 1000)
            y2 = gaussian(x2, second_best_result['amplitude'],
                          second_best_result['continuum'], second_best_result['center'],
                          second_best_result['sigma'])
            ax2.plot(x2, y2, 'r-', linewidth=2, label='Gaussian Fit')
            hd2 = second_best_result['continuum'] - second_best_result['amplitude']/2
            hw2 = second_best_result['sigma'] * np.sqrt(2 * np.log(2))
            lh2 = second_best_result['center'] - hw2
            rh2 = second_best_result['center'] + hw2
            ax2.axvline(x=second_best_result['center'], linestyle='--')
            ax2.axhline(y=hd2, linestyle='--')
            ax2.axvline(x=lh2, linestyle=':')
            ax2.axvline(x=rh2, linestyle=':')
            ax2.text(0.05, 0.05,
                     f"R: {second_best_result['resolving_power']:.0f}\nR²: {second_best_result['r_squared']:.3f}",
                     transform=ax2.transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.8))
            ax2.set_xlabel('Wavelength (Å)')
            ax2.set_ylabel('Flux')
            ax2.set_title('Second Best Absorption Line')
            ax2.legend()
            ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('best_absorption_lines.png', dpi=300)

        # summary of high-quality fits
        valid_results = [r for r in all_results if r['r_squared'] > 0.9]
        if valid_results:
            best_valid = valid_results[0]
            print(f"\nBest high-quality resolving power: {best_valid['resolving_power']:.0f}")
            print(f"At wavelength: {best_valid['center']:.6f} Å with FWHM: {best_valid['fwhm']:.6f} Å")
            print(f"Fit quality (R²): {best_valid['r_squared']:.3f}")
        if len(valid_results) > 1:
            avg_r = np.mean([r['resolving_power'] for r in valid_results])
            std_r = np.std([r['resolving_power'] for r in valid_results])
            print(f"\nAverage resolving power from {len(valid_results)} validated lines: {avg_r:.0f} ± {std_r:.0f}")
    else:
        print("No suitable absorption lines found for resolving power calculation.")

# Example usage function
def complete_spectrum_analysis(obs_filename, syn_filename='m05t4000g3.0_0.spc.r00'):
    """
    Complete spectrum analysis workflow
    
    Parameters:
    obs_filename: str - filename of observed spectra
    syn_filename: str - filename of synthetic spectra (optional)
    
    Returns:
    tuple: (avg_resolving_power, all_results, best_results_per_section)
    """
    # Load and process spectra
    wl_obs, flux_obs, normalized_flux, continuum, spline_fit, wl_syn, flux_syn = load_and_process_spectra(obs_filename, syn_filename)
    
    # Plot continuum normalization
    plot_continuum_normalization(wl_obs, flux_obs, continuum, normalized_flux)
    
    # Run resolution analysis
    avg_R, all_results, best_results = run_resolution_analysis(wl_obs, normalized_flux, spline_fit)
    
    # Plot best absorption lines
    if all_results:
        plot_best_absorption_lines(all_results, gaussian)
    
    return avg_R, all_results, best_results