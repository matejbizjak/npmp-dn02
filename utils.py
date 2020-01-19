# simulation
import scipy.signal as scipy


# peaks and valleys detection
def find_peaks(sim_results):
    peaks = scipy.find_peaks(sim_results)[0]  # polozaji
    minimums = scipy.find_peaks(list(map(lambda x: x * (-1), sim_results)))[0]  # polozaji
    peaks_vals = list(map(lambda x: sim_results[x], peaks))  # vrednosti
    minimums_vals = list(map(lambda x: sim_results[x], minimums))  # vrednosti

    return peaks, minimums, peaks_vals, minimums_vals


# check if oscilating
def oscilating(peak_vals, tolerance):
    try:
        return abs(peak_vals[-2] - peak_vals[-3]) < tolerance
    except:
        return False


# find amplitude and periode
def eval_signal(peaks, min_vals, peak_vals):
    amplitude = (peak_vals[-2] - min_vals[-2]) / 2
    periode = (peaks[-2] - peaks[-3]) / 10  # 10 = n/t_end
    return amplitude, periode
