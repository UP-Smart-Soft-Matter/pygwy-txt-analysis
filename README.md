# README

## Overview
This script provides tools for analyzing and visualizing surface profiles exported from **Gwyddion** `.txt` files. It allows you to create height maps, profile plots, perform peak/valley detection, and compute statistical metrics and model fits from surface data.

---

## Classes and Initialization

### `PygwyTxt`
Handles reading, analysis, and visualization of a single Gwyddion `.txt` scan file.

**Initialization parameters:**
- `file_path: str` – Path to the input `.txt` file containing surface data.  
- `scan_size_x: float` – Horizontal scan size in micrometers.  
- `scan_size_y: float` – Vertical scan size in micrometers.  
- `name: str, optional` – Custom name for the dataset. Defaults to the filename.  
- `peak_finder_settings: PeakFinderSettings, optional` – Settings controlling peak/valley detection.

**Key methods:**
- `plot_scan()` – Creates a heatmap of the full scan.  
- `plot_profile()` – Plots the height profile along the central scan line.  
- `plot_profile_section(start, stop, line)` – Plots a selected section of a chosen scan line.  
- `plot_debug(line)` – Visualizes detected peaks and valleys for inspection.  
- `export_stats()` – Saves calculated statistics as a JSON file.  

**Automatically computed values:**
- Mean height and mean period  
- Standard deviation, minimum, maximum  

---

### `StatJson`
Collects and processes multiple JSON statistic files to visualize aggregated results or fits.

**Initialization parameters:**
- `base_path: str` – Directory containing JSON statistic files.

**Key methods:**
- `plot(plot_type, x_label, x_unit, plot_name_appendix='', model=None, params=None, show_title=True)` – Plots mean and standard deviation of height (`plot_type=0`) or period (`plot_type=1`). Optionally applies a model fit.  
- `export_plot_data(plot_type)` – Exports the plotted data (including fits) to a CSV file.  

---

### `PeakFinderSettings`
Defines configurable parameters for the `scipy.signal.find_peaks` algorithm to control peak and valley detection.

**Initialization parameters (all optional):**
- `height` – Required height of peaks.  
- `threshold` – Required vertical difference between peaks and neighbors.  
- `distance` – Minimum horizontal distance between peaks.  
- `prominence` – Required prominence of peaks.  
- `width` – Required width of peaks.  
- `wlen` – Window length for peak prominence evaluation.  
- `rel_height` – Relative height at which the peak width is measured.  
- `plateau_size` – Range of flat peak plateaus.  
