# README

## Overview
This script provides tools for analyzing and visualizing surface profiles exported from **Gwyddion** `.txt` files. It allows you to create height maps, profile plots, perform peak/valley detection, and compute statistical metrics and model fits from surface data.

---

## Main Components

### `PygwyTxt`
Handles reading, analysis, and visualization of a single Gwyddion `.txt` scan file.

**Key methods:**
- **`plot_scan()`** – Creates a heatmap of the full scan.  
- **`plot_profile()`** – Plots the height profile along the central scan line.  
- **`plot_profile_section(start, stop, line)`** – Plots a selected section of a chosen scan line.  
- **`plot_debug(line)`** – Visualizes detected peaks and valleys for inspection.  
- **`export_stats()`** – Saves calculated statistics as a JSON file.  

**Automatically computed values:**
- Mean height and mean period  
- Standard deviation, minimum, maximum  

---

### `StatJson`
Collects and processes multiple JSON statistic files to visualize aggregated results or fits.

**Key methods:**
- **`plot(plot_type, x_label, x_unit, ...)`** – Plots mean and standard deviation of height (`plot_type=0`) or period (`plot_type=1`). Optionally applies a model fit.  
- **`export_plot_data(plot_type)`** – Exports the plotted data (including fits) to a CSV file.  

---

### `PeakFinderSettings`
Defines configurable parameters for the `scipy.signal.find_peaks` algorithm, such as:
- `height`, `threshold`, `distance`, `prominence`, `width`, `wlen`, `rel_height`, `plateau_size`

These parameters are used to fine-tune peak and valley detection during surface analysis.
