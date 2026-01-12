# README

## Overview

This script provides tools for analyzing and visualizing surface profiles exported from **Gwyddion** `.txt` files. It allows you to create height maps, profile plots, perform peak/valley detection, and compute statistical metrics and model fits from surface data.

---
## Installation

`pip install pygwy_txt_analysis`

---
### Example: `batch_analysis.py`

### Important File Naming Convention
To ensure proper ordering in the aggregated plots, the filenames of the `.txt` files should include a numeric designator at the start representing their position or parameter (e.g., time). This number is used to sort the files before processing. Example:


```
01_sample.txt
02_sample.txt
03_sample.txt
```


---


## How to Use
1. Select the folder containing the `.txt` scan files when prompted.
2. Set scan sizes (`scan_size_x`, `scan_size_y`) in micrometers.
3. Optionally configure peak detection using `PeakFinderSettings`. Example:


```python
peak_finder_settings = pygwy_txt_analysis.PeakFinderSettings(prominence=0.3e-7)
```


4. The script will:
- Load each `.txt` file in order.
- Create and display heatmaps (`plot_scan`) and central profiles (`plot_profile`).
- Detect peaks and valleys (`plot_debug`).
- Export statistics as JSON in an `export` folder.


5. Aggregate statistics are plotted using `StatJson` and exported to CSV for height and period:


```python
stats = pygwy_txt_analysis.StatJson(os.path.join(basepath, 'export'))
stats.plot(0, x_label, x_unit) # height
stats.plot(1, x_label, x_unit) # period
stats.export_plot_data(0)
stats.export_plot_data(1)
```
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
