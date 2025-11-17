import numbers
from collections.abc import Sequence
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from astropy import units as u
import math
import os
from scipy.signal import find_peaks
import json
import tkinter as tk
from tkinter import filedialog
import glob2
import re
import csv
from lmfit import Model

def homogenize_array(array):
    """
    Pads each row of a 2D array with NaN values so all rows have equal length.

    Parameters
    ----------
    array : list of lists or list of np.ndarray
        Nested list or list of arrays containing numeric values of varying lengths.

    Returns
    -------
    np.ndarray
        Rectangular 2D NumPy array padded with NaN where needed.
    """

    if isinstance(array[0], np.ndarray):
        array = [row.tolist() for row in array]

    max_len = max(len(row) for row in array)
    return np.array([row + [np.nan] * (max_len - len(row)) for row in array])

def get_file_path():
    """
    Opens a directory selection dialog and returns the chosen directory path.

    Returns
    -------
    str
        Absolute path of the selected directory.
    """
    root = tk.Tk()
    root.withdraw()
    return filedialog.askdirectory()

def calculate_optimal_exponent(array):
    """
    Converts an Astropy Quantity array to an SI unit with an appropriate metric prefix
    based on the order of magnitude of its mean value.

    Parameters
    ----------
    array : astropy.units.Quantity
        Quantity array with a physical unit (e.g. meters).

    Returns
    -------
    astropy.units.Quantity
        The same array converted to a matching metric subunit
        (mm, µm, nm, pm, or fm).

    Raises
    ------
    Exception
        If the computed optimal exponent is outside the supported range (-15 to -3).
    """
    mean_exponent = int(math.floor(math.log10(array.value.mean())))
    optimal_unit_exponent = 3 * round(mean_exponent / 3)
    if optimal_unit_exponent == -3:
        return array.to(u.mm)
    elif optimal_unit_exponent == -6:
        return array.to(u.um)
    elif optimal_unit_exponent == -9:
        return array.to(u.nm)
    elif optimal_unit_exponent == -12:
        return array.to(u.pm)
    elif optimal_unit_exponent == -15:
        return array.to(u.fm)
    else:
        raise Exception('Optimal exponent out of range (exp < -15 or exp > -3)')


class PygwyTxt:
    """
    Handles reading, analyzing, and visualizing surface profiles
    from Gwyddion-exported `.txt` files.
    """

    def __init__(self, file_path: str, scan_size_x: float, scan_size_y: float, name: str = None, peak_finder_settings = None):
        """
        Initializes a PygwyTxt instance.

        Parameters
        ----------
        file_path : str
            Path to the input `.txt` file containing surface data.
        scan_size_x : float
            Horizontal scan size in µm.
        scan_size_y : float
            Vertical scan size in µm.
        name : str, optional
            Custom name for the dataset. Defaults to the filename.
        peak_finder_settings : PeakFinderSettings, optional
            Settings controlling the peak/valley detection algorithm.
        """
        self.__file_path = file_path
        if name is None:
            self.__name = os.path.basename(os.path.splitext(self.__file_path)[0])
        else:
            self.__name = name
        if peak_finder_settings is None:
            self.__peak_finder_settings = PeakFinderSettings()
        else:
            self.__peak_finder_settings = peak_finder_settings
        self.__scan_size_x = scan_size_x
        self.__scan_size_y = scan_size_y
        self.__scan = np.genfromtxt(file_path, delimiter='\t') * u.m
        self.__profile_line = int(self.__scan.value.shape[0]/2)
        self.__distance_per_index_x = self.__scan_size_x / self.__scan.value.shape[1]
        self.__distance_per_index_y = self.__scan_size_y / self.__scan.value.shape[0]
        self.__peak_array = None
        self.__valley_array = None
        self.__height_map, self.__period_map = self.__generate_height_and_period_map()

        self.__export_path = os.path.join(os.path.dirname(self.__file_path), 'export')
        if not os.path.exists(self.__export_path):
            os.mkdir(self.__export_path)

        self.__scan = calculate_optimal_exponent(self.__scan)
        self.__stats = self.__calculate_stats()

    def plot_scan(self, show_plot_line=True, cmap='viridis', show_title=True):
        """
        Creates and saves a heatmap of the full scan.

        Parameters
        ----------
        show_plot_line : bool, optional
            If True, marks the central profile line on the heatmap.
        cmap : str, optional
            Matplotlib colormap for visualization.
        show_title : bool, optional
            Whether to display the plot title.
        """
        fig, ax = plt.subplots()
        if show_plot_line:
            plt.axhline(y=self.__profile_line * self.__distance_per_index_y, color='red', linewidth=0.5)
        if show_title:
            ax.set_title(f'{self.__name}')
        plt.imshow(self.__scan.value, cmap=cmap, extent=(0, self.__scan_size_x, self.__scan_size_y, 0), interpolation='nearest')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(label=f'height [{str(self.__scan.unit).replace("u", "µ")}]', cax=cax)
        ax.xaxis.set_label_position('top')
        ax.xaxis.set_ticks_position('top')
        ax.text(-1.7, -0.53, "µm", ha='left', va='bottom')
        plt.tight_layout()
        plt.show()
        fig.savefig(os.path.join(self.__export_path, f'{self.__name}_heatmap.png'), bbox_inches='tight', pad_inches=0.05, dpi=300)

    def plot_profile(self, show_title=True):
        """
        Plots and saves the height profile along the central horizontal line.

        Parameters
        ----------
        show_title : bool, optional
            Whether to display the plot title.

        """
        ls = np.linspace(0, self.__scan_size_x, self.__scan.value.shape[1])
        fig, ax = plt.subplots()
        if show_title:
            ax.set_title(f'{self.__name}')
        plt.plot(ls, self.__scan[self.__profile_line])
        plt.xlabel("width [µm]")
        plt.ylabel(f"height [{str(self.__scan.unit).replace('u', 'µ')}]")
        plt.show()
        fig.savefig(os.path.join(self.__export_path, f'{self.__name}_profile.png'), bbox_inches='tight', pad_inches=0.05, dpi=300)

    def plot_profile_section(self, start: int, stop: int, line: int, show_title=True):
        """
        Plots and saves a specific section of a selected scan line.

        Parameters
        ----------
        start : int
            Starting index of the profile section.
        stop : int
            Ending index of the profile section.
        line : int
            Line index in the scan array.
        show_title : bool, optional
            Whether to display the plot title.
        """
        plot_line = self.__scan[line][start:stop+1]
        plot_line_length = len(plot_line) * (self.__scan_size_x / self.__scan.shape[1])
        ls = np.linspace(0, plot_line_length, len(plot_line))
        fig, ax = plt.subplots()
        if show_title:
            ax.set_title(f'{self.__name}')
        plt.plot(ls, plot_line)
        plt.xlabel("width [µm]")
        plt.ylabel(f"height [{str(self.__scan.unit).replace('u', 'µ')}]")
        plt.show()
        fig.savefig(os.path.join(self.__export_path, f'{self.__name}_profile_line_{line}_from_{start}_to_{stop}.png'), bbox_inches='tight', pad_inches=0.05, dpi=300)

    def plot_debug(self, line: int=None):
        """
        Visualizes peak and valley detection results for a given line.

        Parameters
        ----------
        line : int, optional
            Line index to visualize. Defaults to the central profile line.
        """
        if line is None:
            line = self.__profile_line

        plot_line = self.__scan.to(u.m).value[line]

        peaks = self.__peak_array[line][~np.isnan(self.__peak_array[line])].astype(int)
        valleys = self.__valley_array[line][~np.isnan(self.__valley_array[line])].astype(int)

        fig, ax = plt.subplots(2)
        fig.suptitle(f'{self.__name}')
        ax[0].plot(plot_line)
        ax[0].plot(peaks, plot_line[peaks], 'x')
        ax[0].plot(valleys, plot_line[valleys], 'x')
        for i, height in enumerate(self.__height_map[line]):
            if i > len(peaks) -1 or i > len(valleys)-1:
                break
            ax[0].vlines((peaks[i] + valleys[i]) / 2, plot_line[valleys[i]], height + plot_line[valleys[i]], color='red')

        ax[1].plot(plot_line)
        ax[1].plot(peaks, plot_line[peaks], 'x')
        for i, period in enumerate(self.__period_map[line]):
            if i > len(peaks) -1:
                break
            plt.hlines(plot_line[peaks[i]], peaks[i], peaks[i] + (period / self.__distance_per_index_x), color='green')

        fig.tight_layout()
        plt.show()
        fig.savefig(os.path.join(self.__export_path, f'{self.__name}_debug.png'),
                bbox_inches='tight', pad_inches=0.05, dpi=300)


    def __generate_height_and_period_map(self):
        """
        Computes height and period maps from detected peaks and valleys.

        Returns
        -------
        tuple of np.ndarray
            (height_map, period_map) — both padded with NaN values.
        """
        height_list = []
        period_list = []
        peak_list = []
        valley_list = []

        for line in self.__scan.value:
            peaks, peak_metadata = find_peaks(line,
                                              self.__peak_finder_settings.height,
                                              self.__peak_finder_settings.threshold,
                                              self.__peak_finder_settings.distance,
                                              self.__peak_finder_settings.prominence,
                                              self.__peak_finder_settings.width,
                                              self.__peak_finder_settings.wlen,
                                              self.__peak_finder_settings.rel_height,
                                              self.__peak_finder_settings.plateau_size)
            valleys, valley_metadata = find_peaks(line * -1,
                                                  self.__peak_finder_settings.height,
                                                  self.__peak_finder_settings.threshold,
                                                  self.__peak_finder_settings.distance,
                                                  self.__peak_finder_settings.prominence,
                                                  self.__peak_finder_settings.width,
                                                  self.__peak_finder_settings.wlen,
                                                  self.__peak_finder_settings.rel_height,
                                                  self.__peak_finder_settings.plateau_size)

            peak_list.append(peaks)
            valley_list.append(valleys)

            heights = []
            for i, peak in enumerate(peaks):
                if i < len(valleys):
                    height = line[peak] - line[valleys[i]]
                    heights.append(height)

            periods = []
            for i in range(len(peaks)):
                if i + 1 < len(peaks):
                    period = (peaks[i + 1] - peaks[i]) * self.__distance_per_index_x
                    periods.append(period)

            height_list.append(heights)
            period_list.append(periods)

        self.__peak_array = homogenize_array(peak_list)
        self.__valley_array = homogenize_array(valley_list)

        height_map = homogenize_array(height_list)
        period_map = homogenize_array(period_list)
        return height_map, period_map

    def __calculate_stats(self):
        """
        Computes descriptive statistics for the height and period data.

        Returns
        -------
        dict
            Statistical metrics (mean, std, min, max) in meters.
        """
        mean_height = (np.nanmean(self.__height_map) * u.m).to(self.__scan.unit)
        std_height = (np.nanstd(self.__height_map) * u.m).to(self.__scan.unit)
        mean_period = np.nanmean(self.__period_map) * u.um
        std_period = np.nanstd(self.__period_map) * u.um
        min_height = (float(np.nanmin(self.__height_map)) * u.m).to(self.__scan.unit)
        max_height = (float(np.nanmax(self.__height_map)) * u.m).to(self.__scan.unit)
        min_period = float(np.nanmin(self.__period_map)) * u.um
        max_period = float(np.nanmax(self.__period_map)) * u.um

        header = f"========== {self.__name} =========="
        footer = "=" * len(header)
        body = (f"{header}\n"
              f"height: {mean_height:.2f} +/- {std_height:.2f} \n"
              f"period: {mean_period:.2f} +/- {std_period:.2f} \n\n"
              f"min height: {min_height:.2f} \n"
              f"max height: {max_height:.2f} \n\n"
              f"min period: {min_period:.2f} \n"
              f"max period: {max_period:.2f} \n"
              f"{footer}\n")
        print(body.replace('u', 'µ'))

        stats = {
            "name": self.__name,
            "mean_height": float(mean_height.to(u.m).value),
            "std_height": float(std_height.to(u.m).value),
            "mean_period": float(mean_period.to(u.m).value),
            "std_period": float(std_period.to(u.m).value),
            "min_height": float(min_height.to(u.m).value),
            "max_height": float(max_height.to(u.m).value),
            "min_period": float(min_period.to(u.m).value),
            "max_period": float(max_period.to(u.m).value)
        }

        return stats

    def export_stats(self):
        """
        Saves the computed statistics to a JSON file inside the export directory.
        """
        export_path = os.path.join(self.__export_path, f'{self.__name}_stats.json')
        with open(export_path, 'w') as f:
            json.dump(self.__stats, f)

    def plot_heatmap(self, plot_type:int, cmap='hot'):
        """
        Creates and saves a heatmap of the height or period map.

        Parameters
        ----------
        plot_type : int
            0 for height map, 1 for period map.
        cmap : str, optional
            Matplotlib colormap for visualization (default is 'hot').
        """
        fig, ax = plt.subplots()
        if plot_type == 0:
            ax.set_title(f"{self.__name} Height Heatmap")
            plt.imshow(self.__height_map, cmap=cmap, interpolation='none', aspect='auto')
            name_appendix = 'height_heat_map'
        elif plot_type == 1:
            ax.set_title(f"{self.__name} Period Heatmap")
            plt.imshow(self.__period_map, cmap=cmap, interpolation='none', aspect='auto')
            name_appendix = 'period_heat_map'
        plt.colorbar()
        plt.show()
        fig.tight_layout()
        fig.savefig(os.path.join(self.__export_path, f'{self.__name}_{name_appendix}.png'),
                    bbox_inches='tight', pad_inches=0.05, dpi=300)


class StatJson:
    """
    Loads, aggregates, and plots statistical results from multiple JSON files
    generated by `PygwyTxt`.
    """


    def __init__(self, base_path: str):
        """
        Loads all JSON statistic files from the given directory.

        Parameters
        ----------
        base_path : str
            Path containing the JSON files.
        """
        self.__base_path = base_path
        self.__export_path = os.path.join(base_path, 'stat_plots')
        if not os.path.exists(self.__export_path):
            os.makedirs(self.__export_path)

        self.__lamda = lambda x: int(re.findall(r'\d+', os.path.split(x)[-1])[0])
        file_list = sorted(glob2.glob(os.path.join(base_path, '*[!exclude]*.json')), key=self.__lamda)

        self.__stat_list = []
        for file_path in file_list:
            with open(file_path, 'r') as file:
                data = json.load(file)
                self.__stat_list.append(data)

        self.__plot_data_height = None
        self.__plot_data_period = None

    def plot(self, plot_type: int, x_lable: str, x_unit: str, plot_name_appendix='', model=None, params=None, show_title=True):
        """
        Creates a plot showing mean and standard deviation of height or period data.
        Optionally applies a fitted model to the data.

        Parameters
        ----------
        plot_type : int
            0 for height, 1 for period.
        x_lable : str
            Label for the x-axis.
        x_unit : str
            Unit for the x-axis.
        plot_name_appendix : str, optional
            String appended to the filename and plot title.
        model : lmfit.Model, optional
            Fit model to apply to the mean values.
        params : lmfit.Parameters, optional
            Parameters for the fit model.
        show_title : bool, optional
            Whether to display the plot title.
        """
        if plot_type == 0:
            plot_name = "mean height"
            mean_key = 'mean_height'
            std_key = 'std_height'
            y_label = 'height'
        elif plot_type == 1:
            plot_name = "mean period"
            mean_key = 'mean_period'
            std_key = 'std_period'
            y_label = 'period'
        else:
            raise ValueError('plot_type must be 0 or 1')

        mean = []
        std = []
        x_values = []
        for stat in self.__stat_list:
            mean.append(stat[mean_key])
            std.append(stat[std_key])
            x_values.append(self.__lamda(stat['name']))

        if plot_type == 0:
            self.__plot_data_height = [x_values, mean, std]
        elif plot_type == 1:
            self.__plot_data_period = [x_values, mean, std]

        mean = calculate_optimal_exponent(mean * u.m)
        std = (std * u.m).to(mean.unit)

        fig, ax = plt.subplots()

        if model is not None and params is not None:
            assert isinstance(model, Model) == True
            result = model.fit(mean.value, params, x=x_values, weights=1/std.value)
            print(result.fit_report())
            plt.plot(x_values, result.best_fit, color='red', label='best fit')

            if plot_type == 0:
                self.__plot_data_height.append((result.best_fit * mean.unit).to(u.m).value)
            elif plot_type == 1:
                self.__plot_data_period.append((result.best_fit * mean.unit).to(u.m).value)

            with open(os.path.join(self.__export_path, f'fit_report_{plot_name}_{plot_name_appendix}.txt'), 'w') as file:
                file.write(result.fit_report())

        if show_title:
            ax.set_title(f'{plot_name} {plot_name_appendix}')
        ax.plot(x_values, mean, 'o', zorder=2, label='mean')
        ax.errorbar(x_values, mean, yerr=std, fmt='none', capsize=5, ecolor='black', elinewidth=1, zorder=1, label='std')
        plt.xlabel(f"{x_lable} [{x_unit}]")
        plt.ylabel(f"{y_label} [{str(mean.unit).replace('u', 'µ')}]")
        plt.legend()
        plt.show()
        if plot_name_appendix == '':
            fig.savefig(os.path.join(self.__export_path, f'{plot_name}.png'), bbox_inches='tight',
                        pad_inches=0.05, dpi=300)
        else:
            fig.savefig(os.path.join(self.__export_path, f'{plot_name}_{plot_name_appendix}.png'), bbox_inches='tight',
                    pad_inches=0.05, dpi=300)

    def export_plot_data(self, plot_type: int):
        """
        Exports the x-values, mean, standard deviation, and fit (if available)
        used in the plots to a CSV file.

        Parameters
        ----------
        plot_type : int
            0 for height data, 1 for period data.
        """
        if plot_type == 0:
            plot_data = self.__plot_data_height
            name = 'data_height_plot.csv'
        elif plot_type == 1:
            plot_data = self.__plot_data_period
            name = 'data_period_plot.csv'
        else:
            raise ValueError('plot_type must be 0 or 1')

        if plot_data is not None:
            with open(os.path.join(self.__export_path, name), 'w', newline='') as file:
                writer = csv.writer(file)
                if len(plot_data) == 3:
                    writer.writerow(['x', 'y', 'std'])
                    for i in range(len(plot_data[0])):
                        writer.writerow([plot_data[0][i], plot_data[1][i], plot_data[2][i]])
                elif len(plot_data) == 4:
                    writer.writerow(['x', 'y', 'std', 'best fit'])
                    for i in range(len(plot_data[0])):
                        writer.writerow([plot_data[0][i], plot_data[1][i], plot_data[2][i], plot_data[3][i]])

class PeakFinderSettings:
    """
    Defines configurable parameters for `scipy.signal.find_peaks`.
    """
    def __init__(self,
                 height=None,
                 threshold=None,
                 distance=None,
                 prominence=None,
                 width=None,
                 wlen=None,
                 rel_height=None,
                 plateau_size=None):
        """
        Initializes the peak-finding parameter object.

        Parameters
        ----------
        height : float or sequence, optional
            Required height of peaks.
        threshold : float or sequence, optional
            Required vertical difference between peaks and neighbors.
        distance : int, optional
            Minimum horizontal distance between peaks.
        prominence : float or sequence, optional
            Required prominence of peaks.
        width : float or sequence, optional
            Required width of peaks.
        wlen : int, optional
            Window length for peak prominence evaluation.
        rel_height : float, optional
            Relative height at which the peak width is measured.
        plateau_size : float or sequence, optional
            Range of flat peak plateaus.
        """
        assert isinstance(height, (numbers.Number, np.ndarray, Sequence)) or height is None
        assert isinstance(threshold, (numbers.Number, np.ndarray, Sequence)) or threshold is None
        assert isinstance(distance, numbers.Number) or distance is None
        assert isinstance(prominence, (numbers.Number, np.ndarray, Sequence)) or prominence is None
        assert isinstance(width, (numbers.Number, np.ndarray, Sequence)) or width is None
        assert wlen is int or wlen is None
        assert rel_height is float or rel_height is None
        assert isinstance(plateau_size, (numbers.Number, np.ndarray, Sequence)) or plateau_size is None

        self.height = height
        self.threshold = threshold
        self.distance = distance
        self.prominence = prominence
        self.width = width
        self.wlen = wlen
        self.rel_height = rel_height
        self.plateau_size = plateau_size
