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
    Makes a 2D array rectangular by padding shorter rows with NaN values.

    Parameters
    ----------
    array : list of lists
        Nested list containing numeric values of varying lengths.

    Returns
    -------
    np.ndarray
        2D NumPy array with equal row lengths, padded with NaN values.
    """
    max_len = max(len(row) for row in array)
    return np.array([row + [np.nan] * (max_len - len(row)) for row in array])

def get_file_path():
    """
    Opens a directory selection dialog and returns the chosen path.

    Returns
    -------
    str
        Path to the selected directory.
    """
    root = tk.Tk()
    root.withdraw()
    return filedialog.askdirectory()

def calculate_optimal_exponent(array):
    """
    Converts an Astropy Quantity array to an appropriate SI unit prefix
    based on the order of magnitude of its mean value.

    Parameters
    ----------
    array : astropy.units.Quantity
        Input array with a physical unit.

    Returns
    -------
    astropy.units.Quantity
        The same array converted to a suitable metric unit
        (mm, µm, nm, pm, or fm).

    Raises
    ------
    Exception
        If the optimal exponent is outside the supported range (-15 to -3).
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
    Class for analyzing and visualizing surface profiles from .txt files
    (e.g., Gwyddion exports).
    """

    def __init__(self, file_path: str, scan_size_x: float, scan_size_y: float, name: str = None, peak_finder_settings = None):
        """
        Initializes a PygwyTxt instance.

        Parameters
        ----------
        file_path : str
            Path to the input .txt file.
        scan_size_x : float
            Scan width in µm.
        scan_size_y : float
            Scan height in µm.
        name : str, optional
            Custom name for the object.
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
        self.__height_map, self.__period_map = self.__generate_height_and_period_map()

        self.__export_path = os.path.join(os.path.dirname(self.__file_path), 'export')
        if not os.path.exists(self.__export_path):
            os.mkdir(self.__export_path)

        self.__scan = calculate_optimal_exponent(self.__scan)
        self.__stats = self.__calculate_stats()

    def plot_scan(self, show_plot_line=True, cmap='viridis'):
        """
        Creates and saves a heatmap of the scan.

        Parameters
        ----------
        show_plot_line : bool, optional
            Whether to display the profile line in the center of the scan (default: True).
        cmap : str, optional
            Colormap used for visualization.
        """
        fig, ax = plt.subplots()
        if show_plot_line:
            plt.axhline(y=self.__profile_line * self.__distance_per_index_y, color='red', linewidth=0.5)
        ax.set_title(f'{self.__name}')
        plt.imshow(self.__scan.value, cmap=cmap, extent=(0, self.__scan_size_x, self.__scan_size_y, 0), interpolation='nearest')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(label=f'height [{str(self.__scan.unit).replace("u", "µ")}]', cax=cax)
        ax.set_title(f'{self.__name}')
        ax.xaxis.set_label_position('top')
        ax.xaxis.set_ticks_position('top')
        ax.text(-1.7, -0.53, "µm", ha='left', va='bottom')
        plt.tight_layout()
        plt.show()
        fig.savefig(os.path.join(self.__export_path, f'{self.__name}_heatmap.png'), bbox_inches='tight', pad_inches=0.05, dpi=300)

    def plot_profile(self):
        """
        Plots and saves the height profile along the central scan line.
        """
        ls = np.linspace(0, self.__scan_size_x, self.__scan.value.shape[1])
        fig, ax = plt.subplots()
        ax.set_title(f'{self.__name}')
        plt.plot(ls, self.__scan[self.__profile_line])
        plt.xlabel("width [µm]")
        plt.ylabel(f"height [{str(self.__scan.unit).replace('u', 'µ')}]")
        plt.show()
        fig.savefig(os.path.join(self.__export_path, f'{self.__name}_profile.png'), bbox_inches='tight', pad_inches=0.05, dpi=300)

    def plot_profile_section(self, start: int, stop: int, line: int):
        """
        Plots and saves a defined section of a specific profile line.

        Parameters
        ----------
        start : int
            Start index of the section.
        stop : int
            End index of the section.
        line : int
            Line index in the scan.
        """
        plot_line = self.__scan[line][start:stop+1]
        plot_line_length = len(plot_line) * (self.__scan_size_x / self.__scan.shape[1])
        ls = np.linspace(0, plot_line_length, len(plot_line))
        fig, ax = plt.subplots()
        ax.set_title(f'{self.__name}')
        plt.plot(ls, plot_line)
        plt.xlabel("width [µm]")
        plt.ylabel(f"height [{str(self.__scan.unit).replace('u', 'µ')}]")
        plt.show()
        fig.savefig(os.path.join(self.__export_path, f'{self.__name}_profile_line_{line}_from_{start}_to_{stop}.png'), bbox_inches='tight', pad_inches=0.05, dpi=300)

    def __generate_height_and_period_map(self):
        """
        Generates maps of peak heights and periods for each scan line.

        Returns
        -------
        tuple of np.ndarray
            (height_map, period_map), both padded with NaN values.
        """
        height_list = []
        period_list = []

        for line in self.__scan.value:
            peaks, peak_metadata = find_peaks(line, self.__peak_finder_settings.height, )
            valleys, valley_metadata = find_peaks(line * -1)

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

        height_map = np.array(homogenize_array(height_list))
        period_map = np.array(homogenize_array(period_list))
        return height_map, period_map

    def __calculate_stats(self):
        """
        Calculates statistical metrics for height and period data
        (mean, standard deviation, min, max).

        Returns
        -------
        dict
            Dictionary containing the computed values in meters.
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
        Exports the computed statistics as a JSON file in the export directory.
        """
        export_path = os.path.join(self.__export_path, f'{self.__name}_stats.json')
        with open(export_path, 'w') as f:
            json.dump(self.__stats, f)


class StatJson:
    """
    Class for aggregating and analyzing statistics from multiple JSON files
    produced by PygwyTxt.
    """

    def __init__(self, base_path: str):
        """
        Initializes a StatJson instance and loads all JSON statistic files.

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

    def plot(self, plot_type: int, x_lable: str, x_unit: str, plot_name_appendix='', model=None, params=None):
        """
        Creates a plot for mean and standard deviation of either height or period.
        Optionally applies a fit model to the data.

        Parameters
        ----------
        plot_type : int
            0 = height, 1 = period.
        x_lable : str
            X-axis label.
        x_unit : str
            X-axis unit.
        plot_name_appendix : str, optional
            Optional text to append to the plot name and title.
        model : lmfit.Model, optional
            Model used for curve fitting.
        params : lmfit.Parameters, optional
            Parameters of the fit model.
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
        Exports the plot data (x, mean, std, and optional fit) to a CSV file.

        Parameters
        ----------
        plot_type : int
            0 = height, 1 = period.
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
    def __init__(self,
                 height=None,
                 threshold=None,
                 distance=None,
                 prominence=None,
                 width=None,
                 wlen=None,
                 rel_height=None,
                 plateau_size=None):
        assert isinstance(height, (numbers.Number, np.ndarray, Sequence)) or height is None
        assert isinstance(threshold, (numbers.Number, np.ndarray, Sequence)) or threshold is None
        assert isinstance(distance, numbers.Number) or distance is None
        assert isinstance(prominence, (numbers.Number, np.ndarray, Sequence)) or prominence is None
        assert isinstance(width, (numbers.Number, np.ndarray, Sequence)) or width is None
        assert wlen is int or wlen is None
        assert rel_height is float or rel_height is None
        assert isinstance(plateau_size, (numbers.Number, np.ndarray, Sequence)) or plateau_size is None

        self.__height = height
        self.__threshold = threshold
        self.__distance = distance
        self.__prominence = prominence
        self.__width = width
        self.__wlen = wlen
        self.__rel_height = rel_height
        self.__plateau_size = plateau_size
