import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from astropy import units as u
import math
import os
from scipy.signal import find_peaks
import json


class PygwyTxt:
    def __init__(self, file_path: str, scan_size_x: float, scan_size_y: float, name:str=None):
        self.__file_path = file_path
        if name is None:
            self.__name = os.path.basename(os.path.splitext(self.__file_path)[0])
        else:
            self.__name = name
        self.__scan_size_x = scan_size_x
        self.__scan_size_y = scan_size_y
        self.__scan = np.genfromtxt(file_path, delimiter='\t') * u.m
        self.__profile_line = int(self.__scan.value.shape[0]/2)
        self.__distance_per_index_x = self.__scan_size_x / self.__scan.value.shape[1]
        self.__distance_per_index_y = self.__scan_size_y / self.__scan.value.shape[0]
        self.__height_map, self.__period_map = self.__generate_height_and_period_map()

        mean_exponent = int(math.floor(math.log10(self.__scan.value.mean())))
        optimal_unit_exponent = 3 * round(mean_exponent / 3)
        if optimal_unit_exponent == -3:
            self.__scan = self.__scan.to(u.mm)
        if optimal_unit_exponent == -6:
            self.__scan = self.__scan.to(u.um)
        if optimal_unit_exponent == -9:
            self.__scan = self.__scan.to(u.nm)
        if optimal_unit_exponent == -12:
            self.__scan = self.__scan.to(u.pm)
        if optimal_unit_exponent == -15:
            self.__scan = self.__scan.to(u.fm)

        self.__stats = self.__calculate_stats()

    def plot_scan(self, show_plot_line=True, cmap='viridis'):
        fig, ax = plt.subplots()
        if show_plot_line:
            plt.axhline(y=self.__profile_line * self.__distance_per_index_y, color='red', linewidth=0.5)
        ax.set_title(f'{self.__name}')
        plt.imshow(self.__scan.value, cmap=cmap, extent=(0, self.__scan_size_x, self.__scan_size_y, 0), interpolation='nearest')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(label=f'height [{str(self.__scan.unit)}]', cax=cax)
        ax.set_title(f'{self.__name}')
        ax.xaxis.set_label_position('top')
        ax.xaxis.set_ticks_position('top')
        ax.text(-1.7, -0.53, "µm", ha='left', va='bottom')
        plt.tight_layout()
        plt.show()
        fig.savefig(os.path.join(os.path.dirname(self.__file_path), f'{self.__name}_heatmap.png'), bbox_inches='tight', pad_inches=0.05, dpi=300)

    def plot_profile(self):
        ls = np.linspace(0, self.__scan_size_x, self.__scan.value.shape[1])
        fig, ax = plt.subplots()
        ax.set_title(f'{self.__name}')
        plt.plot(ls, self.__scan[self.__profile_line])
        plt.xlabel("width [µm]")
        plt.ylabel(f"height [{str(self.__scan.unit)}]")
        plt.show()
        fig.savefig(os.path.join(os.path.dirname(self.__file_path), f'{self.__name}_profile.png'), bbox_inches='tight', pad_inches=0.05, dpi=300)

    def plot_profile_section(self, start:int, stop:int, line:int):
        plot_line = self.__scan[line][start:stop+1]
        plot_line_length = len(plot_line) * (self.__scan_size_x / self.__scan.shape[1])

        ls = np.linspace(0, plot_line_length, len(plot_line))
        fig, ax = plt.subplots()
        ax.set_title(f'{self.__name}')
        plt.plot(ls, plot_line)
        plt.xlabel("width [µm]")
        plt.ylabel(f"height [{str(self.__scan.unit)}]")
        plt.show()
        fig.savefig(os.path.join(os.path.dirname(self.__file_path), f'{self.__name}_profile_line_{line}_from_{start}_to_{stop}.png'), bbox_inches='tight', pad_inches=0.05, dpi=300)

    def __generate_height_and_period_map(self):
        height_list = []
        period_list = []
        for line in self.__scan.value:
            peaks, peak_metadata = find_peaks(line)
            valleys, valley_metadata = find_peaks(line * -1)
            heights = line[peaks] - line[valleys]
            periods = []
            for i in range(len(peaks)):
                if i + 1 < len(peaks):
                    period = (peaks[i + 1] - peaks[i]) * self.__distance_per_index_x
                    periods.append(period)
            height_list.append(heights)
            period_list.append(periods)

        height_map = np.array(height_list)
        period_map = np.array(period_list)

        return height_map, period_map

    def __calculate_stats(self):
        mean_height = (np.mean(self.__height_map) * u.m).to(self.__scan.unit)
        std_height = (np.std(self.__height_map) * u.m).to(self.__scan.unit)

        mean_period = np.mean(self.__period_map) * u.um
        std_period = np.std(self.__period_map) * u.um

        min_height = (float(np.min(self.__height_map)) * u.m).to(self.__scan.unit)
        max_height = (float(np.max(self.__height_map)) * u.m).to(self.__scan.unit)

        min_period = float(np.min(self.__period_map)) * u.um
        max_period = float(np.max(self.__period_map)) * u.um

        header = f"========== {self.__name} =========="
        footer = "=" * len(header)

        print(f"{header}\n"
              f"height: {mean_height:.2f} +/- {std_height:.2f} \n"
              f"period: {mean_period:.2f} +/- {std_period:.2f} \n\n"
              f"min height: {min_height:.2f} \n"
              f"max height: {max_height:.2f} \n\n"
              f"min period: {min_period:.2f} \n"
              f"max period: {max_period:.2f} \n"
              f"{footer}\n\n\n")

        stats = {
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
        export_path = os.path.join(os.path.dirname(self.__file_path), f'stats_{self.__name}.json')
        with open(export_path, 'w') as f:
            json.dump(self.__stats, f)


scan = PygwyTxt(r"C:\Users\Mika Music\Data\251029_WNE_pygwy\gwy\ref.txt", 20, 5)

scan.export_stats()

