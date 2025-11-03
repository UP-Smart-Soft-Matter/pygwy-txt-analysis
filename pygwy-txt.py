import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from astropy import units as u
import math
import os
from scipy.signal import find_peaks


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
        self.__calculate_stats()

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
        height_mean = []
        height_std = []
        for line in self.__height_map:
            mean = np.mean(line)
            std = np.std(line, ddof=1)
            height_mean.append(mean)
            height_std.append(std)

        period_mean = []
        period_std = []
        for line in self.__period_map:
            mean = np.mean(line)
            std = np.std(line, ddof=1)
            period_mean.append(mean)
            period_std.append(std)

        global_height_mean = np.mean(height_mean)
        global_period_mean = np.mean(period_mean)

        print(f'Global mean height: {global_height_mean}, Global mean period: {global_period_mean}')


scan = PygwyTxt(r"C:\Users\Mika Music\Data\251029_WNE_pygwy\gwy\ref.txt", 20, 5)

