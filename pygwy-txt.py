import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from astropy import units as u
import math
import os


class PygwyTxt:
    def __init__(self, file_path: str, scan_size_x: float, scan_size_y: float, name:str=None, cmap='afmhot'):
        self.__file_path = file_path
        if name is None:
            self.__name = os.path.basename(os.path.splitext(self.__file_path)[0])
        else:
            self.__name = name
        self.__scan_size_x = scan_size_x
        self.__scan_size_y = scan_size_y
        self.__cmap = cmap
        self.__scan = np.genfromtxt(file_path, delimiter='\t') * u.m

        print(self.__scan.value.shape)

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

    def plot_scan(self):
        fig, ax = plt.subplots()
        ax.set_title(f'{self.__name}')
        plt.imshow(self.__scan.value, cmap=self.__cmap, extent=(0, self.__scan_size_x, self.__scan_size_y, 0), interpolation='nearest')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(label=f'height [{str(self.__scan.unit)}]', cax=cax)
        ax.set_title(f'{self.__name}')
        ax.xaxis.set_label_position('top')
        ax.xaxis.set_ticks_position('top')
        ax.text(-1.7, -0.53, "µm", ha='left', va='bottom')
        plt.tight_layout()
        plt.show()
        fig.savefig(os.path.join(os.path.dirname(self.__file_path), f'{self.__name}.png'), bbox_inches='tight', pad_inches=0.05, dpi=300)



scan = PygwyTxt(r"C:\Users\Mika Music\Data\251029_WNE_pygwy\gwy\ref.txt", 20, 5, cmap='viridis')

scan.plot_scan()
