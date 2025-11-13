import pygwy_analysis
import glob2
import os
from lmfit import Model

basepath = pygwy_analysis.get_file_path()
peak_finder_settings = pygwy_analysis.PeakFinderSettings(prominence=0.3e-7)

file_list_txt = glob2.glob(os.path.join(basepath, '*[!exclude]*.txt'))
for file_path in file_list_txt:
    scan = pygwy_analysis.PygwyTxt(file_path, 25, 25, peak_finder_settings=peak_finder_settings)
    scan.plot_scan()
    scan.plot_profile()
    scan.export_stats()
    scan.plot_debug()

stats = pygwy_analysis.StatJson(os.path.join(basepath, 'export'))
stats.plot(0, 'x_max', 'px')
stats.plot(1, 'x_max', 'px')
stats.export_plot_data(0)
stats.export_plot_data(1)
