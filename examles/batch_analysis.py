import pygwy_txt_analysis
import glob2
import os

scan_size_x = 20
scan_size_y = 20
x_lable = 'time'
x_unit = 's'
peak_finder_settings = pygwy_txt_analysis.PeakFinderSettings(prominence=0.3e-7)

basepath = pygwy_txt_analysis.get_folder_path()
file_list_txt = glob2.glob(os.path.join(basepath, '*[!exclude]*.txt'))
for file_path in file_list_txt:
    scan = pygwy_txt_analysis.PygwyTxt(file_path, scan_size_x, scan_size_y, peak_finder_settings=peak_finder_settings)
    scan.plot_scan()
    scan.plot_profile()
    scan.export_stats()
    scan.plot_histogram(0)
    scan.plot_histogram(1)
    scan.plot_debug()
    scan.plot_heatmap(0)
    scan.plot_heatmap(1)

stats = pygwy_txt_analysis.StatJson(os.path.join(basepath, 'export'))
stats.plot(0, x_lable, x_unit)
stats.plot(1, x_lable, x_unit)
stats.export_plot_data(0)
stats.export_plot_data(1)
