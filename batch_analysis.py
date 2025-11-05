import pygwy_analysis
import glob2
import os

basepath = pygwy_analysis.get_file_path()

file_list = glob2.glob(os.path.join(basepath, '*[!exclude]*.txt'))
for file_path in file_list:
    scan = pygwy_analysis.PygwyTxt(file_path, 20, 5)
    scan.plot_scan()
    scan.plot_profile()
    scan.export_stats()