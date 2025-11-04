import pygwytxt
import glob2

basepath = r''

file_list = glob2.glob('[!stats]*.txt')

for file_path in file_list:
    scan = pygwytxt.PygwyTxt(file_path, 20, 5)
    scan.plot_scan()
    scan.plot_profile()
    scan.export_stats()