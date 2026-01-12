import pygwy_txt_analysis


scan_size_x = 25
scan_size_y = 25
peak_finder_settings = pygwy_txt_analysis.PeakFinderSettings(prominence=0.3e-7)

basepath = pygwy_txt_analysis.get_file_path()
scan = pygwy_txt_analysis.PygwyTxt(basepath, scan_size_x, scan_size_y, peak_finder_settings=peak_finder_settings)
scan.plot_scan()
scan.plot_profile()
scan.export_stats()
scan.plot_debug()
scan.plot_heatmap(0)
scan.plot_heatmap(1)

