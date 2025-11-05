import pygwytxt
import glob2
import tkinter as tk
from tkinter import filedialog
import os


root = tk.Tk()
root.withdraw()

basepath = filedialog.askdirectory()

file_list = glob2.glob(os.path.join(basepath, '*[!exclude]*.txt'))

for file_path in file_list:
    scan = pygwytxt.PygwyTxt(file_path, 20, 5)
    scan.plot_scan()
    scan.plot_profile()
    scan.export_stats()