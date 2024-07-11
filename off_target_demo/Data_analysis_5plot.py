"""
Analyze the generated docking scores for visualization
"""

import csv
import numpy as np
import pandas as pd

data_list = []
with open(r'./vina_scores.csv') as f:
    reader = csv.reader(f)
    headers = next(reader)  # this will jump the 1st row
    for row in reader:
        # source pocket id, generated molecular id, target pocket id, if the pair, score
        data_list.append((row[0], row[1], row[2], row[6], row[7]))
    # print(data_list)
    # ('26', '10', '25', '0', '-9.2') is ele in data_list

# we take 20 source pocket, to generated 100 +- mols for each of them
exp_data_list = []
for data_item in data_list:
    if int(data_item[0]) < 20:
        exp_data_list.append((int(data_item[0]), int(data_item[1]), int(data_item[2]),
                              int(data_item[3]), float(data_item[4])))
# exp_data_list, simulate the animal env.

"""scores for heatmap, rows are all mols, cols are all pockets, 100 pockets"""
heatmap = []
for i in range(100):
    heatmap.append([])

for data_item in exp_data_list:
    pock_index = int(data_item[2])
    score = data_item[4]
    heatmap[pock_index].append(score)

print(heatmap)
# print(111)


csv_path = r'./off_target_5plot.csv'

first_row = ['x0_pocket', 'x1_pocket', 'x2_pocket', 'x3_pocket', 'x4_pocket', 'x5_pocket', 'x6_pocket', 'x7_pocket', 'x8_pocket', 'x9_pocket',
             'x10_pocket', 'x11_pocket', 'x12_pocket', 'x13_pocket', 'x14_pocket', 'x15_pocket', 'x16_pocket', 'x17_pocket', 'x18_pocket', 'x19_pocket',
             'x20_pocket', 'x21_pocket', 'x22_pocket', 'x23_pocket', 'x24_pocket', 'x25_pocket', 'x26_pocket', 'x27_pocket', 'x28_pocket', 'x29_pocket',
             'x30_pocket', 'x31_pocket', 'x32_pocket', 'x33_pocket', 'x34_pocket', 'x35_pocket', 'x36_pocket', 'x37_pocket', 'x38_pocket', 'x39_pocket',
             'x40_pocket', 'x41_pocket', 'x42_pocket', 'x43_pocket', 'x44_pocket', 'x45_pocket', 'x46_pocket', 'x47_pocket', 'x48_pocket', 'x49_pocket',
             'x50_pocket', 'x51_pocket', 'x52_pocket', 'x53_pocket', 'x54_pocket', 'x55_pocket', 'x56_pocket', 'x57_pocket', 'x58_pocket', 'x59_pocket',
             'x60_pocket', 'x61_pocket', 'x62_pocket', 'x63_pocket', 'x64_pocket', 'x65_pocket', 'x66_pocket', 'x67_pocket', 'x68_pocket', 'x69_pocket',
             'x70_pocket', 'x71_pocket', 'x72_pocket', 'x73_pocket', 'x74_pocket', 'x75_pocket', 'x76_pocket', 'x77_pocket', 'x78_pocket', 'x79_pocket',
             'x80_pocket', 'x81_pocket', 'x82_pocket', 'x83_pocket', 'x84_pocket', 'x85_pocket', 'x86_pocket', 'x87_pocket', 'x88_pocket', 'x89_pocket',
             'x90_pocket', 'x91_pocket', 'x92_pocket', 'x93_pocket', 'x94_pocket', 'x95_pocket', 'x96_pocket', 'x97_pocket', 'x98_pocket', 'x99_pocket']
with open(csv_path, "a", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(first_row)
    # for i in range(max_num):
    for i in range(len(heatmap[0])):
        row = []
        for j in range(0, 100):
            row.append(heatmap[j][i])
        writer.writerow(row)

