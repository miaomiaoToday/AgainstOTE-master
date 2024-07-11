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

"""1. & 2. target scores"""
x0_target_list = []
x1_other_list = []
for data_item in exp_data_list:
    if data_item[3] == 1:
        # print(data_item[4])
        x0_target_list.append(data_item[4])
    else:
        if data_item[4] <= 0:
            x1_other_list.append(data_item[4])
        else:
            print(data_item[4])
# print(x0_target_list)

for i in range(len(x1_other_list)):
    try:
        tt = x0_target_list[i]
    except:
        x0_target_list.append(' ')

print(len(x0_target_list), len(x1_other_list))

csv_path = r'./off_target_4plot.csv'

first_row = ['x0_target_list', 'x1_other_list']
with open(csv_path, "a", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(first_row)
    # for i in range(max_num):
    for i in range(len(x1_other_list)):
        row = [x0_target_list[i], x1_other_list[i]]
        writer.writerow(row)

