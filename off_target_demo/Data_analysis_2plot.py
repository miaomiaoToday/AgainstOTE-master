"""
Analyze the generated docking scores for visualization
"""
import csv
import numpy as np


def pocket_each_ana(x0_list):
    x0_target = []
    for data_item in x0_list:
        if data_item[3] == 1:
            x0_target.append(data_item[4])
    print(x0_target)

    """7. Top 10 mean, ..."""
    temp_list = []
    x10_max_list = []
    for data_item in x0_list:
        if data_item[0] == data_item[1] == data_item[2] == 0:
            # temp_list.append(data_item[4])
            pass
        elif 0 < data_item[2] < 100:
            if data_item[3] == 1:
                continue
            temp_list.append(data_item[4])
        elif data_item[2] == 0:
            if len(temp_list) == 0 and data_item[3] == 0:
                temp_list.append(data_item[4])
                continue
            if len(temp_list) == 0 and data_item[3] == 1:
                continue
            temp_list = sorted(temp_list)
            max_val = (temp_list[0] + temp_list[1] + temp_list[2] + temp_list[3] + temp_list[4] + temp_list[5] + temp_list[6] + temp_list[7] + temp_list[8] + temp_list[9]) / 10
            x10_max_list.append(max_val)
            temp_list = []
            # temp_list.append(data_item[4])

    temp_list = sorted(temp_list)
    max_val = (temp_list[0] + temp_list[1] + temp_list[2] + temp_list[3] + temp_list[4] + temp_list[5] + temp_list[6] + temp_list[7] + temp_list[8] + temp_list[9]) / 10
    x10_max_list.append(max_val)
    print(x10_max_list)

    """8. Top 20 mean, ..."""
    temp_list = []
    x20_max_list = []
    for data_item in x0_list:
        if data_item[0] == data_item[1] == data_item[2] == 0:
            # temp_list.append(data_item[4])
            pass
        elif 0 < data_item[2] < 100:
            if data_item[3] == 1:
                continue
            temp_list.append(data_item[4])
        elif data_item[2] == 0:
            if len(temp_list) == 0 and data_item[3] == 0:
                temp_list.append(data_item[4])
                continue
            if len(temp_list) == 0 and data_item[3] == 1:
                continue
            temp_list = sorted(temp_list)
            max_val = 0
            for i in range(20):
                max_val += temp_list[i]
            max_val = max_val / 20
            # max_val = (temp_list[0] + temp_list[1] + temp_list[2] + temp_list[3] + temp_list[4] + temp_list[5] + temp_list[6] + temp_list[7] + temp_list[8] + temp_list[9]) / 10
            x20_max_list.append(max_val)
            temp_list = []
            # temp_list.append(data_item[4])

    temp_list = sorted(temp_list)
    max_val = 0
    for i in range(20):
        max_val += temp_list[i]
    max_val = max_val / 20
    x20_max_list.append(max_val)
    print(x20_max_list)

    return x0_target, x10_max_list, x20_max_list


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

# split as mols of 20 pockets
mols_pock = []
x0_list, x1_list, x2_list, x3_list, x4_list, x5_list, x6_list, x7_list, x8_list, x9_list, x10_list, x11_list, \
    x12_list, x13_list, x14_list, x15_list, x16_list, x17_list, x18_list, x19_list = [], [], [], [], [], [], [], [], [], [], \
    [], [], [], [], [], [], [], [], [], [] # pocket0_target
for data_item in exp_data_list:
    if data_item[0] == 0:
        x0_list.append(data_item)
    if data_item[0] == 1:
        x1_list.append(data_item)
    if data_item[0] == 2:
        x2_list.append(data_item)
    if data_item[0] == 3:
        x3_list.append(data_item)
    if data_item[0] == 4:
        x4_list.append(data_item)
    if data_item[0] == 5:
        x5_list.append(data_item)
    if data_item[0] == 6:
        x6_list.append(data_item)
    if data_item[0] == 7:
        x7_list.append(data_item)
    if data_item[0] == 8:
        x8_list.append(data_item)
    if data_item[0] == 9:
        x9_list.append(data_item)
    if data_item[0] == 10:
        x10_list.append(data_item)
    if data_item[0] == 11:
        x11_list.append(data_item)
    if data_item[0] == 12:
        x12_list.append(data_item)
    if data_item[0] == 13:
        x13_list.append(data_item)
    if data_item[0] == 14:
        x14_list.append(data_item)
    if data_item[0] == 15:
        x15_list.append(data_item)
    if data_item[0] == 16:
        x16_list.append(data_item)
    if data_item[0] == 17:
        x17_list.append(data_item)
    if data_item[0] == 18:
        x18_list.append(data_item)
    if data_item[0] == 19:
        x19_list.append(data_item)

x0_target, x0_x10_max_list, x0_x20_max_list = pocket_each_ana(x0_list)
x1_target, x1_x10_max_list, x1_x20_max_list = pocket_each_ana(x1_list)
x2_target, x2_x10_max_list, x2_x20_max_list = pocket_each_ana(x2_list)
x3_target, x3_x10_max_list, x3_x20_max_list = pocket_each_ana(x3_list)
x4_target, x4_x10_max_list, x4_x20_max_list = pocket_each_ana(x4_list)
x5_target, x5_x10_max_list, x5_x20_max_list = pocket_each_ana(x5_list)
x6_target, x6_x10_max_list, x6_x20_max_list = pocket_each_ana(x6_list)
x7_target, x7_x10_max_list, x7_x20_max_list = pocket_each_ana(x7_list)
x8_target, x8_x10_max_list, x8_x20_max_list = pocket_each_ana(x8_list)
x9_target, x9_x10_max_list, x9_x20_max_list = pocket_each_ana(x9_list)
x10_target, x10_x10_max_list, x10_x20_max_list = pocket_each_ana(x10_list)
x11_target, x11_x10_max_list, x11_x20_max_list = pocket_each_ana(x11_list)
x12_target, x12_x10_max_list, x12_x20_max_list = pocket_each_ana(x12_list)
x13_target, x13_x10_max_list, x13_x20_max_list = pocket_each_ana(x13_list)
x14_target, x14_x10_max_list, x14_x20_max_list = pocket_each_ana(x14_list)
x15_target, x15_x10_max_list, x15_x20_max_list = pocket_each_ana(x15_list)
x16_target, x16_x10_max_list, x16_x20_max_list = pocket_each_ana(x16_list)
x17_target, x17_x10_max_list, x17_x20_max_list = pocket_each_ana(x17_list)
x18_target, x18_x10_max_list, x18_x20_max_list = pocket_each_ana(x18_list)
x19_target, x19_x10_max_list, x19_x20_max_list = pocket_each_ana(x19_list)


"""summary 1-8. 8 list"""
# x0_list, x1_max_list, x2_max_list, x3_max_list, x4_max_list, x5_max_list, x10_max_list, x20_max_list, x99_max_list
"""Write into a csv."""
# check the length of them.
print(len(x0_target), len(x0_x10_max_list), len(x0_x20_max_list))   # 111 111 111
print(len(x1_target), len(x1_x10_max_list), len(x1_x20_max_list))   # 112 112 112
max_num = max(len(x0_target), len(x1_target), len(x2_target), len(x3_target), len(x4_target), len(x5_target), len(x6_target), len(x7_target), len(x8_target), len(x9_target),
              len(x10_target), len(x11_target), len(x12_target), len(x13_target), len(x14_target), len(x15_target), len(x16_target), len(x17_target), len(x18_target), len(x19_target))
print(max_num)
min_num = min(len(x0_target), len(x1_target), len(x2_target), len(x3_target), len(x4_target), len(x5_target), len(x6_target), len(x7_target), len(x8_target), len(x9_target),
              len(x10_target), len(x11_target), len(x12_target), len(x13_target), len(x14_target), len(x15_target), len(x16_target), len(x17_target), len(x18_target), len(x19_target))
print(min_num)
# OK

csv_path = r'./off_target_2plot.csv'

first_row = ['x0_target', 'x0_x10_max_list', 'x0_x20_max_list',
             'x1_target', 'x1_x10_max_list', 'x1_x20_max_list',
             'x2_target', 'x2_x10_max_list', 'x2_x20_max_list',
             'x3_target', 'x3_x10_max_list', 'x3_x20_max_list',
             'x4_target', 'x4_x10_max_list', 'x4_x20_max_list',
             'x5_target', 'x5_x10_max_list', 'x5_x20_max_list',
             'x6_target', 'x6_x10_max_list', 'x6_x20_max_list',
             'x7_target', 'x7_x10_max_list', 'x7_x20_max_list',
             'x8_target', 'x8_x10_max_list', 'x8_x20_max_list',
             'x9_target', 'x9_x10_max_list', 'x9_x20_max_list',
             'x10_target', 'x10_x10_max_list', 'x10_x20_max_list',
             'x11_target', 'x11_x10_max_list', 'x11_x20_max_list',
             'x12_target', 'x12_x10_max_list', 'x12_x20_max_list',
             'x13_target', 'x13_x10_max_list', 'x13_x20_max_list',
             'x14_target', 'x14_x10_max_list', 'x14_x20_max_list',
             'x15_target', 'x15_x10_max_list', 'x15_x20_max_list',
             'x16_target', 'x16_x10_max_list', 'x16_x20_max_list',
             'x17_target', 'x17_x10_max_list', 'x17_x20_max_list',
             'x18_target', 'x18_x10_max_list', 'x18_x20_max_list',
             'x19_target', 'x19_x10_max_list', 'x19_x20_max_list']
with open(csv_path, "a", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(first_row)
    for i in range(max_num):
        row = [x0_target[i], x0_x10_max_list[i], x0_x20_max_list[i],
               x1_target[i], x1_x10_max_list[i], x1_x20_max_list[i],
               x2_target[i], x2_x10_max_list[i], x2_x20_max_list[i],
               x3_target[i], x3_x10_max_list[i], x3_x20_max_list[i],
               x4_target[i], x4_x10_max_list[i], x4_x20_max_list[i],
               x5_target[i], x5_x10_max_list[i], x5_x20_max_list[i],
               x6_target[i], x6_x10_max_list[i], x6_x20_max_list[i],
               x7_target[i], x7_x10_max_list[i], x7_x20_max_list[i],
               x8_target[i], x8_x10_max_list[i], x8_x20_max_list[i],
               x9_target[i], x9_x10_max_list[i], x9_x20_max_list[i],
               x10_target[i], x10_x10_max_list[i], x10_x20_max_list[i],
               x11_target[i], x11_x10_max_list[i], x11_x20_max_list[i],
               x12_target[i], x12_x10_max_list[i], x12_x20_max_list[i],
               x13_target[i], x13_x10_max_list[i], x13_x20_max_list[i],
               x14_target[i], x14_x10_max_list[i], x14_x20_max_list[i],
               x15_target[i], x15_x10_max_list[i], x15_x20_max_list[i],
               x16_target[i], x16_x10_max_list[i], x16_x20_max_list[i],
               x17_target[i], x17_x10_max_list[i], x17_x20_max_list[i],
               x18_target[i], x18_x10_max_list[i], x18_x20_max_list[i],
               x19_target[i], x19_x10_max_list[i], x19_x20_max_list[i]
               ]
        writer.writerow(row)

