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
    # print(x0_target)

    """7. Top 10 mean, ..."""
    temp_list = []
    x10_max_list = []
    counts = []
    index = 0
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
            count_off = 0
            for score in temp_list:
                if score <= x0_target[index]:
                    count_off += 1
            counts.append(count_off)
            temp_list = []
            index += 1
            # temp_list.append(data_item[4])


    temp_list = sorted(temp_list)
    count_off = 0
    for score in temp_list:
        if score <= x0_target[index]:
            count_off += 1
    counts.append(count_off)
    print(counts)

    return counts


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

x0_counts = pocket_each_ana(x0_list)
x1_counts = pocket_each_ana(x1_list)
x2_counts = pocket_each_ana(x2_list)
x3_counts = pocket_each_ana(x3_list)
x4_counts = pocket_each_ana(x4_list)
x5_counts = pocket_each_ana(x5_list)
x6_counts = pocket_each_ana(x6_list)
x7_counts = pocket_each_ana(x7_list)
x8_counts = pocket_each_ana(x8_list)
x9_counts = pocket_each_ana(x9_list)
x10_counts = pocket_each_ana(x10_list)
x11_counts = pocket_each_ana(x11_list)
x12_counts = pocket_each_ana(x12_list)
x13_counts = pocket_each_ana(x13_list)
x14_counts = pocket_each_ana(x14_list)
x15_counts = pocket_each_ana(x15_list)
x16_counts = pocket_each_ana(x16_list)
x17_counts = pocket_each_ana(x17_list)
x18_counts = pocket_each_ana(x18_list)
x19_counts = pocket_each_ana(x19_list)


"""summary 1-8. 8 list"""
# x0_list, x1_max_list, x2_max_list, x3_max_list, x4_max_list, x5_max_list, x10_max_list, x20_max_list, x99_max_list
"""Write into a csv."""
# check the length of them.
print(len(x0_counts))   # 111 111 111
print(len(x0_counts))   # 112 112 112
max_num = max(len(x0_counts), len(x1_counts), len(x2_counts), len(x3_counts), len(x4_counts), len(x5_counts), len(x6_counts), len(x7_counts), len(x8_counts), len(x9_counts),
              len(x10_counts), len(x11_counts), len(x12_counts), len(x13_counts), len(x14_counts), len(x15_counts), len(x16_counts), len(x17_counts), len(x18_counts), len(x19_counts))
print(max_num)
min_num = min(len(x0_counts), len(x1_counts), len(x2_counts), len(x3_counts), len(x4_counts), len(x5_counts), len(x6_counts), len(x7_counts), len(x8_counts), len(x9_counts),
              len(x10_counts), len(x11_counts), len(x12_counts), len(x13_counts), len(x14_counts), len(x15_counts), len(x16_counts), len(x17_counts), len(x18_counts), len(x19_counts))
print(min_num)
# OK
x0_counts_m = np.mean(x0_counts)
x1_counts_m = np.mean(x1_counts)
x2_counts_m = np.mean(x2_counts)
x3_counts_m = np.mean(x3_counts)
x4_counts_m = np.mean(x4_counts)
x5_counts_m = np.mean(x5_counts)
x6_counts_m = np.mean(x6_counts)
x7_counts_m = np.mean(x7_counts)
x8_counts_m = np.mean(x8_counts)
x9_counts_m = np.mean(x9_counts)
x10_counts_m = np.mean(x10_counts)
x11_counts_m = np.mean(x11_counts)
x12_counts_m = np.mean(x12_counts)
x13_counts_m = np.mean(x13_counts)
x14_counts_m = np.mean(x14_counts)
x15_counts_m = np.mean(x15_counts)
x16_counts_m = np.mean(x16_counts)
x17_counts_m = np.mean(x17_counts)
x18_counts_m = np.mean(x18_counts)
x19_counts_m = np.mean(x19_counts)
print(x0_counts_m, x1_counts_m, x2_counts_m, x3_counts_m, x4_counts_m, x5_counts_m, x6_counts_m, x7_counts_m, x8_counts_m, x9_counts_m,
x10_counts_m, x11_counts_m, x12_counts_m, x13_counts_m, x14_counts_m, x15_counts_m, x16_counts_m, x17_counts_m, x18_counts_m, x19_counts_m)


csv_path = r'./off_target_3plot.csv'

first_row = ['x0_counts_m', 'x1_counts_m', 'x2_counts_m',
             'x3_counts_m', 'x4_counts_m', 'x5_counts_m',
             'x6_counts_m', 'x7_counts_m', 'x8_counts_m',
             'x9_counts_m', 'x10_counts_m', 'x11_counts_m',
             'x12_counts_m', 'x13_counts_m', 'x14_counts_m',
             'x15_counts_m', 'x16_counts_m', 'x17_counts_m',
             'x18_counts_m', 'x19_counts_m']
with open(csv_path, "a", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(first_row)
    # for i in range(max_num):
    row = [x0_counts_m, x1_counts_m, x2_counts_m, x3_counts_m, x4_counts_m, x5_counts_m, x6_counts_m, x7_counts_m, x8_counts_m, x9_counts_m,
x10_counts_m, x11_counts_m, x12_counts_m, x13_counts_m, x14_counts_m, x15_counts_m, x16_counts_m, x17_counts_m, x18_counts_m, x19_counts_m
           ]
    writer.writerow(row)

