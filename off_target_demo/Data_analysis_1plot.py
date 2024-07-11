"""
Analyze the generated docking scores for visualization
"""
import csv
import numpy as np


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

"""1. Mol-Target, 2000+ Scores"""
x0_list = []
for data_item in exp_data_list:
    if data_item[3] == 1:
        # print(data_item[4])
        x0_list.append(data_item[4])
print(x0_list)

"""2. Top 1, Mol-OtherTarget, For each mol with 100 pockets, 2000+ Scores"""
temp_list = []
x1_max_list = []
for data_item in exp_data_list:
    if data_item[0] == data_item[1] == data_item[2] == 0:
        # temp_list.append(data_item[4])
        pass
    elif 0 < data_item[2] < 100:
        if data_item[3] == 1:
            continue
        temp_list.append(data_item[4])
    elif data_item[2] == 0:
        max_val = min(temp_list)
        x1_max_list.append(max_val)
        temp_list = []
        # temp_list.append(data_item[4])
    else:
        print('error: out of 100 pockets')

max_val = min(temp_list)    # the last one max mol-other_pocket docking score
x1_max_list.append(max_val)
print(x1_max_list)

"""3. Top 2 mean, """
temp_list = []
x2_max_list = []
for data_item in exp_data_list:
    if data_item[0] == data_item[1] == data_item[2] == 0:
        # temp_list.append(data_item[4])
        pass
    elif 0 < data_item[2] < 100:
        if data_item[3] == 1:
            continue
        temp_list.append(data_item[4])
    elif data_item[2] == 0:
        temp_list = sorted(temp_list)
        max_val = (temp_list[0] + temp_list[1]) / 2
        x2_max_list.append(max_val)
        temp_list = []
        # temp_list.append(data_item[4])

temp_list = sorted(temp_list)
max_val = (temp_list[0] + temp_list[1]) / 2
x2_max_list.append(max_val)
print(x2_max_list)

"""4.  Top 3 mean, ..."""
temp_list = []
x3_max_list = []
for data_item in exp_data_list:
    if data_item[0] == data_item[1] == data_item[2] == 0:
        # temp_list.append(data_item[4])
        pass
    elif 0 < data_item[2] < 100:
        if data_item[3] == 1:
            continue
        temp_list.append(data_item[4])
    elif data_item[2] == 0:
        temp_list = sorted(temp_list)
        max_val = (temp_list[0] + temp_list[1] + temp_list[2]) / 3
        x3_max_list.append(max_val)
        temp_list = []
        # temp_list.append(data_item[4])

temp_list = sorted(temp_list)
max_val = (temp_list[0] + temp_list[1] + temp_list[2]) / 3
x3_max_list.append(max_val)
print(x3_max_list)

"""5. Top 4 mean, ..."""
temp_list = []
x4_max_list = []
for data_item in exp_data_list:
    if data_item[0] == data_item[1] == data_item[2] == 0:
        # temp_list.append(data_item[4])
        pass
    elif 0 < data_item[2] < 100:
        if data_item[3] == 1:
            continue
        temp_list.append(data_item[4])
    elif data_item[2] == 0:
        temp_list = sorted(temp_list)
        max_val = (temp_list[0] + temp_list[1] + temp_list[2] + temp_list[3]) / 4
        x4_max_list.append(max_val)
        temp_list = []
        # temp_list.append(data_item[4])

temp_list = sorted(temp_list)
max_val = (temp_list[0] + temp_list[1] + temp_list[2] + temp_list[3]) / 4
x4_max_list.append(max_val)
print(x4_max_list)

"""6. Top 5 mean, ..."""
temp_list = []
x5_max_list = []
for data_item in exp_data_list:
    if data_item[0] == data_item[1] == data_item[2] == 0:
        # temp_list.append(data_item[4])
        pass
    elif 0 < data_item[2] < 100:
        if data_item[3] == 1:
            continue
        temp_list.append(data_item[4])
    elif data_item[2] == 0:
        temp_list = sorted(temp_list)
        max_val = (temp_list[0] + temp_list[1] + temp_list[2] + temp_list[3] + temp_list[4]) / 5
        x5_max_list.append(max_val)
        temp_list = []
        # temp_list.append(data_item[4])

temp_list = sorted(temp_list)
max_val = (temp_list[0] + temp_list[1] + temp_list[2] + temp_list[3] + temp_list[4]) / 5
x5_max_list.append(max_val)
print(x5_max_list)

"""7. Top 10 mean, ..."""
temp_list = []
x10_max_list = []
for data_item in exp_data_list:
    if data_item[0] == data_item[1] == data_item[2] == 0:
        # temp_list.append(data_item[4])
        pass
    elif 0 < data_item[2] < 100:
        if data_item[3] == 1:
            continue
        temp_list.append(data_item[4])
    elif data_item[2] == 0:
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
for data_item in exp_data_list:
    if data_item[0] == data_item[1] == data_item[2] == 0:
        # temp_list.append(data_item[4])
        pass
    elif 0 < data_item[2] < 100:
        if data_item[3] == 1:
            continue
        temp_list.append(data_item[4])
    elif data_item[2] == 0:
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

"""9. Top 99 mean, ..."""
temp_list = []
x99_max_list = []
for data_item in exp_data_list:
    if data_item[0] == data_item[1] == data_item[2] == 0:
        # temp_list.append(data_item[4])
        pass
    elif 0 < data_item[2] < 100:
        if data_item[3] == 1:
            continue
        temp_list.append(data_item[4])
    elif data_item[2] == 0:
        temp_list = sorted(temp_list)
        max_val = 0
        for i in range(len(temp_list)):
            max_val += temp_list[i]
        max_val = max_val / len(temp_list)
        # max_val = (temp_list[0] + temp_list[1] + temp_list[2] + temp_list[3] + temp_list[4] + temp_list[5] + temp_list[6] + temp_list[7] + temp_list[8] + temp_list[9]) / 10
        x99_max_list.append(max_val)
        temp_list = []
        # temp_list.append(data_item[4])

temp_list = sorted(temp_list)
max_val = 0
for i in range(len(temp_list)):
    max_val += temp_list[i]
max_val = max_val / len(temp_list)
x99_max_list.append(max_val)
print(x99_max_list)

"""summary 1-8. 8 list"""
# x0_list, x1_max_list, x2_max_list, x3_max_list, x4_max_list, x5_max_list, x10_max_list, x20_max_list, x99_max_list
"""Write into a csv."""
# check the length of them.
print(len(x0_list), len(x1_max_list), len(x2_max_list), len(x3_max_list), len(x4_max_list), len(x5_max_list), len(x10_max_list), len(x20_max_list), len(x99_max_list))
# OK

csv_path = r'./off_target_1plot.csv'

first_row = ['x0_list', 'x1_max_list', 'x2_max_list', 'x3_max_list', 'x4_max_list', 'x5_max_list', 'x10_max_list', 'x20_max_list', 'x99_max_list']
with open(csv_path, "a", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(first_row)
    for i in range(len(x0_list)):
        row = [x0_list[i], x1_max_list[i], x2_max_list[i], x3_max_list[i], x4_max_list[i], x5_max_list[i], x10_max_list[i], x20_max_list[i], x99_max_list[i]]
        writer.writerow(row)

