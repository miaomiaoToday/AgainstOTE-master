"""
Analyze the generated docking scores for visualization
"""
import csv
import numpy as np
import math


def sigmoid_reverse(x):
    return 1 / (1 + math.exp(x))


def pocket_each_ana(x0_list, k=5):
    x0_target = [[], [], []]
    for data_item in x0_list:
        if data_item[3] == 1:
            x0_target[0].append(data_item[0])
            x0_target[1].append(data_item[1])
            x0_target[2].append(data_item[4])

    if len(x0_target[2]) < k:
        k = min(len(x0_target[2]), k) - 1

    # print(x0_target)

    # take only top k target scores and other scores
    max_five_indices = np.argpartition(x0_target[2], k)[:k]
    x0_target_top5 = [x0_target[2][index] for index in max_five_indices]
    # print('top k target scores ')

    x0_others = []
    for index in max_five_indices:
        exp_id = x0_target[0][index]
        mol_id = x0_target[1][index]
        # get the corresponding other scores.
        for data_item in x0_list:
            if data_item[0] == exp_id and data_item[1] == mol_id and data_item[3] != 1:
                x0_others.append(data_item[4])


    # """7. Top 10 mean, ..."""
    # """7. All other mean, ..."""
    # temp_list = []
    # x_other_list = []
    # for data_item in x0_list:
    #     if data_item[0] == data_item[1] == data_item[2] == 0:
    #         # temp_list.append(data_item[4])
    #         pass
    #     elif 0 < data_item[2] < 100:
    #         if data_item[3] == 1:
    #             continue
    #         temp_list.append(data_item[4])
    #     elif data_item[2] == 0:
    #         if len(temp_list) == 0 and data_item[3] == 0:
    #             temp_list.append(data_item[4])
    #             continue
    #         if len(temp_list) == 0 and data_item[3] == 1:
    #             continue
    #         # temp_list = sorted(temp_list)
    #         # max_val = (temp_list[0] + temp_list[1] + temp_list[2] + temp_list[3] + temp_list[4] + temp_list[5] + temp_list[6] + temp_list[7] + temp_list[8]) / 9
    #         # x_other_list = x_other_list + temp_list
    #         # for ele in temp_list:
    #
    #         x_other_list.append(temp_list)
    #         temp_list = []
    #         # temp_list.append(data_item[4])
    #
    # # temp_list = sorted(temp_list)
    # # max_val = (temp_list[0] + temp_list[1] + temp_list[2] + temp_list[3] + temp_list[4] + temp_list[5] + temp_list[6] + temp_list[7] + temp_list[8]) / 9
    # # x_other_list = x_other_list + temp_list
    # x_other_list.append(temp_list)
    # print(x_other_list)
    #
    # # """8. Top 20 mean, ..."""
    # # temp_list = []
    # # x20_max_list = []
    # # for data_item in x0_list:
    # #     if data_item[0] == data_item[1] == data_item[2] == 0:
    # #         # temp_list.append(data_item[4])
    # #         pass
    # #     elif 0 < data_item[2] < 100:
    # #         if data_item[3] == 1:
    # #             continue
    # #         temp_list.append(data_item[4])
    # #     elif data_item[2] == 0:
    # #         if len(temp_list) == 0 and data_item[3] == 0:
    # #             temp_list.append(data_item[4])
    # #             continue
    # #         if len(temp_list) == 0 and data_item[3] == 1:
    # #             continue
    # #         temp_list = sorted(temp_list)
    # #         max_val = 0
    # #         for i in range(20):
    # #             max_val += temp_list[i]
    # #         max_val = max_val / 20
    # #         # max_val = (temp_list[0] + temp_list[1] + temp_list[2] + temp_list[3] + temp_list[4] + temp_list[5] + temp_list[6] + temp_list[7] + temp_list[8] + temp_list[9]) / 10
    # #         x20_max_list.append(max_val)
    # #         temp_list = []
    # #         # temp_list.append(data_item[4])
    # #
    # # temp_list = sorted(temp_list)
    # # max_val = 0
    # # for i in range(20):
    # #     max_val += temp_list[i]
    # # max_val = max_val / 20
    # # x20_max_list.append(max_val)
    # # print(x20_max_list)

    return x0_target_top5, x0_others


data_list = []
with open(r'./vina_scores_2.4_formal.csv') as f:
    reader = csv.reader(f)
    headers = next(reader)  # this will jump the 1st row
    for row in reader:
        # source pocket id, generated molecular id, target pocket id, if the pair, score
        data_list.append((row[0], row[1], row[2], row[6], row[7]))
    # print(data_list)
    # ('26', '10', '25', '0', '-9.2') is ele in data_list

# we take 20 source pocket, to generated 100 +- mols for each of them
# we take 10 source pocket, to generated 100 +- mols for each of them
exp_data_list = []
for data_item in data_list:
    if int(data_item[0]) < 10:
        exp_data_list.append((int(data_item[0]), int(data_item[1]), int(data_item[2]),
                              int(data_item[3]), float(data_item[4])))
# exp_data_list, simulate the animal env.

# split as mols of 20 pockets
# split as mols of 10 pockets
mols_pock = []
x0_list, x1_list, x2_list, x3_list, x4_list, x5_list, x6_list, x7_list, x8_list, x9_list = [], [], [], [], [], [], [], [], [], []  # pocket0_target
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

k = 5

x0_target, x0_x10_other_list = pocket_each_ana(x0_list, k=k)
x1_target, x1_x10_other_list = pocket_each_ana(x1_list, k=k)
x2_target, x2_x10_other_list = pocket_each_ana(x2_list, k=k)
x3_target, x3_x10_other_list = pocket_each_ana(x3_list, k=k)
x4_target, x4_x10_other_list = pocket_each_ana(x4_list, k=k)
x5_target, x5_x10_other_list = pocket_each_ana(x5_list, k=k)
x6_target, x6_x10_other_list = pocket_each_ana(x6_list, k=k)
x7_target, x7_x10_other_list = pocket_each_ana(x7_list, k=k)
x8_target, x8_x10_other_list = pocket_each_ana(x8_list, k=k)
x9_target, x9_x10_other_list = pocket_each_ana(x9_list, k=k)

# delete nan
x0_target = [x for x in x0_target if not isinstance(x, float) or not math.isnan(x)]
x0_x10_other_list = [x for x in x0_x10_other_list if not isinstance(x, float) or not math.isnan(x)]
x1_target = [x for x in x1_target if not isinstance(x, float) or not math.isnan(x)]
x1_x10_other_list = [x for x in x1_x10_other_list if not isinstance(x, float) or not math.isnan(x)]
x2_target = [x for x in x2_target if not isinstance(x, float) or not math.isnan(x)]
x2_x10_other_list = [x for x in x2_x10_other_list if not isinstance(x, float) or not math.isnan(x)]
x3_target = [x for x in x3_target if not isinstance(x, float) or not math.isnan(x)]
x3_x10_other_list = [x for x in x3_x10_other_list if not isinstance(x, float) or not math.isnan(x)]
x4_target = [x for x in x4_target if not isinstance(x, float) or not math.isnan(x)]
x4_x10_other_list = [x for x in x4_x10_other_list if not isinstance(x, float) or not math.isnan(x)]
x5_target = [x for x in x5_target if not isinstance(x, float) or not math.isnan(x)]
x5_x10_other_list = [x for x in x5_x10_other_list if not isinstance(x, float) or not math.isnan(x)]
x6_target = [x for x in x6_target if not isinstance(x, float) or not math.isnan(x)]
x6_x10_other_list = [x for x in x6_x10_other_list if not isinstance(x, float) or not math.isnan(x)]
x7_target = [x for x in x7_target if not isinstance(x, float) or not math.isnan(x)]
x7_x10_other_list = [x for x in x7_x10_other_list if not isinstance(x, float) or not math.isnan(x)]
x8_target = [x for x in x8_target if not isinstance(x, float) or not math.isnan(x)]
x8_x10_other_list = [x for x in x8_x10_other_list if not isinstance(x, float) or not math.isnan(x)]
x9_target = [x for x in x9_target if not isinstance(x, float) or not math.isnan(x)]
x9_x10_other_list = [x for x in x9_x10_other_list if not isinstance(x, float) or not math.isnan(x)]

# do @ top K, K = 5
# do top k at other list, not all other scores.
k = 5
x0_x10_other_list = sorted(x0_x10_other_list)
x0_x10_other_list = x0_x10_other_list[:k]
x1_x10_other_list = sorted(x1_x10_other_list)
x1_x10_other_list = x1_x10_other_list[:k]
x2_x10_other_list = sorted(x2_x10_other_list)
x2_x10_other_list = x2_x10_other_list[:k]
x3_x10_other_list = sorted(x3_x10_other_list)
x3_x10_other_list = x3_x10_other_list[:k]
x4_x10_other_list = sorted(x4_x10_other_list)
x4_x10_other_list = x4_x10_other_list[:k]
x5_x10_other_list = sorted(x5_x10_other_list)
x5_x10_other_list = x5_x10_other_list[:k]
x6_x10_other_list = sorted(x6_x10_other_list)
x6_x10_other_list = x6_x10_other_list[:k]
x7_x10_other_list = sorted(x7_x10_other_list)
x7_x10_other_list = x7_x10_other_list[:k]
x8_x10_other_list = sorted(x8_x10_other_list)
x8_x10_other_list = x8_x10_other_list[:k]
x9_x10_other_list = sorted(x9_x10_other_list)
x9_x10_other_list = x9_x10_other_list[:k]


x0_target_mean, x0_target_std, x0_other_mean, x0_other_std = np.mean(x0_target), np.std(x0_target), np.mean(x0_x10_other_list), np.std(x0_x10_other_list)
x1_target_mean, x1_target_std, x1_other_mean, x1_other_std = np.mean(x1_target), np.std(x1_target), np.mean(x1_x10_other_list), np.std(x1_x10_other_list)
x2_target_mean, x2_target_std, x2_other_mean, x2_other_std = np.mean(x2_target), np.std(x2_target), np.mean(x2_x10_other_list), np.std(x2_x10_other_list)
x3_target_mean, x3_target_std, x3_other_mean, x3_other_std = np.mean(x3_target), np.std(x3_target), np.mean(x3_x10_other_list), np.std(x3_x10_other_list)
x4_target_mean, x4_target_std, x4_other_mean, x4_other_std = np.mean(x4_target), np.std(x4_target), np.mean(x4_x10_other_list), np.std(x4_x10_other_list)
x5_target_mean, x5_target_std, x5_other_mean, x5_other_std = np.mean(x5_target), np.std(x5_target), np.mean(x5_x10_other_list), np.std(x5_x10_other_list)
x6_target_mean, x6_target_std, x6_other_mean, x6_other_std = np.mean(x6_target), np.std(x6_target), np.mean(x6_x10_other_list), np.std(x6_x10_other_list)
x7_target_mean, x7_target_std, x7_other_mean, x7_other_std = np.mean(x7_target), np.std(x7_target), np.mean(x7_x10_other_list), np.std(x7_x10_other_list)
x8_target_mean, x8_target_std, x8_other_mean, x8_other_std = np.mean(x8_target), np.std(x8_target), np.mean(x8_x10_other_list), np.std(x8_x10_other_list)
x9_target_mean, x9_target_std, x9_other_mean, x9_other_std = np.mean(x9_target), np.std(x9_target), np.mean(x9_x10_other_list), np.std(x9_x10_other_list)

# TSR0 = sigmoid_reverse((x0_target_mean - x0_other_mean)/x0_other_std)
# TSR1 = sigmoid_reverse((x1_target_mean - x1_other_mean)/x1_other_std)
# TSR2 = sigmoid_reverse((x2_target_mean - x2_other_mean)/x2_other_std)
# TSR3 = sigmoid_reverse((x3_target_mean - x3_other_mean)/x3_other_std)
# TSR4 = sigmoid_reverse((x4_target_mean - x4_other_mean)/x4_other_std)
# TSR5 = sigmoid_reverse((x5_target_mean - x5_other_mean)/x5_other_std)
# TSR6 = sigmoid_reverse((x6_target_mean - x6_other_mean)/x6_other_std)
# TSR7 = sigmoid_reverse((x7_target_mean - x7_other_mean)/x7_other_std)
# TSR8 = sigmoid_reverse((x8_target_mean - x8_other_mean)/x8_other_std)
# TSR9 = sigmoid_reverse((x9_target_mean - x9_other_mean)/x9_other_std)

TSR0 = sigmoid_reverse((x0_target_mean - x0_other_mean)/1)
TSR1 = sigmoid_reverse((x1_target_mean - x1_other_mean)/1)
TSR2 = sigmoid_reverse((x2_target_mean - x2_other_mean)/1)
TSR3 = sigmoid_reverse((x3_target_mean - x3_other_mean)/1)
TSR4 = sigmoid_reverse((x4_target_mean - x4_other_mean)/1)
TSR5 = sigmoid_reverse((x5_target_mean - x5_other_mean)/1)
TSR6 = sigmoid_reverse((x6_target_mean - x6_other_mean)/1)
TSR7 = sigmoid_reverse((x7_target_mean - x7_other_mean)/1)
TSR8 = sigmoid_reverse((x8_target_mean - x8_other_mean)/1)
TSR9 = sigmoid_reverse((x9_target_mean - x9_other_mean)/1)

total_target_mean = (x0_target_mean + x1_target_mean + x2_target_mean + x3_target_mean + x4_target_mean + x5_target_mean + x6_target_mean + x7_target_mean + x8_target_mean + x9_target_mean)/10
total_other_mean = (x0_other_mean + x1_other_mean + x2_other_mean + x3_other_mean + x4_other_mean + x5_other_mean + x6_other_mean + x7_other_mean + x8_other_mean + x9_other_mean)/10
total_TSR_mean = (TSR0 + TSR1 + TSR2 + TSR3 + TSR4 + TSR5 + TSR6 + TSR7 + TSR8 + TSR9)/10
print(TSR0, TSR1, TSR2, TSR3, TSR4, TSR5, TSR6, TSR7, TSR8, TSR9)
print(total_TSR_mean)
# total_TSR_mean = 0
"""summary 1-8. 8 list"""
# x0_list, x1_max_list, x2_max_list, x3_max_list, x4_max_list, x5_max_list, x10_max_list, x20_max_list, x99_max_list
"""Write into a csv."""
# check the length of them.

# OK

csv_path = r'./eval_off_target_2.4_formal.csv'

first_row = ['x0_target_mean', 'x0_target_std', 'x0_other_mean', 'x0_other_std', 'TSR0',
             'x1_target_mean', 'x1_target_std', 'x1_other_mean', 'x1_other_std', 'TSR1',
             'x2_target_mean', 'x2_target_std', 'x2_other_mean', 'x2_other_std', 'TSR2',
             'x3_target_mean', 'x3_target_std', 'x3_other_mean', 'x3_other_std', 'TSR3',
             'x4_target_mean', 'x4_target_std', 'x4_other_mean', 'x4_other_std', 'TSR4',
             'x5_target_mean', 'x5_target_std', 'x5_other_mean', 'x5_other_std', 'TSR5',
             'x6_target_mean', 'x6_target_std', 'x6_other_mean', 'x6_other_std', 'TSR6',
             'x7_target_mean', 'x7_target_std', 'x7_other_mean', 'x7_other_std', 'TSR7',
             'x8_target_mean', 'x8_target_std', 'x8_other_mean', 'x8_other_std', 'TSR8',
             'x9_target_mean', 'x9_target_std', 'x9_other_mean', 'x9_other_std', 'TSR9',
             'total_target_mean', 'total_other_mean', 'total_TSR_mean']

with open(csv_path, "a", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(first_row)
    # for i in range(10):
    row = [x0_target_mean, x0_target_std, x0_other_mean, x0_other_std, TSR0,
                x1_target_mean, x1_target_std, x1_other_mean, x1_other_std, TSR1,
                x2_target_mean, x2_target_std, x2_other_mean, x2_other_std, TSR2,
                x3_target_mean, x3_target_std, x3_other_mean, x3_other_std, TSR3,
                x4_target_mean, x4_target_std, x4_other_mean, x4_other_std, TSR4,
                x5_target_mean, x5_target_std, x5_other_mean, x5_other_std, TSR5,
                x6_target_mean, x6_target_std, x6_other_mean, x6_other_std, TSR6,
                x7_target_mean, x7_target_std, x7_other_mean, x7_other_std, TSR7,
                x8_target_mean, x8_target_std, x8_other_mean, x8_other_std, TSR8,
                x9_target_mean, x9_target_std, x9_other_mean, x9_other_std, TSR9,
           total_target_mean, total_other_mean, total_TSR_mean
           ]
    writer.writerow(row)

# first_row = ['x0_target_mean', 'x0_other_mean', 'x0_other_std', 
#              'x1_target_mean', 'x1_other_mean', 'x1_other_std', 
#              'x2_target_mean', 'x2_other_mean', 'x2_other_std', 
#              'x3_target_mean', 'x3_other_mean', 'x3_other_std', 
#              'x4_target_mean', 'x4_other_mean', 'x4_other_std', 
#              'x5_target_mean', 'x5_other_mean', 'x5_other_std', 
#              'x6_target_mean', 'x6_other_mean', 'x6_other_std', 
#              'x7_target_mean', 'x7_other_mean', 'x7_other_std', 
#              'x8_target_mean', 'x8_other_mean', 'x8_other_std', 
#              'x9_target_mean', 'x9_other_mean', 'x9_other_std', 
#              'total_target_mean', 'total_other_mean']

# with open(csv_path, "a", newline='') as f:
#     writer = csv.writer(f)
#     writer.writerow(first_row)
#     # for i in range(10):
#     row = [x0_target_mean, x0_other_mean, x0_other_std, 
#                 x1_target_mean, x1_other_mean, x1_other_std, 
#                 x2_target_mean, x2_other_mean, x2_other_std, 
#                 x3_target_mean, x3_other_mean, x3_other_std, 
#                 x4_target_mean, x4_other_mean, x4_other_std, 
#                 x5_target_mean, x5_other_mean, x5_other_std, 
#                 x6_target_mean, x6_other_mean, x6_other_std, 
#                 x7_target_mean, x7_other_mean, x7_other_std, 
#                 x8_target_mean, x8_other_mean, x8_other_std, 
#                 x9_target_mean, x9_other_mean, x9_other_std, 
#            total_target_mean, total_other_mean
#            ]
#     writer.writerow(row)
    