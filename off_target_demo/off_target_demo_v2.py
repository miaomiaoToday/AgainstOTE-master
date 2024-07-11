"""
Analyze the off-target effect.
Provider: Yechao Xu
"""
# import sys
# sys.path.append('..\evaluation')
# sys.path.append(r'..\utils')
import os
import torch
from tqdm.auto import tqdm
# from evaluation.docking import *
from docking_v2 import *
import csv
import re
from utils.reconstruct import reconstruct_from_generated_with_edges
import numpy as np
import datetime

"""
Obtain all target pocket pdb path, the fixed list is in out_ligand_receptor_list = [] and 
saved in out_ligand_receptor_list.csv
"""
protein_root = '../data/crossdocked_pocket10'
root_gen_dir = './outputs_2.4'
exp_dirs = os.listdir(root_gen_dir)

pattern = r'(?<=)\d+'
regex = re.compile(pattern)
exp_dirs.sort(key=lambda name: int(re.search(regex, name).group()))
# for a in exp_dirs:
#     print(a.split('_')[1])
# print(exp_dirs)

out_ligand_receptor_list = []
if not os.path.exists(r'./out_ligand_receptor_list.csv'):
    for exp_dir in exp_dirs:
        save_path = os.path.join(exp_dir, 'samples_all.pt')
        save_path = os.path.join(root_gen_dir, save_path)
        samples = torch.load(save_path, map_location='cpu')

        for i, data in enumerate(tqdm(samples.finished, desc='All')):
            ligand_name = data.ligand_filename
            receptor_name = os.path.join(
                os.path.dirname(data.ligand_filename),
                # os.path.basename(data.ligand_filename)[:10] + '.pdb'
                os.path.basename(data.ligand_filename)[:].split('.')[0] + '_pocket10.pdb'
            )
            # out_ligand_receptor_list.append((exp_dir, ligand_name, receptor_name))
            # print('Generated output: {}; Ligand: {}; Receptor: {}'.format(exp_dir, ligand_name, receptor_name))

            """here we want to read the centor position of target binding location."""
            # ligand_path = os.path.join('./data/crossdocked_pocket10', data.ligand_filename)
            ligand4pos = reconstruct_from_generated_with_edges(data)
            ligand_position = Chem.AddHs(ligand4pos, addCoords=True)
            try:
                not_converge = 10
                while not_converge > 0:
                    flag = UFFOptimizeMolecule(ligand_position)
                    not_converge = min(not_converge - 1, flag * 10)
            except RuntimeError:
                pass
            pos = ligand_position.GetConformer(0).GetPositions()
            center = (pos.max(0) + pos.min(0)) / 2
            x, y, z = center[0], center[1], center[2]
            # print(center, x, y, z)

            out_ligand_receptor_list.append((exp_dir, ligand_name, receptor_name, x, y, z))
            break
            # vina_task = QVinaDockingTask.from_generated_data(data, protein_root=protein_root)

    # generate a fixed gallery list of pockets, [(exp_dir, ligand_name, receptor_name), (), ...]
    print(out_ligand_receptor_list)
    # write into a csv, saving time
    csv_path = r'./out_ligand_receptor_list.csv'
    with open(csv_path, "w", newline='') as f:
        writer = csv.writer(f)
        for row in out_ligand_receptor_list:
            writer.writerow(row)
else:
    with open(r'./out_ligand_receptor_list.csv') as f:
        reader = csv.reader(f)
        # headers = next(reader)  # this will jump the 1st row
        for row in reader:
            # print(row)
            out_ligand_receptor_list.append((row[0], row[1], row[2], row[3], row[4], row[5]))


"""
Calculate the qvina2 score of current pdb, and other pdbs
"""
csv_path = r'./vina_scores_2.4_formal.csv'
first_row = ['exp_id', 'sdf_id', 'pdb_id', 'exp_dir', 'ligand', 'receptor', 'original pair', 'Best affinity', 'RMSD', 'data_i', 'pock_id']
with open(csv_path, "a", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(first_row)
    f.close()
exp_index = -1
# exp_index = 6
total_count = 0
for exp_dir in exp_dirs:
    exp_index += 1
    # if exp_index != 7:
        # continue
    # if exp_dir == 'sample_39_2023_03_17__02_10_07':
    #     print('ok')
    # else:
    #     continue
    save_path = os.path.join(exp_dir, 'samples_all.pt')
    save_path = os.path.join(root_gen_dir, save_path)
    samples = torch.load(save_path, map_location='cpu')

    for i, data in enumerate(tqdm(samples.finished, desc='All')):
        smiles = data.smiles
        # here we want to do 100 tasks, including 99 other ones and 1 matched one.
        for pock_id in range(0, 10):
            total_count += 1
            
            # if total_count <= 266711:
            	# continue
            # if total_count > 300000:
            # 	break
            now = datetime.datetime.now()
            print('Time:{:0>4d}-{:0>2d}-{:0>2d}, {:0>2d}: {:0>2d}: {:0>2d},  {} / 100,000 is ok.'.format(now.year,
                                                                                                now.month,
                                                                                                now.day,
                                                                                                now.hour,
                                                                                                now.minute,
                                                                                                now.second,
                                                                                                total_count))
            pair_idx = 0
             # here ! when generating task, also take into the ligand data matched with protein_root, for correct position at pocket.
            """ ---------- """
            out_ligand_receptor = out_ligand_receptor_list[pock_id]
            cexp_dir = out_ligand_receptor[0]
            ligand = out_ligand_receptor[1]
            receptor = out_ligand_receptor[2]
            center = [out_ligand_receptor[3], out_ligand_receptor[4], out_ligand_receptor[5]]
            center = np.array(center, dtype=np.float64)
            """ -----just need to feed center position-----"""
            vina_task = QVinaDockingTask.from_generated_data_mod(data, protein_root=protein_root, pock_id=pock_id, out_ligand_receptor_list=out_ligand_receptor_list, center=center)

            # out_ligand_receptor = out_ligand_receptor_list[pock_id]
            # cexp_dir = out_ligand_receptor[0]
            # ligand = out_ligand_receptor[1]
            # receptor = out_ligand_receptor[2]
            # where is a original pair of ligand and receptor
            # print('Here is the original pair {} and {} is matched'.format(ligand, data.ligand_filename))
            if ligand == data.ligand_filename:
                pair_idx = 1
            if pair_idx == 1:
                print('Here is the original pair {} and {} is matched'.format(ligand, data.ligand_filename))

            vina_score = vina_task.run_sync()
            try:
            	affinity = vina_score[0]['affinity']
            	rmsd_ref = vina_score[0]['rmsd_ref']
            except:
            	affinity = np.nan
            	rmsd_ref = np.nan
            	print('affinity nan.')
            # vina_score = 10
            # save csv: ligand, receptor, original pair, affinity, rmsd_ref, data_i, pock_id

            with open(csv_path, "a", newline='') as f:
                writer = csv.writer(f)
                row = [exp_index, i, pock_id, exp_dir, ligand, receptor, pair_idx, affinity, rmsd_ref, i, pock_id]
                writer.writerow(row)

        # ligand_name = data.ligand_filename
        # receptor_name = os.path.join(
        #     os.path.dirname(data.ligand_filename),
        #     # os.path.basename(data.ligand_filename)[:10] + '.pdb'
        #     os.path.basename(data.ligand_filename)[:].split('.')[0] + '_pocket10.pdb'
        # )
        # out_ligand_receptor_list.append((exp_dir, ligand_name, receptor_name))

print('finished')
