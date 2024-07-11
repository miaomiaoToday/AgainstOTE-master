import os
import argparse
from copy import deepcopy
import torch
from tqdm.auto import tqdm
from rdkit.Chem.QED import qed
import sys

sys.path.append('.')
import numpy as np
from utils.reconstruct import reconstruct_from_generated_with_edges
from sascorer import compute_sa_score
from evaluation.docking import *
from utils.misc import *
from evaluation.scoring_func import *
# for calculate MW and other properties for QED
from rdkit.Chem import rdMolDescriptors as rdmd  # for MW
from rdkit.Chem import MolSurf  # for TPSA
import csv

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp_name', type=str)
    parser.add_argument('--result_root', type=str, default='./outputs')
    parser.add_argument('--protein_root', type=str, default='../data/crossdocked_pocket10')
    args = parser.parse_args()

    ours = r'xxx\sample_12_2023_11_26__03_04_48'

    exp_dir = ours
    save_path = os.path.join(exp_dir, 'samples_all.pt')  # get_saved_filename(exp_dir))

    logger = get_logger('evaluate', exp_dir)
    logger.info(args)
    logger.info(save_path)

    samples = torch.load(save_path, map_location='cpu')

    sim_with_train = SimilarityWithTrain()
    results = []

    # write down the chemical performance for plotting
    csv_path = r'.\molpropertiesMW_Ours7B.csv'
    first_row = ['mols_root', 'hba', 'hbd', 'MW', 'TPSA']
    with open(csv_path, "a", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(first_row)
        f.close()

    for i, data in enumerate(tqdm(samples.finished, desc='All')):
        try:
            temp_res = []
            mol = reconstruct_from_generated_with_edges(data)
            _, _, _, hba, hbd = get_chem(mol)
            MW = rdmd._CalcMolWt(mol)
            TPSA = MolSurf.TPSA(mol)

            # record chemical performance
            with open(csv_path, "a", newline='') as f:
                writer = csv.writer(f)
                row = [exp_dir, hba, hbd, MW, TPSA]
                writer.writerow(row)

            results.append({
                'smiles': data.smiles,
                'qed': qed(mol),
                'sa': compute_sa_score(mol),
                'lipinski': obey_lipinski(mol),
                'logp': get_logp(mol),
                'hba': hba,
                'hbd': hbd,
                'MW': MW,
                'TPSA': TPSA,
            })
        except Exception as e:
            logger.warning('Failed %d' % i)
            logger.warning(e)

        logger.info('Number of results: %d' % len(results))

    hba_list, hbd_list, MW_list, TPSA_list = [], [], [], []
    for result in results:
        hba_list.append(result['hba'])
        hbd_list.append(result['hbd'])
        MW_list.append(result['MW'])
        TPSA_list.append(result['TPSA'])
    hba_mean, hba_std = np.mean(hba_list), np.std(hba_list)
    hbd_mean, hbd_std = np.mean(hbd_list), np.std(hbd_list)
    MW_mean, MW_std = np.mean(MW_list), np.std(MW_list)
    TPSA_mean, TPSA_std = np.mean(TPSA_list), np.std(TPSA_list)

    print('exp_dir: {}; hba mean: {}; hba std: {}; hbd mean: {}; hbd std: {}; '
          'MW mean: {}; MW std: {}; TPSA mean: {}; TPSA std: {}'.format(exp_dir,
                                                                               hba_mean, hba_std,
                                                                               hbd_mean, hbd_std,
                                                                               MW_mean, MW_std,
                                                                               TPSA_mean, TPSA_std))




