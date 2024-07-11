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

    """get multiple exp_dirs, ! need checking"""
    mols_root = r'./outputs'
    mols_dirs = os.listdir(mols_root)  # confused sequences of dirs. ignored.
    exp_dirs = [os.path.join(mols_root, mols_dir) for mols_dir in mols_dirs]

    qed_mean_list, sa_mean_list, lipinski_mean_list, logp_mean_list = [], [], [], []
    """do each exp_dir."""
    for exp_dir in exp_dirs:
        save_path = os.path.join(exp_dir, 'samples_all.pt')  # get_saved_filename(exp_dir))

        logger = get_logger('evaluate', exp_dir)
        logger.info(args)
        logger.info(save_path)
        print(exp_dir)
        samples = torch.load(save_path, map_location='cpu')

        sim_with_train = SimilarityWithTrain()
        results = []

        for i, data in enumerate(tqdm(samples.finished, desc='All')):
            try:
                mol = reconstruct_from_generated_with_edges(data)
                _, _, _, hba, hbd = get_chem(mol)
                MW = rdmd._CalcMolWt(mol)
                TPSA = MolSurf.TPSA(mol)
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

        qed_list, sa_list, lipinski_list, logp_list = [], [], [], []
        for result in results:
            qed_list.append(result['qed'])
            sa_list.append(result['sa'])
            lipinski_list.append(result['lipinski'])
            logp_list.append(result['logp'])
        qed_mean, qed_std = np.mean(qed_list), np.std(qed_list)
        sa_mean, sa_std = np.mean(sa_list), np.std(sa_list)
        lipinski_mean, lipinski_std = np.mean(lipinski_list), np.std(lipinski_list)
        logp_mean, logp_std = np.mean(logp_list), np.std(logp_list)

        print('exp_dir: {}; qed mean: {}; qed std: {}; sa mean: {}; sa std: {}; '
              'lipinski mean: {}; lip std: {}; logp mean: {}; logp std: {}'.format(exp_dir,
                                                                                   qed_mean, qed_std,
                                                                                   sa_mean, sa_std,
                                                                                   lipinski_mean, lipinski_std,
                                                                                   logp_mean, logp_std))
        result_path = os.path.join(exp_dir, 'results.pt')
        torch.save(results, result_path)

        qed_mean_list += qed_list
        sa_mean_list += sa_list
        lipinski_mean_list += lipinski_list
        logp_mean_list += logp_list


    qed_mean_mean, qed_mean_std = np.mean(qed_mean_list), np.std(qed_mean_list)
    sa_mean_mean, sa_mean_std = np.mean(sa_mean_list), np.std(sa_mean_list)
    lipinski_mean_mean, lipinski_mean_std = np.mean(lipinski_mean_list), np.std(lipinski_mean_list)
    logp_mean_mean, logp_mean_std = np.mean(logp_mean_list), np.std(logp_mean_list)
    print('FINAL!!! mols_root: {}; qed mean: {}+-{}; sa mean: {}+-{}; lipinski mean: {}+-{}; logp mean: {}+-{}'.format(
        mols_root,
        qed_mean_mean, qed_mean_std,
        sa_mean_mean, sa_mean_std,
        lipinski_mean_mean, lipinski_mean_std,
        logp_mean_mean, logp_mean_std))
