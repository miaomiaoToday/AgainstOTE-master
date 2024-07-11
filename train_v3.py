# import sys
# sys.path.append('.')
import os
import shutil
import argparse
from tqdm.auto import tqdm
import torch
from torch.nn.utils import clip_grad_norm_
import torch.utils.tensorboard
# import torch_geometric
# assert not torch_geometric.__version__.startswith('2'), 'Please use torch_geometric lower than version 2.0.0'
from torch_geometric.loader import DataLoader

from models.maskfill import MaskFillModelVN
from utils.datasets import *
from utils.transforms import *
from utils.misc import *
from utils.train import *
from models.pca_linear import PCALinear


def h_abnormaldim(lists):
    dim_sca_list, dim_vec_list = [], []
    for h_list in lists:
        # assuming num aug > 1
        abdim_sca = [each_h[0].size()[1] for each_h in h_list]
        abdim_vec = [each_h[1].size()[1] for each_h in h_list]
        min_dim_sca = min(abdim_sca)
        min_dim_vec = min(abdim_vec)
        dim_sca_list.append(min_dim_sca)
        dim_vec_list.append(min_dim_vec)
    min_dim_sca = min(dim_sca_list)
    min_dim_vec = min(dim_vec_list)
    for h_list in lists:
        for each_h in h_list:
            each_h[0] = each_h[0][:, :min_dim_sca]
            each_h[1] = each_h[1][:, :min_dim_vec]

    return lists


def nt_xent_loss(queries, keys, temperature=0.1):
    b = queries.shape[0]

    n = b * 2
    projs = torch.cat((queries, keys))
    logits = projs @ projs.t()

    mask = torch.eye(n).bool()
    logits = logits[~mask].reshape(n, n - 1)
    logits /= temperature

    labels = torch.cat(((torch.arange(b) + b - 1), torch.arange(b)), dim=0)
    loss = F.cross_entropy(logits, labels, reduction='sum')
    loss /= n
    return loss


import torch
from torch import nn
import torch.nn.functional as F


class ContrastiveLossELI5(nn.Module):
    def __init__(self, batch_size, temperature=0.5, verbose=True):
        super().__init__()
        self.batch_size = batch_size
        self.register_buffer("temperature", torch.tensor(temperature))
        self.verbose = verbose

    def forward(self, emb_i, emb_j):
        """
        emb_i and emb_j are batches of embeddings, where corresponding indices are pairs
        z_i, z_j as per SimCLR paper
        """
        z_i = F.normalize(emb_i, dim=1)
        z_j = F.normalize(emb_j, dim=1)

        representations = torch.cat([z_i, z_j], dim=0)
        similarity_matrix = F.cosine_similarity(representations.unsqueeze(1), representations.unsqueeze(0), dim=2)
        if self.verbose: print("Similarity matrix\n", similarity_matrix, "\n")

        def l_ij(i, j):
            z_i_, z_j_ = representations[i], representations[j]
            sim_i_j = similarity_matrix[i, j]
            if self.verbose: print(f"sim({i}, {j})={sim_i_j}")

            numerator = torch.exp(sim_i_j / self.temperature)
            one_for_not_i = torch.ones((2 * self.batch_size,)).to(emb_i.device).scatter_(0, torch.tensor([i]).to(
                emb_i.device), 0.0)
            if self.verbose: print(f"1{{k!={i}}}", one_for_not_i)

            denominator = torch.sum(
                one_for_not_i * torch.exp(similarity_matrix[i, :] / self.temperature)
            )
            if self.verbose: print("Denominator", denominator)

            loss_ij = -torch.log(numerator / denominator)
            if self.verbose: print(f"loss({i},{j})={loss_ij}\n")

            return loss_ij.squeeze(0)

        N = self.batch_size
        loss = 0.0
        for k in range(0, N):
            loss += l_ij(k, k + N) + l_ij(k + N, k)
        return 1.0 / (2 * N) * loss


def model_forward(model, batch, compose_noise):
    loss, loss_frontier, loss_pos, loss_cls, loss_edge, loss_real, loss_fake, loss_surf, h_compose, h_ligand, h_protein = model.get_loss(

        pos_real=batch.pos_real,
        y_real=batch.cls_real.long(),
        # p_real = batch.ind_real.float(),    # Binary indicators: float
        pos_fake=batch.pos_fake,

        edge_index_real=torch.stack([batch.real_compose_edge_index_0, batch.real_compose_edge_index_1], dim=0),
        edge_label=batch.real_compose_edge_type,

        index_real_cps_edge_for_atten=batch.index_real_cps_edge_for_atten,
        tri_edge_index=batch.tri_edge_index,
        tri_edge_feat=batch.tri_edge_feat,

        compose_feature=batch.compose_feature.float(),
        compose_pos=batch.compose_pos + compose_noise,
        idx_ligand=batch.idx_ligand_ctx_in_compose,
        idx_protein=batch.idx_protein_in_compose,

        y_frontier=batch.ligand_frontier,
        idx_focal=batch.idx_focal_in_compose,
        pos_generate=batch.pos_generate,
        idx_protein_all_mask=batch.idx_protein_all_mask,
        y_protein_frontier=batch.y_protein_frontier,

        compose_knn_edge_index=batch.compose_knn_edge_index,
        compose_knn_edge_feature=batch.compose_knn_edge_feature,
        real_compose_knn_edge_index=torch.stack(
            [batch.real_compose_knn_edge_index_0, batch.real_compose_knn_edge_index_1], dim=0),
        fake_compose_knn_edge_index=torch.stack(
            [batch.fake_compose_knn_edge_index_0, batch.fake_compose_knn_edge_index_1], dim=0),
    )

    return loss, loss_frontier, loss_pos, loss_cls, loss_edge, loss_real, loss_fake, loss_surf, h_compose, h_ligand, h_protein


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, default='./configs/train.yml')
    parser.add_argument('--device', type=str, default='cuda')
    parser.add_argument('--logdir', type=str, default='./logs')
    args = parser.parse_args()

    # Load configs
    config = load_config(args.config)
    config_name = os.path.basename(args.config)[:os.path.basename(args.config).rfind('.')]
    seed_all(config.train.seed)
    if config.train.use_apex:
        from apex import amp

    # Logging
    log_dir = get_new_log_dir(args.logdir, prefix=config_name)
    ckpt_dir = os.path.join(log_dir, 'checkpoints')
    ckpt_dir2 = os.path.join(log_dir, 'checkpoints_proj')
    os.makedirs(ckpt_dir, exist_ok=True)
    os.makedirs(ckpt_dir2, exist_ok=True)
    logger = get_logger('train', log_dir)
    writer = torch.utils.tensorboard.SummaryWriter(log_dir)
    logger.info(args)
    logger.info(config)
    shutil.copyfile(args.config, os.path.join(log_dir, os.path.basename(args.config)))
    shutil.copytree('./models', os.path.join(log_dir, 'models'))

    # Transforms
    protein_featurizer = FeaturizeProteinAtom()
    ligand_featurizer = FeaturizeLigandAtom()
    masking = get_mask(config.train.transform.mask)
    composer = AtomComposer(protein_featurizer.feature_dim, ligand_featurizer.feature_dim, config.model.encoder.knn)

    edge_sampler = EdgeSample(config.train.transform.edgesampler)
    cfg_ctr = config.train.transform.contrastive
    contrastive_sampler = ContrastiveSample(cfg_ctr.num_real, cfg_ctr.num_fake, cfg_ctr.pos_real_std,
                                            cfg_ctr.pos_fake_std, config.model.field.knn)
    transform = Compose([
        RefineData(),
        LigandCountNeighbors(),
        protein_featurizer,
        ligand_featurizer,
        masking,
        composer,

        FocalBuilder(),
        edge_sampler,
        contrastive_sampler,
    ])

    # Datasets and loaders
    logger.info('Loading dataset...')
    dataset, subsets = get_dataset(
        config=config.dataset,
        transform=transform,
    )
    train_set, val_set = subsets['train'], subsets['test']
    follow_batch = []
    collate_exclude_keys = ['ligand_nbh_list']
    """num_workers, not okay in windows, it should be 0"""
    train_iterator = inf_iterator(DataLoader(
        train_set,
        batch_size=config.train.batch_size,
        shuffle=True,
        # num_workers = config.train.num_workers,
        num_workers=0,
        pin_memory=config.train.pin_memory,
        follow_batch=follow_batch,
        exclude_keys=collate_exclude_keys,
    ))
    val_loader = DataLoader(val_set, config.train.batch_size, shuffle=False, follow_batch=follow_batch,
                            exclude_keys=collate_exclude_keys, )

    # Model
    logger.info('Building model...')
    if config.model.vn == 'vn':
        model = MaskFillModelVN(
            config.model,
            num_classes=contrastive_sampler.num_elements,
            num_bond_types=edge_sampler.num_bond_types,
            protein_atom_feature_dim=protein_featurizer.feature_dim,
            ligand_atom_feature_dim=ligand_featurizer.feature_dim,
        ).to(args.device)
    print('Num of parameters is', np.sum([p.numel() for p in model.parameters()]))

    # Optimizer and scheduler
    optimizer = get_optimizer(config.train.optimizer, model)
    scheduler = get_scheduler(config.train.scheduler, optimizer)
    if config.train.use_apex:
        model, optimizer = amp.initialize(model, optimizer, opt_level='O1')

    if config.train.finetune:
        ckpt = torch.load(config.model.checkpoint, map_location=args.device)
        # pretrained_dict = {k: v for k, v in ckpt.items() if k in model.state_dict()}
        model.load_state_dict(ckpt['model'])

    """loss from torch"""
    L1_Loss = torch.nn.L1Loss()
    MSE_Loss = torch.nn.MSELoss()

    # another model
    model_proj = PCALinear(indim=128, outdim=1, blocknum=4).to(args.device)
    # norm_obj = config.train.normalize_obj
    # pcalinear = nn.Linear(128, 1).to(args.device)

    def train(it):
        # model.train()  has been moved to the end of validation function
        optimizer.zero_grad()

        """T times of model forward for T times of structural displacement
        Return: 
            [loss_list]: the list of each loss for backward from model forward
            [h_list]: the list of each projected features of compose
            ...
        """
        pos_batch = next(train_iterator).to(args.device)
        # pos_loss_list, pos_h_compose_list, pos_h_ligand_list, pos_h_protein_list = [], [], [], []
        pos_loss_list, pos_h_ligand_list = [], []
        for it_aug in range(config.train.num_aug):
            compose_noise = torch.randn_like(pos_batch.compose_pos) * config.train.pos_noise_std
            loss, loss_frontier, loss_pos, loss_cls, loss_edge, loss_real, loss_fake, loss_surf, h_compose, h_ligand, h_protein = model_forward(
                model, pos_batch, compose_noise)

            # continue model
            h_ligand = model_proj(h_ligand)

            pos_loss_list.append(loss)
            # pos_h_compose_list.append(h_compose)
            pos_h_ligand_list.append(h_ligand)
            # pos_h_protein_list.append(h_protein)  # [(list len: aug num) [(list len: 2) [1, 256][1, 192]]]

        neg_batch = next(train_iterator).to(args.device)
        # neg_loss_list, neg_h_compose_list, neg_h_ligand_list, neg_h_protein_list = [], [], [], []
        neg_loss_list, neg_h_ligand_list = [], []
        for it_aug in range(config.train.num_aug):
            compose_noise = torch.randn_like(neg_batch.compose_pos) * config.train.pos_noise_std
            loss, loss_frontier, loss_pos, loss_cls, loss_edge, loss_real, loss_fake, loss_surf, h_compose, h_ligand, h_protein = model_forward(
                model, neg_batch, compose_noise)

            # continue model
            h_ligand = model_proj(h_ligand)

            neg_loss_list.append(loss)
            # neg_h_compose_list.append(h_compose)
            neg_h_ligand_list.append(h_ligand)
            # neg_h_protein_list.append(h_protein)

        # the predicton loss original
        loss1 = 0
        for loss in pos_loss_list:
            loss1 += loss
        for loss in neg_loss_list:
            loss1 += loss
        loss1 = loss1 / (2 * config.train.num_aug)

        """"""
        """2. projector to project h_composes to unified space."""
        # finished in maskfill for now
        """"""

        """2. New loss, may be we here should only thake the 1st h_compose and use MLP or other tech to project"""
        [pos_h_ligand_list, neg_h_ligand_list] = h_abnormaldim([pos_h_ligand_list, neg_h_ligand_list])


        # ligand, obtain scas and vecs for ...
        pos_ligand_scas, pos_ligand_vecs = [], []
        neg_ligand_scas, neg_ligand_vecs = [], []
        for pos_h_ligand in pos_h_ligand_list:
            pos_ligand_scas.append(pos_h_ligand[0])
            pos_ligand_vecs.append(pos_h_ligand[1])
        for neg_h_ligand in neg_h_ligand_list:
            neg_ligand_scas.append(neg_h_ligand[0])
            neg_ligand_vecs.append(neg_h_ligand[1])

        pos_ligand_scas = torch.cat(pos_ligand_scas, dim=0)     # [2, 256]
        pos_ligand_vecs = torch.cat(pos_ligand_vecs, dim=0)     # [2, 192]
        neg_ligand_scas = torch.cat(neg_ligand_scas, dim=0)     # [2, 256]
        neg_ligand_vecs = torch.cat(neg_ligand_vecs, dim=0)     # [2, 192]

        # construct batch for contrastive loss, for self and other differece, lack self loss
        # loss6 = nt_xent_loss(queries, keys, temperature=0.1)
        loss_eli5 = ContrastiveLossELI5(batch_size=config.train.num_aug, temperature=1.0, verbose=False)
        loss3 = loss_eli5(pos_ligand_scas, neg_ligand_scas) + loss_eli5(pos_ligand_vecs, neg_ligand_vecs)

        # centroid loss
        # loss5 = MSE_Loss(pos_ligand_scas, neg_ligand_scas) + MSE_Loss(pos_ligand_vecs, neg_ligand_vecs) * 1e2
        loss5 = MSE_Loss(pos_ligand_scas, neg_ligand_scas) + MSE_Loss(pos_ligand_vecs, neg_ligand_vecs)

        # self normalize loss
        # using cycle cl loss
        loss7 = loss_eli5(neg_ligand_scas, pos_ligand_scas) + loss_eli5(neg_ligand_vecs, pos_ligand_vecs)

        loss = loss1 + loss3 - loss5 + loss7
        # loss = loss + loss2 - loss3
        """"""
        if config.train.use_apex:
            with amp.scale_loss(loss, optimizer) as scaled_loss:
                scaled_loss.backward()
        else:
            loss.backward()
        orig_grad_norm = clip_grad_norm_(model.parameters(), config.train.max_grad_norm,
                                         error_if_nonfinite=True)  # 5% running time, from true to false to disable  Runtime Error The total norm of order 2.0 for gradients from `parameters` is non-finite, so it cannot be clipped. To disable this error and scale the gradients by the non-finite norm anyway, set `error_if_nonfinite=False`
        optimizer.step()

        logger.info(
            '[Train] Iter %d | Loss %.6f | Loss1 %.6f | Loss2 %.6f | Loss3 %.6f | Loss4 %.6f | Loss5 %.6f | Loss6 %.6f | Loss7 %.6f' % (
                it, loss.item(), loss1.item(), loss3.item(), loss3.item(), loss5.item(), loss5.item(),
                loss7.item(), loss7.item()
            ))
        writer.add_scalar('train/loss', loss, it)
        writer.add_scalar('train/loss_fron', loss_frontier, it)
        writer.add_scalar('train/loss_pos', loss_pos, it)
        writer.add_scalar('train/loss_cls', loss_cls, it)
        writer.add_scalar('train/loss_edge', loss_edge, it)
        writer.add_scalar('train/loss_real', loss_real, it)
        writer.add_scalar('train/loss_fake', loss_fake, it)
        writer.add_scalar('train/loss_surf', loss_surf, it)
        writer.add_scalar('train/lr', optimizer.param_groups[0]['lr'], it)
        writer.add_scalar('train/grad', orig_grad_norm, it)
        writer.flush()


    def validate(it):
        sum_loss, sum_n = np.zeros(5 + 2 + 1), 0  # num of loss
        with torch.no_grad():
            model.eval()
            for batch in tqdm(val_loader, desc='Validate'):
                batch = batch.to(args.device)
                loss, loss_frontier, loss_pos, loss_cls, loss_edge, loss_real, loss_fake, loss_surf, h_compose, h_ligand, h_protein = model.get_loss(
                    pos_real=batch.pos_real,
                    y_real=batch.cls_real.long(),
                    pos_fake=batch.pos_fake,

                    edge_index_real=torch.stack([batch.real_compose_edge_index_0, batch.real_compose_edge_index_1],
                                                dim=0),
                    edge_label=batch.real_compose_edge_type,

                    index_real_cps_edge_for_atten=batch.index_real_cps_edge_for_atten,
                    tri_edge_index=batch.tri_edge_index,
                    tri_edge_feat=batch.tri_edge_feat,

                    compose_feature=batch.compose_feature.float(),
                    compose_pos=batch.compose_pos,
                    idx_ligand=batch.idx_ligand_ctx_in_compose,
                    idx_protein=batch.idx_protein_in_compose,

                    y_frontier=batch.ligand_frontier,
                    idx_focal=batch.idx_focal_in_compose,
                    pos_generate=batch.pos_generate,
                    idx_protein_all_mask=batch.idx_protein_all_mask,
                    y_protein_frontier=batch.y_protein_frontier,

                    compose_knn_edge_index=batch.compose_knn_edge_index,
                    compose_knn_edge_feature=batch.compose_knn_edge_feature,
                    real_compose_knn_edge_index=torch.stack(
                        [batch.real_compose_knn_edge_index_0, batch.real_compose_knn_edge_index_1], dim=0),
                    fake_compose_knn_edge_index=torch.stack(
                        [batch.fake_compose_knn_edge_index_0, batch.fake_compose_knn_edge_index_1], dim=0),
                )
                loss_list = [loss, loss_frontier, loss_pos, loss_cls, loss_edge, loss_real, loss_fake, loss_surf]
                sum_loss = sum_loss + np.array([torch.nan_to_num(l).item() for l in loss_list])
                sum_n += 1
        avg_loss = sum_loss / sum_n

        if config.train.scheduler.type == 'plateau':
            scheduler.step(avg_loss[0])
        elif config.train.scheduler.type == 'warmup_plateau':
            scheduler.step_ReduceLROnPlateau(avg_loss[0])
        else:
            scheduler.step()

        logger.info(
            '[Validate]  Iter %d | Loss %.6f | Loss(Fron) %.6f | Loss(Pos) %.6f | Loss(Cls) %.6f | Loss(Edge) %.6f | Loss(Real) %.6f | Loss(Fake) %.6f  | Loss(Surf) %.6f' % (
                it, *avg_loss,
            ))
        writer.add_scalar('val/loss', avg_loss[0], it)
        writer.add_scalar('val/loss_fron', avg_loss[1], it)
        writer.add_scalar('val/loss_pos', avg_loss[2], it)
        writer.add_scalar('val/loss_cls', avg_loss[3], it)
        writer.add_scalar('val/loss_edge', avg_loss[4], it)
        writer.add_scalar('val/loss_real', avg_loss[5], it)
        writer.add_scalar('val/loss_fake', avg_loss[6], it)
        writer.add_scalar('val/loss_surf', avg_loss[7], it)
        writer.flush()
        return avg_loss


    try:
        model.train()
        model_proj.train()
        for it in range(1, config.train.max_iters + 1):
            try:
                train(it)
            except RuntimeError as e:
                logger.error('Runtime Error ' + str(e))
            if it % config.train.val_freq == 0 or it == config.train.max_iters:
                validate(it)
                ckpt_path = os.path.join(ckpt_dir, '%d.pt' % it)
                torch.save({
                    'config': config,
                    'model': model.state_dict(),
                    'optimizer': optimizer.state_dict(),
                    'scheduler': scheduler.state_dict(),
                    'iteration': it,
                }, ckpt_path)
                ckpt_path2 = os.path.join(ckpt_dir2, '%d.pt' % it)
                torch.save({
                    'config': config,
                    'model': model_proj.state_dict(),
                    'optimizer': optimizer.state_dict(),
                    'scheduler': scheduler.state_dict(),
                    'iteration': it,
                }, ckpt_path2)
                model.train()
                # model_proj.train()
    except KeyboardInterrupt:
        logger.info('Terminating...')
