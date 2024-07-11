import torch
from torch.nn import Module
from torch.nn import functional as F

from .encoders import get_encoder_vn
from .fields import get_field_vn
from .common import *
from .embedding import AtomEmbedding
from .frontier import FrontierLayerVN
from .position import PositionPredictor
# from .debug import check_true_bonds_len, check_pred_bonds_len
from utils.misc import unique


def MolPCA(X, k=128, center=True):
    """PCA for mol feature"""
    n = X.size()[0]
    ones = torch.ones(n).view([n, 1])
    h = ((1 / n) * torch.mm(ones, ones.t())) if center else torch.zeros(n * n).view([n, n])
    H = torch.eye(n) - h
    H = H.cuda()
    X_center = torch.mm(H.double(), X.double())
    u, s, v = torch.svd(X_center)
    # components = v[:k].t()
    components = v[:k]
    # explained_variance = torch.mul(s[:k], s[:k])/(n-1)
    return components


def h_norm(h):
    """norm the feature for learning"""
    sca = h[0]
    vec = h[1]
    try:
        a, b, c = vec.size()
        vec = torch.reshape(vec, [a, b * c])
    except:
        pass
    sca = MolPCA(sca)
    vec = MolPCA(vec)
    return [sca, vec]


class PCALinear(Module):
    def __init__(self, indim=128, outdim=1, blocknum=4):
        super().__init__()
        self.linearout = nn.Linear(indim, outdim)
        self.linearin = nn.Linear(indim, indim)
        self.layernorm = nn.LayerNorm(indim, eps=1e-05, elementwise_affine=True)
        self.relu = nn.ReLU(inplace=False)
        self.blocknum = blocknum

    def forward(self, h_ligand):

        # h_compose = h_norm(h_compose)
        h_ligand = h_norm(h_ligand)         #  [128, 256] [128, 192]   # after transpose, [256, 128] [192, 128]
        # h_protein = h_norm(h_protein)

        for i in range(self.blocknum):
            h_ligand[0] = self.layernorm(h_ligand[0].transpose(-2, -1).float()).transpose(-2, -1)
            h_ligand[1] = self.layernorm(h_ligand[1].transpose(-2, -1).float()).transpose(-2, -1)

            h_ligand[0] = self.linearin(h_ligand[0].transpose(-2, -1).float()).transpose(-2, -1)
            h_ligand[1] = self.linearin(h_ligand[1].transpose(-2, -1).float()).transpose(-2, -1)

            h_ligand[0] = self.relu(h_ligand[0].transpose(-2, -1).float()).transpose(-2, -1)
            h_ligand[1] = self.relu(h_ligand[1].transpose(-2, -1).float()).transpose(-2, -1)

        h_ligand[0] = self.linearout(h_ligand[0].transpose(-2, -1).float()).transpose(-2, -1)
        h_ligand[1] = self.linearout(h_ligand[1].transpose(-2, -1).float()).transpose(-2, -1)

        h_ligand[0] = self.relu(h_ligand[0].transpose(-2, -1).float()).transpose(-2, -1)
        h_ligand[1] = self.relu(h_ligand[1].transpose(-2, -1).float()).transpose(-2, -1)

        return h_ligand







