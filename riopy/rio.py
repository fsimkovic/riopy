
import numpy as np

from iotbx import pdb
from riopy.misc import count_res_in_frag 
from riopy.neighbors import Neighbors


class Rio(object):

    def __init__(self, h_model, h_target, sg, uc):
        self.h_model = h_model
        self.h_target = h_target
        self.space_group = sg
        self.unit_cell = uc

    def compute(self, max_dist=1.5, cells=1):
        max_dist_sq = max_dist**2

        rio_score = 0
        
        h_model_atoms = self.h_model.atoms()
        h_model_atoms_xyz = np.array(h_model_atoms.extract_xyz())

        for rotop, transop in Neighbors.ops(self.space_group, cells):
            h_target_cur = self.h_target.deep_copy()
            h_target_cur_atoms = h_target_cur.atoms()

            rt_a = self.unit_cell.fractionalize(sites_cart=h_target_cur_atoms.extract_xyz()) * rotop + transop
            h_target_cur_atoms_xyz = np.array(self.unit_cell.orthogonalize(sites_frac=rt_a))

            fragmented = set()
            for h_target_cur_atom_xyz in h_target_cur_atoms_xyz:
                diff = h_model_atoms_xyz - h_target_cur_atom_xyz
                dist_sq = np.sum(np.power(diff, 2), axis=1)
                for k in np.flatnonzero(dist_sq <= max_dist_sq):
                    fragmented.add(Rio.get_resseq(h_model_atoms[k]))

            rio_score += count_res_in_frag(sorted(list(fragmented)))

        return rio_score

    @staticmethod
    def get_resseq(a):
        return a.parent().parent().resseq_as_int()

