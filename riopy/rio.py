#  BSD 3-Clause License
#
#  Copyright (c) 2018, University of Liverpool
#  All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
__author__ = "Felix Simkovic"
__date__ = "16 Jan 2018"
__version__ = "1.0"

import numpy as np
import os

from iotbx import pdb, reflection_file_reader
from scitbx.array_family import flex

from riopy.neighbors import Neighbors
from riopy.misc import calculate_origin_shift, count_res_in_frag, move_model_to_target


class RioScore(object):
    """RIO score data container"""
    total = 0
    norm_model = 0.0
    norm_target = 0.0


class Rio(object):
    """RIO score calculator"""

    def __init__(self, mtz, target, model):
        self._mtz = None
        self._target = None
        self._model = None

        self._unit_cell = None
        self._space_group = None
        self._hierarchy_target = None
        self._hierarchy_model = None

        self.mtz = mtz
        self.target = target
        self.model = model

    def compute(self, atom="CA", max_dist=1.5, ncells=1):
        """Compute the RIO score

        Parameters
        ----------
        atom : str
           The atom type for which we calculate the RIO score [default: CA]
        max_dist : float
           The distance for defining a contact (``d < max_dist``) [default: 1.5]
        ncels : int
           The number of neighboring unit cells to search in [default: 1]

        Returns
        -------
        int
           The RIO score

        """
        org_model = self._model
        self.model = self.correct_model_origin()

        selection = "name " + atom.lower()
        h_mod_atm = Rio.select_from_hierarchy(self._hierarchy_model, selection)
        h_tar_atm = Rio.select_from_hierarchy(self._hierarchy_target, selection)

        rio = RioScore()

        h_mod_xyz = np.array(h_mod_atm.atoms().extract_xyz())
        for rotop, transop in Neighbors.ops(self._space_group, ncells):
            # TODO: orthogonalize rotop and transop instead of all xyz
            rt_a = self._unit_cell.fractionalize(sites_cart=h_tar_atm.atoms().extract_xyz()) * rotop + transop
            h_tar_xyz = np.array(self._unit_cell.orthogonalize(sites_frac=rt_a))

            con = []
            for o in h_tar_xyz:
                dists = self.calculate_dists_2d(h_mod_xyz, o)
                con.extend(self.find_contacts(dists, max_dist))
            con_as_resseqs = list(Rio.convert_contacts_to_resseqs(h_mod_atm, con))

            rio.total += count_res_in_frag(con_as_resseqs)

        os.unlink(self._model)
        self.model = org_model

        rio.norm_model = round(rio.total / float(h_mod_xyz.shape[0]), 3)
        rio.norm_target = round(rio.total / float(h_tar_xyz.shape[0]), 3)
        return rio

    def correct_model_origin(self):
        """Correct the origin of the model according to the Xtal data and the target"""
        origin_shift = calculate_origin_shift(self.mtz, self.target, self.model)
        return move_model_to_target(self.target, self.model, origin_shift)

    @property
    def mtz(self):
        return self._mtz

    @mtz.setter
    def mtz(self, mtz):
        content = reflection_file_reader.any_reflection_file(file_name=mtz).file_content()
        self._space_group = content.space_group_name()
        self._unit_cell = content.crystals()[0].unit_cell()
        self._mtz = mtz

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        self._hierarchy_model = pdb.pdb_input(model).construct_hierarchy()
        self._model = model

    @property
    def target(self):
        return self._target

    @target.setter
    def target(self, target):
        self._hierarchy_target = pdb.pdb_input(target).construct_hierarchy()
        self._target = target

    @staticmethod
    def calculate_dists_2d(a, b, decimals=3):
        # TODO: Research alternatives by using np.dot
        return np.sqrt(np.sum(np.power(a - b, 2), axis=1)).round(decimals=decimals)

    @staticmethod
    def find_contacts(a, dist):
        return np.flatnonzero(a <= dist)

    @staticmethod
    def select_from_hierarchy(h, sel):
        return h.select(h.atom_selection_cache().selection(sel))

    @staticmethod
    def convert_contacts_to_resseqs(h, con):
        for c in list(sorted(set(con))):
            yield h.atoms()[c].parent().parent().resseq_as_int()
