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

import logging
import os

logger = logging.getLogger(__name__)


def calculate_origin_shift(mtz, target, model):
    """Calculate the origin shift between a model and a target structure"""
    from ample.util.shelxe import MRinfo
    logger.debug("Calculating origin shift for %s given %s and %s", model, target, mtz)
    mrinfo = MRinfo('shelxe', target, mtz)
    mrinfo.analyse(model)
    try:
        map(os.remove, ["shelxe-input.ent", "shelxe-input.hkl"])
    except:
        pass
    logger.debug("Origin shift calculated to be: %s", mrinfo.originShift)
    return mrinfo.originShift


def move_model_to_target(target, model, origin_shift):
    """Move a model to the target using ``csymmatch```"""
    from ample.util.csymmatch import Csymmatch
    logger.debug("Translating %s to %s with origin shift %s", target, model, origin_shift)
    return Csymmatch().wrapModelToNative(model, target, origin=origin_shift)


def count_res_in_frag(residues, min_frag_len=3):
    """Count the number of fragments in a list
    
    Parameters
    ----------
    residues : list
       A list of residue indexes
    min_frag_len : int
       The minimum length in residues of a fragment [default: 3]

    Returns
    -------
    int
       The number of fragments

    """
    if len(residues) < 1:
        return 0

    def valid_frag_len(c):
        return c >= min_frag_len

    res_count = 0
    frag_len = 1

    previous = residues[0]
    for current in residues[1:]:
        if current - 1 == previous:
            frag_len += 1
        else:
            if valid_frag_len(frag_len):
                res_count += frag_len
            frag_len = 1
        previous = current

    if valid_frag_len(frag_len):
        res_count += frag_len

    return res_count
