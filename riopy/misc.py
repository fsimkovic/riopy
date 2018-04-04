import logging
import os
import uuid

logger = logging.getLogger(__name__)


def calculate_origin_shift(mtz, target, model):
    """Calculate the origin shift between a model and a target structure"""
    from ample.util.shelxe import MRinfo
    logger.debug("Calculating origin shift for %s given %s and %s", model, target, mtz)
    mrinfo = MRinfo('shelxe', target, mtz)
    mrinfo.analyse(model)
    logger.debug("Origin shift calculated to be: %s", mrinfo.originShift)
    return mrinfo.originShift


def move_model_to_target(target, model, origin_shift):
    """Move a model to the target using ``csymmatch```"""
    from ample.util.csymmatch import Csymmatch
    logger.debug("Translating %s to %s with origin shift %s", target, model, origin_shift)
    csymmatch_pdb = model.replace(".pdb", "_csymmatch_{}.pdb".format(str(uuid.uuid1())))
    return Csymmatch().wrapModelToNative(model, target, origin=origin_shift, csymmatchPdb=csymmatch_pdb)


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
