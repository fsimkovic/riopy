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

import itertools

from halocell import HaloCell
from symops import SymmetryOperator
from scitbx.array_family import flex


class Neighbors(object):

    @staticmethod
    def ops(sg, n_uc):
        """Yield all rotation / translation operations in cartesian space
        
        Parameters
        ----------
        sg : str
           The space group
        n_uc : int
           The number of neighboring cells

        Yields
        ------
        list
           A generator of (rotation, translation) operators

        """
        for symop, cellop in itertools.product(SymmetryOperator.ops(sg), HaloCell.ops(n_uc)):
            rotop = symop.r().as_double()
            transop = tuple(flex.double(symop.t().as_double()) + flex.double(cellop))
            yield (rotop, transop)

