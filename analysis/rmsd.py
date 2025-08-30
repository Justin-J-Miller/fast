# Author: Maxwell I. Zimmerman <mizimmer@wustl.edu>
# Contributors:
# Copywright (C) 2017, Washington University in St. Louis
# All rights reserved.
# Unauthorized copying of this file, via any medium, is strictly prohibited
# Proprietary and confidential


#######################################################################
# imports
#######################################################################


import gc
import glob
import mdtraj as md
import numpy as np
import os
from .base_analysis import base_analysis
from .. import tools


#######################################################################
# code
#######################################################################


class RMSDWrap(base_analysis):
    """Analysis wrapper for calculating state rmsds

    Parameters
    ----------
    base_struct : str or md.Trajectory,
        The base structure to compute RMSD against. This
        topology must match the structures to analyse. Can be provided
        as a pdb location or an md.Trajectory object.
    atom_indices : str or array,
        The atom indices to use for computing RMSD. Can be
        provided as a data file to load or an array.
    ref_struct : str or md.Trajectory,
        The reference structure to compute RMSD against. If topology
        does not match the structures to analyze, ref_atom_indices must
        be proved that match. Can be provided as a path or md.Trajectory.
    ref_atom_indices : str or array,
        The atom indicies to use for computing RMSD. Can be
        provided as a data file to load or an array.

    Attributes
    ----------
    output_name : str,
        The file containing rankings.
    """
    def __init__(
            self, base_struct, atom_indices=None, 
            ref_struct=None, ref_atom_indices=None):
        # determine base_struct
        self.base_struct = base_struct
        if type(base_struct) is md.Trajectory:
            self.base_struct_md = self.base_struct
        else:
            self.base_struct_md = md.load(base_struct)
        # determine atom indices
        self.atom_indices = atom_indices
        if type(atom_indices) is str:
            self.atom_indices_vals = np.genfromtxt(atom_indices, dtype=int)
        else:
            self.atom_indices_vals = self.atom_indices

        #Repeat for ref structure
        self.ref_struct = ref_struct
        if type(ref_struct) is md.Trajectory:
            self.ref_struct_md = self.ref_struct

        else:
            self.ref_struct_md = md.load(ref_struct)

        self.ref_atom_indices = ref_atom_indices
        if type(ref_atom_indices) is str:
            self.ref_atom_indices_vals = np.genfromtxt(ref_atom_indices, dtype=int)
        else:
            self.ref_atom_indices_vals = self.ref_atom_indices_vals            


    @property
    def class_name(self):
        return "RMSDWrap"

    @property
    def config(self):
        return {
            'base_struct': self.base_struct,
            'atom_indices': self.atom_indices,
            'ref_struct' : self.ref_struct,
            'ref_atom_indices' : self.ref_atom_indices
        }

    @property
    def analysis_folder(self):
        return None

    @property
    def base_output_name(self):
        return "rmsd_per_state"

    def run(self):
        # determine if file already exists
        if os.path.exists(self.output_name):
            pass
        else:
            # load centers
            centers = md.load(
                "./data/full_centers.xtc", top=self.base_struct_md)
            # calculate and save rmsds
            if self.ref_struct is None:
                rmsds = md.rmsd(centers, struct_sub, atom_indices=self.atom_indices_vals)
            else:
                rmsds = md.rmsd(centers, reference = self.ref_struct_md,
                        atom_indices = self.atom_indices_vals,
                        ref_atom_indices = self.ref_atom_indices_vals)
            np.save(self.output_name, rmsds)
        
