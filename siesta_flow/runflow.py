from pyDFTutils.ase_utils.frozenphonon import calculate_phonon
from pyDFTutils.siesta.mysiesta import get_species, MySiesta
from ase.io import read, write
import numpy as np
from phonopy import load, Phonopy
import matplotlib.pyplot as plt
import os
import copy
from siesta_flow.pdos import gen_pdos_figure, plot_layer_pdos
from siesta_flow.flow import do_relax_calculation, do_scf_calculation
from ase.units import Ry

def run():
    # generate atom structure
    atoms=read("siesta.vasp")

    ##### 1) spin or no_spin??  number from 0 here because it is python. 
    # number 8 in vasp, here is 7
    m=np.zeros(len(atoms))
    m[7]=2
    m[15]=1
    atoms.set_initial_magnetic_moments(m)

    ##### parameter xc
    xc='PBEsol'
    pseudo_path, species = get_species(atoms, xc=xc, rel='sr')
    fdf_arguments={ 'DM.Tolerance':'0.0001','ElectronicTemperature':'100 K','DM.NumberPulay': '6', 'DM.MixingWeight': '0.1', 'SCF.Mixer.Method':'Pulay','MaxSCFIterations':200}
    ###### 1)neutral or charged
    netcharge=0
    if netcharge !=0:
        fdf_arguments.update({'NetCharge':f'{netcharge}', 'SimulateDoping':'.true.'})
    calc = MySiesta(
                    label='siesta',
                    kpts=[7,7,7],
                    xc=xc,
                    basis_set='DZP',
                    mesh_cutoff=400*Ry,
                    # energy_shift convergence: smaller, better,larger calculation 
                    energy_shift=0.1,
                    species=species,
                    ### 1) spin='non-polarized', 'collinear', 'non-colliner','spin-orbit'
                    spin='collinear',
                    pseudo_path=pseudo_path,
            fdf_arguments=fdf_arguments)


    if not os.path.exists('Results'):
        os.makedirs('Results')

    ####### 1) relax or not; 2)variable cell or not;  3) maxforce: 0.01 or 0.001)
    do_relax_calculation(atoms, calc, VariableCell=False, MaxForceTol=1e-2, MaxStressTol=0.1, NumCGSteps=1000,TypeOfRun='cg')

    ###### 1) UseSaveDM or not;
    do_scf_calculation(atoms, calc, UseDM=True, dos=True, band_structure=True, potential=False)
    
if __name__=='__main__':
    run()
