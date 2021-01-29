from pyDFTutils.ase_utils.frozenphonon import calculate_phonon
from pyDFTutils.siesta.mysiesta import get_species, MySiesta
from ase.io import read, write
import numpy as np
from phonopy import load, Phonopy
import matplotlib.pyplot as plt
import os
import copy

def do_siesta_calculation(xc='PBEsol', use_spin=True, netcharge=0,relax=True, phonon=True):
    atoms = read("pure_rLAO_relaxed_sqrt2_sqrt2_sqrt2.vasp")
    # add spin or not. number from 0 here because it is python. 
    m=np.zeros(len(atoms))
    m[21]=2
    m[7]=1
    if use_spin:
        atoms.set_initial_magnetic_moments(m)
    pseudo_path, species = get_species(atoms, xc=xc, rel='sr')
    fdf_arguments={'Meshcutoff':'800 Ry', 'DM.Tolerance':'0.0001','ElectronicTemperature':'100 K'}
    if netcharge !=0:
        fdf_arguments.update({'NetCharge':f'{netcharge}', 'SimulateDoping':'.true.'})
    calc = MySiesta(
                    label='siesta',
                    kpts=[8,8,8],
                    xc=xc,
                    basis_set='DZP',
                    species=species,
                    pseudo_path=pseudo_path,
		    fdf_arguments=fdf_arguments)
    if not os.path.exists('Results'):
        os.makedirs('Results')
    # relax
    if relax:
        calc.label='relax/siesta'
        atoms=calc.relax(atoms, MaxForceTol=1e-3, MaxStressTol=0.1, NumCGSteps=100)
        write('Results/relaxed.vasp', atoms, vasp5=True) 
    #dos&potential
    dos_calc=copy.deepcopy(calc)
    dos_calc.label='dos/sietsa'
    fdf=dos_calc['fdf_arguments']
    fdf.update(
{'WriteEigenvalues': '.true.', 
'SaveElectrostaticPotential': ' .true.',
'SaveRho': '.true.',
'SaveTotalCharge': '.true.',
'SaveIonicCharge': '.true.',
'SaveDeltaRho': '.true.',
'SaveTotalPotential': ' .true.',
'ProjectedDensityOfStates': ['-70.00 5.0 0.015 3000 eV'],
'BandLinesScale': 'pi/a',
'BandLines':['1  0.0 0.0 0.0 \Gamma',
            '22 1.0 0.0 0.0 X',
            '22 1.0 1.0 0.0 M',
            '33 0.0 0.0 0.0 \Gamma',
            '39 1.0 1.0 1.0 R',
            '33 1.0 0.0 0.0 X']
})
    dos_calc.set_fdf_arguments(fdf)
    dos_calc.calculate(atoms)
    os.system('cp dos/siesta.PDOS Results/siesta.PDOS') 
    os.system('cp dos/siesta.VH Results/siesta.VH') 

    # phonon
    if phonon:
        phonon_calc=copy.deepcopy(calc)
        #calculate_phonon(atoms, calc=phonon_calc, ndim=np.diag([2.,  2.,  2.]),parallel=False,symprec=1e-3)
        #calculate_phonon(atoms, calc=phonon_calc, ndim=np.array([[-1.,  1.,  1.],
        #   [ 1., -1.,  1.],
        #   [ 1.,  1., -1.]]), parallel=False,symprec=1e-3)


if __name__=='__main__':
    do_siesta_calculation()

