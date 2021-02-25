from pyDFTutils.ase_utils.frozenphonon import calculate_phonon
from pyDFTutils.siesta.mysiesta import get_species, MySiesta
from ase.io import read, write
import numpy as np
from phonopy import load, Phonopy
import matplotlib.pyplot as plt
import os
import copy
from siesta_flow.pdos import gen_pdos_figure, plot_layer_pdos 

def do_relax_calculation(atoms, calc, MaxForceTol=1e-2, MaxStressTol=0.1, NumCGSteps=1000,VariableCell=False, TypeOfRun='cg'):
    calc.label='relax/siesta'
    atoms=calc.relax(atoms, MaxForceTol=MaxForceTol, MaxStressTol=MaxStressTol, NumCGSteps=NumCGSteps, VariableCell=VariableCell,TypeOfRun=TypeOfRun)
    write('Results/relaxed.vasp', atoms, vasp5=True) 

    os.system('cp relax/siesta.out Results/siesta_relax.out') 
    os.system('cp relax/siesta.fdf Results/siesta_relax.fdf') 
    os.system('cp relax/siesta.XV Results/siesta.XV') 

def do_scf_calculation(atoms, calc, dos=True, band_structure=True, potential=False, UseDM=True):
    if not os.path.exists('dos'):
        os.makedirs('dos')
    pwd=os.getcwd()
    if os.path.exists("dos/siesta.DM"):
        pass
    elif os.path.exists("relax/siesta.DM"):
        os.system(f"cp {pwd}/relax/siesta.DM dos/siesta.DM")
    dos_calc=copy.deepcopy(calc)
    dos_calc.label='dos/siesta'
    fdf=dos_calc['fdf_arguments']
    ###### DM: 会同步relax/DM file,所以如果算错了，这里要写成false.
    fdf.update({'DM.UseSaveDM':UseDM})
    if dos:
        fdf.update({'WriteEigenvalues': '.true.', 
    		'ProjectedDensityOfStates': ['-70.00 30.0 0.015 3000 eV'],
            'PDOS.kgrid_Monkhorst_Pack': ['7 0 0 0.0',
                                          '0 7 0 0.0',
                                          '0 0 7 0.0']})

    if band_structure:
        fdf.update({'BandLinesScale': 'pi/a',
                   'BandLines':['1  0.0 0.0 0.0 \Gamma',
                               '22 3.0 0.0 0.0 X',
                               '22 3.0 3.0 0.0 M',
                               '33 0.0 0.0 0.0 \Gamma',
                               '39 3.0 3.0 3.0 R',
                               '33 3.0 0.0 0.0 X']
                   })
    if potential:
        fdf.update({'SaveElectrostaticPotential': ' .true.',
                    'SaveRho': '.true.',
                    'SaveTotalCharge': '.true.',
                    'SaveIonicCharge': '.true.',
                    'SaveDeltaRho': '.true.',
                    'SaveTotalPotential': ' .true.',
                    })

    dos_calc.set_fdf_arguments(fdf)
    print(fdf)
    dos_calc.calculate(atoms)
    os.system('cp dos/siesta.out Results/siesta_scf.out') 
    os.system('cp dos/siesta.fdf Results/siesta_scf.fdf') 


    os.system('cp dos/siesta.PDOS Results/siesta.PDOS') 
    os.system('cp dos/siesta.PDOS siesta.PDOS') 
    os.system('cp dos/siesta.DOS Results/siesta.DOS') 
    os.system('cp dos/siesta.VH Results/siesta.VH') 
   
    symbols=atoms.get_chemical_symbols()
    sdict={}
    for s in symbols:
        if s not in sdict:
            sdict[s]=f"{s}.{len(sdict)+1}"
    for s in set(symbols):
        gen_pdos_figure('siesta.PDOS', sdict[s], 0,-1,9 ,output_path='./Results', xlim=(-5,5), ylim=(-20,20))
    os.system('mv pdos*.dat Results/')

    for iatom in range(len(atoms)):
        gen_pdos_figure('siesta.PDOS', iatom+1, 0,-1,9 ,output_path='./Results', xlim=(-5,5), ylim=(-8,8))
    os.system('mv pdos*.dat Results/')
do_phonon_calculation=calculate_phonon

#def do_phonon_calculation(atoms, calc):
#        phonon_calc=copy.deepcopy(calc)
#        calculate_phonon(atoms, calc=phonon_calc, ndim=np.diag([2.,  2.,  2.]),parallel=False,symprec=1e-3)
#        #calculate_phonon(atoms, calc=phonon_calc, ndim=np.array([[-1.,  1.,  1.],
#        #   [ 1., -1.,  1.],
#        #   [ 1.,  1., -1.]]), parallel=False,symprec=1e-3)

