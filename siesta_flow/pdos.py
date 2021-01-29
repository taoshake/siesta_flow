import numpy as np
import matplotlib.pyplot as plt
import os
plt.rcParams.update({'font.size': 15})
ldict={-1:'all',0:'s', 1:'p', 2:'d', 3:'f'}
lmdict={0:{0:'s'},
        1:{-1:'py',
            0:'pz',
            1:'px',
            9:'p-all'},
        2:{-2:'dxy', -1:'dyz', 0:'dz2', 1:'dxz', 2:'dx2-y2',9:'d-all'},
        -1:{9:'all'}
        }

def read_efermi(pdos_fname):
    with open(pdos_fname) as myfile:
        lines=myfile.readlines()
        efermi=float(lines[3].strip()[:-15].split()[2])
    return efermi

def get_pdos_data(pdos_fname, iatom, n, l, m):
    outfile=f"pdos_{iatom}_{n}{lmdict[l][m]}.dat"
    inp=f"""{pdos_fname}
{outfile}
{iatom}
{n}
{l}
{m}
"""
    # For example:
    #LaAlO3_SrTiO3_layer.PDOS
    #LAO_STO_pdos_Ti_3d.dat
    #Ti
    #3
    #2
    #9
    if os.path.exists(outfile):
        os.remove(outfile)
    with open('pdos_tmp_input.txt', 'w') as myfile:
        myfile.write(inp)
    os.system("fmpdos < pdos_tmp_input.txt")
    efermi=read_efermi(pdos_fname)
    return outfile, efermi

#core function
def gen_pdos_figure(pdos_fname, iatom, n, l, m, xlim=(-10, 10), ylim=(0,None)):
    outfile, efermi = get_pdos_data(pdos_fname, iatom, n, l, m)
    figname=f"pdos_{iatom}_{n}{lmdict[l][m]}.png"
    plot_pdos(fname=outfile, figname=figname,efermi=efermi, xlim=xlim, ylim=ylim)


def plot_pdos_ax(fname, efermi, ax=None, conv_n=1, xlim=(-10,10), ylim=(0, None)):
    data=np.loadtxt(fname)
    plt.rc('font', size=16)
    n=conv_n #为了pdos线更平滑
    data[:,1]=np.convolve(data[:,1], np.array([1.0/n]*n),mode='same') #convolution process 
    #d=np.convolve(data[:,1], np.array([1.0/n]*n),mode='same')[:-4] #convolution process 
    ax.plot(data[:,0]-efermi, data[:,1],label=fname)
    #plt.plot(data[4:,0]-efermi, d,label=fname)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    #plt.ylim(0, 15 )
    ax.axvline(color='red') 
    ax.set_xlabel('Energy (eV)') 
    #ax.set_ylabel('DOS')
    #plt.title(figname)
    #plt.tight_layout()
    #plt.savefig(figname)
    #plt.show()
    #plt.close() #plt.show() have a function of close. plt.close() means close the figure.

def plot_pdos(fname, figname, efermi, xlim=(-10, 10), ylim=(0,None)):
    fig, ax=plt.subplots()
    plot_pdos_ax(fname, efermi, ax=ax, xlim=xlim, ylim=ylim)
    plt.title(figname)
    plt.tight_layout()
    plt.savefig(figname)
    #plt.show()
    plt.close() 

def plot_layer_pdos(pdos_fname, figname, iatoms, n, l, m, xlim=(-10, 10), ylim=(0,None)):
    natoms=len(iatoms)
    fig, axes=plt.subplots(natoms,1, sharex=True)
    for i, iatom in enumerate(iatoms):
        outfile, efermi=get_pdos_data(pdos_fname, iatom, n, l, m)
        plot_pdos_ax(outfile, efermi, ax=axes[i], conv_n=5, xlim=xlim, ylim=ylim)
    plt.subplots_adjust(hspace=0.05)
    plt.savefig(figname)

    plt.show()


####111111
#plot all single atoms dos
#for iatom in range(16,17):
#for iatom in [16]:
#    # atom all orbital
#    gen_pdos_figure('LaAlO3.PDOS',iatom=iatom, n=0, l=-1, m=9, ylim=(0,None))
#    #3s (m=0 not 9)
#    gen_pdos_figure('LaAlO3.PDOS',iatom=iatom, n=3, l=0, m=0, ylim=(0,None))
#    #3p
#    gen_pdos_figure('LaAlO3.PDOS',iatom=iatom, n=3, l=1, m=9, ylim=(0,None))
#    #3d-all orbitals
#    gen_pdos_figure('LaAlO3.PDOS',iatom=iatom, n=3, l=2, m=9, ylim=(0,None))

####111111
# plot single element
#gen_pdos_figure('LaAlO3_SrTiO3_layer.PDOS', 'Sr', 0,-1,9 ,ylim=(0,10))
#gen_pdos_figure('LaAlO3_SrTiO3_layer.PDOS', 'Ti', 0,-1,9 ,ylim=(0,40))
#gen_pdos_figure('LaAlO3_SrTiO3_layer.PDOS', 'O', 0,-1,9 )
#gen_pdos_figure('LaAlO3_SrTiO3_layer.PDOS', 'Al', 0,-1,9 ,ylim=(0,10))
#gen_pdos_figure('LaAlO3_SrTiO3_layer.PDOS', 'La', 0,-1,9 )
#gen_pdos_figure('LaAlO3_SrTiO3_layer.PDOS', 'Ti', 3,2,9 )

####111111
#for m in range(-2,3):
#    gen_pdos_figure('LaAlO3_SrTiO3_layer.PDOS', 'Ti', 3,2,m )
#plt.title('Ti_3d.png')
#plt.tight_layout()
#plt.legend(fontsize=10)
#plt.savefig('Ti_3d.png')

####222222 
#cubicLAO_cubicSTO dos of different layers
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_Ti.png',iatoms=[59, 56, 55, 60, 57, 58], n=0, l=-1, m=9)
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_Al.png',iatoms=[2,5,1,4,6,3], n=0, l=-1, m=9)
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_La.png',iatoms=[x + 6 for x in [3,5,2,1,6,4]], n=0, l=-1, m=9)
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_Sr.png',iatoms=[x + 48 for x in [3,2,6,5,4,1]], n=0, l=-1, m=9)
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_O_inplane.png',iatoms=[x + 12 for x in [21,24,27,35,31,19,4,3,8,10,13,28]], n=0, l=-1, m=9)
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_O_outOfplane.png',iatoms=[x + 12 for x in [20,23,26,29,36,32,5,2,16,9,12,15]], n=0, l=-1, m=9)

# R3c-LAO_cubicSTO dos of different layers
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_Ti.png',iatoms=[118,117,113,111,109,115], n=0, l=-1, m=9,ylim=(0,2))
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_Al.png',iatoms=[7,5,8,2,11,10], n=0, l=-1, m=9,ylim=(0,0.5))
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_Sr.png',iatoms=[108,100,105,97,101,103], n=0, l=-1, m=9,ylim=(0,1))
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_O_inplane.png',iatoms=[72,69,82,96,86,37,29,27,57,49,45,60], n=0, l=-1, m=9,ylim=(0,1.2))
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_O_outOfplane.png',iatoms=[79,78,83,66,93,92,25,36,51,56,39,44], n=0, l=-1, m=9,ylim=(0,2))
#plot_layer_pdos('LaAlO3_SrTiO3_layer.PDOS','layeredDOS_La.png',iatoms=[22,14,24,13,20,16], n=0, l=-1, m=9,ylim=(0,2.5))
