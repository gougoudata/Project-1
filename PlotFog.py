
from jasp import *
from enthought.mayavi import mlab
from ase.data import vdw_radii
from ase.data.colors import cpk_colors
from ase.data.molecules import molecule
from ase import Atom, Atoms

# this module plots charge density using volume rendering
#inputs:

# path - location of the DFT files calculated by jasp. (Its input is mandatory)

# UnitCell - determines whether the unit cell is to be plotted. Default is 'False'.

# opacity - contains the minimum and maximum values of the opacity. Default is [0,1].

# view - contains the azimuth and the elevation. It is used to adjust the view. Default is [-90,90].

# SafeFig - Determines whether or not to save the file.

# filepath - stores the path where the image is to be saved. Default is 'images/fog.png'.

def plot_charge_density_fog(path, UnitCell=False, opacity=[0,1], view=[-90,90],SaveFig=True, filepath='images/fog.png'):


    with jasp(path) as calc:
         atoms = calc.get_atoms()
         x, y, z, cd = calc.get_charge_density()

    # making a black background
    mlab.figure(bgcolor=(0, 0, 0))

    # plotting the atoms as spheres,
    for atom in atoms:
        mlab.points3d(atom.x,
                      atom.y,
                      atom.z,
                      scale_factor=vdw_radii[atom.number]/5.,
                      resolution=20,
                      # a tuple is required for the color
                      color=tuple(cpk_colors[atom.number]),
                      scale_mode='none')
    if UnitCell==True:
    # draw the unit cell - there are 8 corners, and 12 connections
        a1, a2, a3 = atoms.get_cell()
        origin = [0, 0, 0]
        cell_matrix = [[origin,  a1],
                       [origin,  a2],
                       [origin,  a3],
                       [a1,      a1 + a2],
                       [a1,      a1 + a3],
                       [a2,      a2 + a1],
                       [a2,      a2 + a3],
                       [a3,      a1 + a3],
                       [a3,      a2 + a3],
                       [a1 + a2, a1 + a2 + a3],
                       [a2 + a3, a1 + a2 + a3],
                       [a1 + a3, a1 + a3 + a2]]

        for p1, p2 in cell_matrix:
            mlab.plot3d([p1[0], p2[0]], # x-positions
                        [p1[1], p2[1]], # y-positions
                        [p1[2], p2[2]], # z-positions
                        tube_radius=0.02)

    # plotting the charge density
    source = mlab.pipeline.scalar_field(x,y,z,cd)
    vmin=cd.min()
    vmax=cd.max()


    vol = mlab.pipeline.volume(source)
 
    # Changing the otf
    from enthought.tvtk.util.ctf import PiecewiseFunction
    otf = PiecewiseFunction()
    otf.add_point(vmin+0.1*(vmax-vmin),0)
    otf.add_point(vmin+0.8*(vmax-vmin),0.8)
    vol._otf = otf
    vol._volume_property.set_scalar_opacity(otf)

    # view adjusted by iteration

    mlab.view(azimuth= view[0], elevation=view[1], distance='auto')

    if SaveFig== True:
        mlab.savefig(filepath)
        mlab.show()
    
