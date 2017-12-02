#!/usr/bin/env python

import os


import numpy as np


from geo_methods import FoldOnThisBox
from read_snapshots import Atom, read_xyz_snapshot, FFormat
from Molecule import Molecule
from getHB import iter_snapshot
from composed import Arrow
from Two_D_classes import Triangle
from vmd_resources import Colors, Materials, Molecule as DrawMol
from hbstrength import Naaas, wrap_strengths


def center_coords(atoms, box):
    center = np.zeros(3)
    for atom in atoms:
        center += np.array(atom.pos)
    center /= len(atoms)
    for i in range(len(atoms)):
        pass
        atoms[i].pos = atoms[i].pos - center + np.array(box) * 0.5
    return atoms


def atomize(atomlist):
    atoms = []
    i = 0
    for atom in atomlist:
        newatom = Atom()
        newatom.name = atom[0]
        newatom.element = atom[0]
        newatom.pos = atom[1:]
        newatom.index = i
        i += 1
        atoms.append(newatom)
    return atoms


def put_in_mols(atoms):
    """
    assume water and thus group every three atoms
    :param atoms: a regular atomlist that make up molecules of one type
    :type: list
    :return: a list of molecules
    :rtype: list
    """
    assert(len(atoms) % 3 == 0)
    molecules = []
    for i in range(0, len(atoms), 3):
        mol = Molecule()
        for j in range(0, 3, 1):
            mol.append(atoms[i + j])
        molecules.append(mol)
    return molecules


def energyradius(energy):
    """
    nasiges
    :param energy:
    :return:
    """
    min = 0
    max = -0.02
    maxrad = 0.75
    if energy < max:
        return maxrad
    elif energy > min:
        return 0
    else:
        return (energy / max) * maxrad


def writearrows(hbonds, outfile, howtofold):
    """
    write all hbonds a arrays to file
    :param hbonds: list of hydrogen bonds
    :type hbonds: HydrogenBonds
    :param outfile: the output file
    :param howtofold: function to fold coordinates
    """
    for hb in hbonds:
        distvec = howtofold(np.array(hb.acceptoratom.pos) - np.array(hb.donoratom.pos))
        radius = energyradius(hb.energy)
        arrow = Arrow(hb.donoratom.pos, np.array(hb.donoratom.pos) + distvec, radius)
        arrow.color = Colors.yellow
        arrow.material = Materials.GlassBubble
        arrow.draw(outfile)
    wfh.write('\n')


def write_proper_coords(atoms, foutput):
    """
    write folded coords back to file
    :param atoms: the list of atoms
    :type atoms: Atom
    :param foutput:
    """
    with open(foutput, 'w') as wfh:
        wfh.write(" {0:d}\n".format(len(atoms)))
        wfh.write("\n")
        for atom in atoms:
            wfh.write(atom.to_file_string(FFormat.xyz))


def add_surface(atoms, wfh):
    """
    write instructions for a surface plane
    :param atoms: list of the atom slabs
    :param wfh: the target file handle
    :type wfh: file
    """
#   get position of highest oxygen
    ymax = 0
    ymin = 1e10
    yoffset = 3
    for atom in atoms:
        if atom.pos[1] > ymax:
            ymax = atom.pos[1]
        if atom.pos[1] < ymin:
            ymin = atom.pos[1]
#   write surface draw instructions
    drawmol = DrawMol(is_new=False, index=1)
    drawmol.material = Materials.Transparent
    drawmol.color = Colors.blue
    drawmol.add_object(Triangle, [100, ymax - yoffset, -20], [-20, ymax - yoffset, 100], [-20, ymax - yoffset, -20])
    drawmol.add_object(Triangle, [100, ymax - yoffset, -20], [-20, ymax - yoffset, 100], [100, ymax - yoffset, 100])
    drawmol.add_object(Triangle, [100, ymin + yoffset, -20], [-20, ymin + yoffset, 100], [100, ymin + yoffset, 100])
    drawmol.add_object(Triangle, [100, ymin + yoffset, -20], [-20, ymin + yoffset, 100], [-20, ymin + yoffset, -20])
    wfh.write(str(drawmol))
    wfh.write('\n')


if __name__ == "__main__":
    geofile = "surface.250000.xyz"
    foldedfile = geofile[:len(geofile) - 4] + "-fold.xyz"
    with open(geofile, "r") as snapshot_fh:
        atomslist = read_xyz_snapshot(snapshot_fh)
    atoms = atomize(atomslist)
    box = np.loadtxt(geofile[:len(geofile) - 4] + ".c")
    atoms = center_coords(atoms, box)
    folder = FoldOnThisBox(box)
    mols = put_in_mols(atoms)
    hbs = iter_snapshot(mols, howtofold=folder.fold)
    print("Found {0:d} hydrogen bonds".format(len(hbs)))
    write_proper_coords(atoms, foldedfile)
    energytable = Naaas(len(mols))
    energytable.create_table("./molecules.lowest.donor", "./molecules.lowest.donor.label")
    wrap_strengths(hbs, energytable)
    with open("./vmd-all-hb.txt", "w") as wfh:
        wfh.write("color Display Background white\n")
        wfh.write("mol new {0:s}\n".format(os.getcwd().replace("\\", "/") + "/" + foldedfile))
        wfh.write("pbc set {{{0:f} {1:f} {2:f}}}\n".format(*box))  # set for both mols
        wfh.write("mol top 0\n")
        wfh.write("pbc set {{{0:f} {1:f} {2:f}}}\n".format(*box))
        wfh.write("draw material GlassBubble\n")
        writearrows(hbs, wfh, howtofold=folder.fold)
        add_surface(atoms, wfh)
