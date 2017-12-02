#!/usr/bin env python

__author__ = """Kristof Karhan"""


import numpy as np


from read_snapshots import FFormat


class Molecule(list):
    __GLOBALID__ = 0

    @staticmethod
    def __getid__():
        Molecule.__GLOBALID__ += 1
        return Molecule.__GLOBALID__ - 1

    def __init__(self):
        super(Molecule, self).__init__()
        self.__id = Molecule.__getid__()


    def printgeo(self, foutput, printformat=FFormat.xyz):
        if printformat != FFormat.xyz:
            raise NotImplementedError("lazily left out format")
        with open(foutput, 'w') as wfh:
            wfh.write(" {0:d}\n".format(len(self)))
            for anatom in self:
                wfh.write("{0:2s}    {1:12.6f} {2:12.6f} {3:12.6f}".format(anatom.element, *anatom.pos))

    def __len__(self):
        return super(Molecule, self).__len__()

    def get_atoms(self):
        return super(Molecule, self)

    def get_id(self):
        return self.__id

    def get_com(self):
        center = np.zeros(3)
        mass = 0
        for atom in self:
            pos = atom.pos
            mass += get_mass(atom.element)
            for i in range(3):
                center[i] += pos[i]
        center /= mass
        return center

    natoms = property(__len__)
    size = property(__len__)
    atoms = property(get_atoms)
    id = property(get_id)
    com = property(get_com)


def get_mass(element):
    element = element.strip()
    if element == 'H':
        return 1
    if element == "O":
        return 16
