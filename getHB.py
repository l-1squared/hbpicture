#!/usr/bin/env python


import math


import numpy as np


from Molecule import Molecule
from read_snapshots import Atom
from geo_methods import vec_len, scalar_product


__author__ = """Kristof Karhan"""
"""
return a list of hydrogen bonded pairs
"""


class HydrogenBond(object):
    """
    store all information about a pair
    """
    def __init__(self, howtofold=None):
#       pointless summary of class members
        self.__donor = Molecule()
        self.__acceptor = Molecule()
        self.__donatingatom = Atom()
        self.__acceptingatom = Atom()
        self.__donorid = 0
        self.__acceptorid = 0
        self.energy = 1
        if howtofold is None:
            self.howtofold = lambda a : a
        else:
            self.howtofold = howtofold

    def add_donor(self, mol, atom):
        """
        set the donor of the bond
        :param mol: the donor molecule
        :param atom: the actual donor atom
        :type atom: Atom, int
        """
        try:
            atom * 3
            self.__donorid = int(atom)
            self.__donatingatom = mol[atom]
        except TypeError:
            self.__donatingatom = atom
            self.__donorid = int(mol.atoms.index(atom))
        self.__donor = mol

    def add_acceptor(self, mol, atom):
        """
        set the acceptor of the Bond
        :param mol: the acceptor molecule
        :param atom: the actual acceptor atom
        :type atom: Atom, int
        """
        try:
            atom * 3 # test whether atom is a number
            self.__acceptorid = int(atom)
            self.__acceptingatom = mol[atom]
        except TypeError:
            self.__acceptingatom = atom
            self.__acceptorid = int(mol.index(atom))
        self.__acceptor = mol

    def __str__(self):
        string = "Hydrogen bond donated from {0:06d} to {1:06d} with length {2:.3f}"
        return string.format(self.__donor.id, self.__acceptor.id, self.length)

    def get_length(self):
        """
        calculate the distance between O and H
        :return: the hydrogen bond length
        """
        length = np.array(self.__acceptingatom.pos) - np.array(self.__donatingatom.pos)
        length = vec_len(self.howtofold(length))
        return length

    def get_donor(self):
        return self.__donor

    def get_acceptor(self):
        return self.__acceptor

    def get_donoratom_id(self):
        return self.__donorid

    def get_acceptoratom_id(self):
        return self.__acceptorid

    def get_donor_atom(self):
        return self.__donatingatom

    def get_acceptor_atom(self):
        return self.__acceptingatom

    length = property(get_length)
    donor = property(get_donor)
    acceptor = property(get_acceptor)
    donorid = property(get_donoratom_id)
    acceptorid = property(get_acceptoratom_id)
    donoratom = property(get_donor_atom)
    acceptoratom = property(get_acceptor_atom)




class CalcHydrogenBond(object):
    """
    set up and store cutoff criteria for geometric hydrogen bonds
    """
    __ANGLES = [3, 9, 13, 15] + [17] * 2 + [19, 21] + [23] * 2 + [25] * 2 + [27] + [29] * 2 + [31] * 3\
               + [33] * 3 + [35] * 5 + [37] * 18 + [35, 33, 31] + [29] * 2 + [25, 21, 19, 15, 7, 3, 1]
    __DISTANCES = [(2.33e0 + 0.02e0 * i) for i in range(1, 57)]

    def __init__(self, custom_set=None, howtofold=None, donorlist=None, acceptorlist=None):
        """
        initialize a new CalcHydrogenBond instance
        :param custom_set: supply an alternative set of distance angle cutoffs
        :param howtofold: supply a folding routine
        :param donorlist: supply a list of atoms that can act as donors, default 'O', 'N', 'F'
        :param acceptorlist: supply a list of atoms that can act ass acceptors, default 'H'
        """
#       unnecessary overview of class variables
        self.__donorlist = []
        self.__acceptorlist = []
        self.__exclusion = []
        if howtofold is None:
            self.howtofold = lambda a : a
        else:
            self.howtofold = howtofold
        if custom_set is None:
            distances = CalcHydrogenBond.__DISTANCES
            angles = CalcHydrogenBond.__ANGLES
        else:
            assert(len(custom_set[0]) == len(custom_set[1]))
            distances = custom_set[0]
            angles = custom_set[1]

        if donorlist is None:
            self.__donorlist = ['O', 'N', 'F']
        else:
            self.__donorlist = donorlist

        if acceptorlist is None:
            self.__acceptorlist = ['H']
        else:
            self.__aceptorlist = acceptorlist
        self.__exclusion = np.zeros((len(distances), 2))
        for i in range(0, len(distances), 1):
            self.__exclusion[i, 0] = distances[i]
            self.__exclusion[i, 1] = angles[i]

    def is_hydrogen_bonded(self, distance, angle):
        """
        check if distance and angle criteria satisfy conditions
        :param distance: the distance between donor and anchor
        :param angle:  the angle \Theta(Nu-H, Nu-D)
        :return: whether there is a hydrogen bond
        """
        for cutoff in self.__exclusion:
            if distance <= cutoff[0]:
                return angle <= cutoff[1]
        return False

    def is_paired(self, donormol, acceptormol):
        """
        test for the existence of a hydrogen bond between two molecules
        :param donormol: the donor molecule
        :param acceptormol:  the acceptor molecules
        :return: all hydrogen bonds connecting the molecules
        """
        potentialdonors = []
        potentialacceptors = []
        hydrogen_pairs = []
        for atom in donormol:
            if atom.element in self.__donorlist:
                potentialdonors.append(atom)
        for atom in acceptormol:
            if atom.element in self.__acceptorlist:
                potentialacceptors.append(atom)
        for donoratom in potentialdonors:
            for acceptoratom in potentialacceptors:
                anchor = self.__anchor(acceptormol, acceptoratom)
                accdon = self.howtofold(np.array(donoratom.pos) - np.array(anchor.pos))
                distance = vec_len(accdon)
                nuh = self.howtofold(np.array(acceptoratom.pos) - np.array(anchor.pos))
                angle = math.acos(scalar_product(accdon, nuh)) * 180 / math.pi
                if self.is_hydrogen_bonded(distance, angle):
                    hb = HydrogenBond()
                    hb.add_donor(donormol, donoratom)
                    hb.add_acceptor(acceptormol, acceptoratom)
                    hydrogen_pairs.append(hb)
        if len(hydrogen_pairs) == 1:
            return hydrogen_pairs[0]
        else:
            return hydrogen_pairs

    def __anchor(self, mol, atom):
        """
        find the atom in a molecule closest to atom
        :param mol: the molecule containing atom
        :type mol: Molecule
        :param atom: a hydrogen atom
        :type atom: Atom
        :return: return the atom anchoring atom
        """
        theatom = None
        thedistance = 1e100
        for anatom in mol:
            if anatom == atom:
                continue
            distance = vec_len(self.howtofold(np.array(anatom.pos) - np.array(atom.pos)))
            if distance < thedistance:
                thedistance = distance
                theatom = anatom
        return theatom


def iter_snapshot(molecules, global_exclusion=10, **kwargs):
    """
    :param molecules: a list of all molecules in a snapshots
    :type molecules: list
    :param global_exclusion: at which com distance to abort H  bond search
    :param kwargs: optional arguments for CalcHydrogenBond
    :return:  all hydrogen bonds in a snapshot
    """
    all_hb = []
    calculator = CalcHydrogenBond(**kwargs)
    for i in range(len(molecules) - 1):
        for j in range(i + 1, len(molecules)):
            if vec_len(kwargs['howtofold'](molecules[i].com - molecules[j].com)) > global_exclusion:
                continue
            for k in range(2):  # k is not important, just run twice
                new_hbs = calculator.is_paired(molecules[i], molecules[j])
                try:
                    iter(new_hbs)  # I just don't want to enter extend to fail
                    all_hb.extend(new_hbs)
                except TypeError:
                    all_hb.append(new_hbs)
                i, j = j, i  # switch donors and acceptors

    return all_hb
