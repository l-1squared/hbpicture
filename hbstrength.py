#!/usr/bin/env python


import numpy as np


class Naaas(object):
    def __init__(self, natoms):
        self.array = np.zeros((natoms, natoms))

    def create_table(self, fenergies, flabels):
        energies = np.loadtxt(fenergies)
        labels = np.loadtxt(flabels)
        for row in labels:
            for i in range(1, len(row), 1):
                self.array[int(row[0]) - 1, int(row[i]) - 1] = energies[int(row[0]) - 1, i]


def assign_strengths(hb, table):
    donorid = hb.donor.id
    acceptorid = hb.acceptor.id
    hb.energy = table.array[donorid, acceptorid]
    return hb


def wrap_strengths(hbs, table):
    for i in range(len(hbs)):
        hbs[i] = assign_strengths(hbs[i], table)
    return hbs
