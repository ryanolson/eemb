from Queue import Queue, PriorityQueue
from threading import Thread
from decimal import Decimal
from copy import copy, deepcopy
from launcher import eemb_launch
import argparse

import os
import re
import time
import itertools


class AtomicCoordinate(object):
    def __init__(self, znuc, x, y, z, label=None):
        self.znuc = int(znuc)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.l = str(label)
        return
    def label(self):
        if self.l:
           return self.l
        return "RMO"
    def xyz(self):
        return (self.x, self.y, self.z)
    def coordinate_string(self):
        return "{l}    {zn}.0   {x:15.10f} {y:15.10f} {z:15.10f}".format(l=self.label(), zn=self.znuc, x=self.x, y=self.y, z=self.z)
    def coordinate_without_znuc(self):
        return "{l}    {x:15.10f} {y:15.10f} {z:15.10f}".format(l=self.label(), x=self.x, y=self.y, z=self.z)

class Monomer(object):
    def __init__(self, description, atoms):
        self.description = description
        self.atoms = deepcopy(atoms)
        if(self.atoms):
           self.efp_rhf = fitQMWaterToEFP1_RHF(self)
           self.efp_dft = fitQMWaterToEFP1_DFT(self)
        return
    def print_coordinates(self,array):
        s = ""
        for a in array:
            if s != "":
               s += "\n"
            s += a.coordinate_string()
        return s
    def print_coordinates_without_znuc(self,array):
        s = ""
        for a in array:
            if s != "":
               s += "\n"
            s += a.coordinate_without_znuc()
        return s
    def coordinates(self):
        return self.print_coordinates(self.atoms)
    def fit_efp1_rhf(self):
        pass
    def fit_efp1_dft(self):
        pass
    def efp_dft_coordinates(self):
        return self.print_coordinates_without_znuc(self.efp_dft)
    def efp_rhf_coordinates(self):
        return self.print_coordinates_without_znuc(self.efp_rhf)


def readCoordinateFile(coord_file):
    m = [ ] 
    monomers = [ ]
    wb = re.compile(r'\s+')
    with open(coord_file, "r") as f:
       for line in f:
           line = wb.sub(' ', line.strip()) 
           (label, znuc, x, y, z) = line.split(' ')
           m.append( AtomicCoordinate(int(Decimal(znuc)), x, y, z, label) )
           if len(m) == 3:
              monomers.append( Monomer("water {0}".format(len(monomers)+1), m) )
              m = [ ]
    return monomers
