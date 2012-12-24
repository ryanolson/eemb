from decimal import Decimal
from copy import copy, deepcopy
from quaternion import *

import re

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
    def distanceBetween(self,other):
        return v_mag( v_sub( self.xyz(), other.xyz() ) )

class Monomer(object):
    def __init__(self, description, atoms):
        self.description = description
        self.atoms = deepcopy(atoms)
        if(self.atoms):
           self.efp_rhf = fitQMWaterToEFP1_RHF(self.atoms)
           self.efp_dft = fitQMWaterToEFP1_DFT(self.atoms)
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

 

def testEFPFit():
    o1 = AtomicCoordinate( 8.0, 1.1208360992, -3.5600913723,  1.2297518475, 'O')
    h2 = AtomicCoordinate( 1.0, 0.1818649382, -3.6467321115,  1.2710434338, 'H')
    h3 = AtomicCoordinate( 1.0, 1.3397706750, -3.4090717327,  0.3241364909, 'H')
    testMonomer = Monomer("test h2o monomer", [ o1, h2, h3 ])
    fitQMWaterToEFP1(testMonomer)

def fitQMWaterToEFP1_RHF(m):
    # EFP-1 H2ORHF Definition
    f_o1 = ( 0.17518,    4.35630,   -0.65086)
    f_h2 = (-0.83927,    3.98667,   -0.44045)
    f_h3 = ( 0.18512,    4.86728,   -1.62492)
    return fitQMWaterToEFP1(m,f_o1,f_h2,f_h3)

def fitQMWaterToEFP1_DFT(m):
    # EFP-1 H2ODFT Definition
    f_o1 = ( 1.1208360992, -3.5600913723,  1.2297518475 )
    f_h2 = ( 0.1818649382, -3.6467321115,  1.2710434338 )
    f_h3 = ( 1.3397706750, -3.4090717327,  0.3241364909 )
    return fitQMWaterToEFP1(m,f_o1,f_h2,f_h3)

def fitQMWaterToEFP1(qm,f_o1,f_h2,f_h3):

    e_o1 = qm[0].xyz()
    e_h2 = qm[1].xyz()
    e_h3 = qm[2].xyz()

#   print "QM Water"
#   print e_o1
#   print e_h2
#   print e_h3

    e_oh1 = v_sub(e_h2, e_o1)
    e_oh2 = v_sub(e_h3, e_o1)

#   print "magnitude of e_oh1 = ", v_mag(e_oh1)
#   print "magnitude of e_oh2 = ", v_mag(e_oh2)

    e_oh1 = normalize(e_oh1)
    e_oh2 = normalize(e_oh2)

#   print "e_oh1 = ", e_oh1
#   print "e_oh2 = ", e_oh2

    e_normal = v_cross( e_oh1, e_oh2 )
#   print "efp normal = ", e_normal

    e_hoh_angle = v_angle( e_oh1, e_oh2 )
#   print "efp h2-o1-h3 angle = ",e_hoh_angle
    q_test = axisangle_to_q( e_normal, e_hoh_angle )
    q_bisect = axisangle_to_q( e_normal, e_hoh_angle/2 )

    # if e_oh1 rotated by e_hoh_angle is equivalent to e_oh2
    # then rotate e_oh1 by half the e_hoh_angle to get the
    # bisecting vector; otherwise rotate e_oh2
    e_test = qv_mult( q_test, e_oh1 )
    e_bisect = None
    if v_AlmostEqual( e_test, e_oh2 ):
       e_bisect = qv_mult( q_bisect, e_oh1 )
    else:
       print "bisect test failed"
       return
#   print "e_bisect = ", e_bisect

#   print "\n molecule to be fit ..."
#   print f_o1
#   print f_h2
#   print f_h3

    f_oh1 = v_sub(f_h2, f_o1)
    f_oh2 = v_sub(f_h3, f_o1)

#   print "f_oh1 = ", f_oh1
#   print "f_oh2 = ", f_oh2

    f_normal = v_cross( f_oh1, f_oh2 )
    n_angle  = v_angle( f_normal, e_normal )
    if(abs(round(n_angle,6)) != 0):
       n_normal = v_cross( f_normal, e_normal )
       q_rotate = axisangle_to_q( n_normal, n_angle )

       f_test = qv_mult( q_rotate, f_normal )
       if( v_AlmostEqual( normalize(f_test), normalize(e_normal) )):
          f_normal = f_test
       else:
          print "failed to properly rotate normal vector"
          return

       f_oh1 = qv_mult( q_rotate, f_oh1 )
       f_oh2 = qv_mult( q_rotate, f_oh2 )


    f_hoh_angle = v_angle( f_oh1, f_oh2 )
#   print "fit h2-o1-h3 angle = ",f_hoh_angle

    q_test = axisangle_to_q( f_normal, f_hoh_angle )
    q_bisect = axisangle_to_q( f_normal, f_hoh_angle/2 )

    # if e_oh1 rotated by e_hoh_angle is equivalent to e_oh2
    # then rotate e_oh1 by half the e_hoh_angle to get the
    # bisecting vector; otherwise rotate e_oh2
    f_oh2n = f_oh2
    f_test = qv_mult( q_test, f_oh1 )
    f_bisect = None
    if v_AlmostEqual( normalize(f_test), normalize(f_oh2n) ):
       f_bisect = qv_mult( q_bisect, f_oh1 ) 
    else:
       print "bisect test failed"
       return

    bisect_normal = v_cross( f_bisect, e_bisect )
    bisect_angle  = v_angle( f_bisect, e_bisect )
    q_rotate = axisangle_to_q( bisect_normal, bisect_angle )
#   print "bisect_angle = ", bisect_angle
   
    f_test = qv_mult(q_rotate,f_bisect)
    if not v_AlmostEqual( normalize(f_test), normalize(e_bisect) ):
       print "rotated f_bisect = ", f_test
       print "bisect rotation failed"
       return

    f_oh1 = qv_mult( q_rotate, f_oh1 )
    f_oh2 = qv_mult( q_rotate, f_oh2 )

#   print "f_oh1 = ", f_oh1
#   print "f_oh2 = ", f_oh2

    f_o1 = e_o1
    f_h2 = v_add(f_oh1,e_o1)
    f_h3 = v_add(f_oh2,e_o1)

    o1 = AtomicCoordinate( 8.0, f_o1[0], f_o1[1], f_o1[2], 'O1')
    h2 = AtomicCoordinate( 1.0, f_h2[0], f_h2[1], f_h2[2], 'H2')
    h3 = AtomicCoordinate( 1.0, f_h3[0], f_h3[1], f_h3[2], 'H3')

    return [ o1, h2, h3 ]

