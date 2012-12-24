from quaternion import *
from eemb import Monomer,AtomicCoordinate

def test():
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
    m.efp_rhf = fitQMWaterToEFP1(m,f_o1,f_h2,f_h3)
    pass

def fitQMWaterToEFP1_DFT(m):
    # EFP-1 H2ODFT Definition
    f_o1 = ( 1.1208360992, -3.5600913723,  1.2297518475 )
    f_h2 = ( 0.1818649382, -3.6467321115,  1.2710434338 )
    f_h3 = ( 1.3397706750, -3.4090717327,  0.3241364909 )
    m.efp_dft = fitQMWaterToEFP1(m,f_o1,f_h2,f_h3)
    pass

def fitQMWaterToEFP1(m,f_o1,f_h2,f_h3):

    e_o1 = m.atoms[0].xyz()
    e_h2 = m.atoms[1].xyz()
    e_h3 = m.atoms[2].xyz()

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

if __name__ == '__main__':
   o1 = AtomicCoordinate( 8.0, 1.1208360992, -3.5600913723,  1.2297518475, 'O')
   h2 = AtomicCoordinate( 1.0, 0.1818649382, -3.6467321115,  1.2710434338, 'H')
   h3 = AtomicCoordinate( 1.0, 1.3397706750, -3.4090717327,  0.3241364909, 'H')
   testMonomer = Monomer("test h2o monomer", [ o1, h2, h3 ])
   fitQMWaterToEFP1_RHF(testMonomer)
   print testMonomer.coordinates()
   print testMonomer.efp_rhf_coordinates()
   
