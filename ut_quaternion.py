import unittest
from math import pi
from quaternion import *

class QuaternionTest(unittest.TestCase):

    def testRotateVectorAroundTheXAxis(self):
        xaxis = (2,0,0)
        theta = pi/2
        vector = (0,2,0)
        q_rotateabout = axisangle_to_q(xaxis,theta)
        self.assertAlmostEqual( qv_mult(q_rotateabout, vector)[0], 0 )
        self.assertAlmostEqual( qv_mult(q_rotateabout, vector)[1], 0 )
        self.assertAlmostEqual( qv_mult(q_rotateabout, vector)[2], 2 )

    def testCrossProductUnitVectors(self):
        xaxis = (1,0,0)
        yaxis = (0,1,0)
        zaxis = (0,0,1)
        self.assertEqual( v_cross(xaxis,yaxis), zaxis )
        self.assertEqual( v_cross(yaxis,zaxis), xaxis )
        self.assertEqual( v_cross(zaxis,xaxis), yaxis )

    def testDotProduct(self):
        v1 = (1,2,3)
        v2 = (2,2,2)
        self.assertEqual( v_dot(v1,v2), 12 )

    def testVectorAngle(self):
        xaxis = (1,0,0)
        yaxis = (0,1,0)
        self.assertAlmostEqual(v_angle(xaxis,yaxis), pi/2)
        

if __name__ == '__main__':
    unittest.main()
