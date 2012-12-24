from math import sqrt, cos, sin, acos

def normalize(v, tolerance=0.00001):
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        mag = sqrt(mag2)
        v = tuple(n / mag for n in v)
    return v

def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return w, x, y, z

def q_conjugate(q):
    q = normalize(q)
    w, x, y, z = q
    return (w, -x, -y, -z)

def qv_mult(q1, v1):
    mag = v_mag(v1)
    v1 = normalize(v1)
    q2 = (0.0,) + v1
    x,y,z = q_mult(q_mult(q1, q2), q_conjugate(q1))[1:]
    return x*mag, y*mag, z*mag

def axisangle_to_q(v, theta):
    v = normalize(v)
    x, y, z = v
    theta /= 2
    w = cos(theta)
    x = x * sin(theta)
    y = y * sin(theta)
    z = z * sin(theta)
    return w, x, y, z

def q_to_axisangle(q):
    w, v = q[0], q[1:]
    theta = acos(w) * 2.0
    return normalize(v), theta

def v_sub(v1,v2):
    a0, a1, a2 = v1
    b0, b1, b2 = v2
    return a0-b0, a1-b1, a2-b2  

def v_add(v1,v2):
    a0, a1, a2 = v1
    b0, b1, b2 = v2
    return a0+b0, a1+b1, a2+b2  

def v_dot(v1,v2):
    dot = sum(v1[i]*v2[i] for i in range(len(v1)))
    return dot

def v_mag(v):
    mag2 = sum(n * n for n in v)
    return sqrt(mag2)

def v_cross(v1,v2):
    a0, a1, a2 = v1
    b0, b1, b2 = v2
    c0 = a1*b2 - a2*b1
    c1 = a2*b0 - a0*b2
    c2 = a0*b1 - a1*b0
    return c0, c1, c2

def v_angle(v1,v2):
    val = v_dot(v1,v2)/(v_mag(v1)*v_mag(v2)) 
    if abs(round(val,8)) == 1: return 0
    return acos( v_dot(v1,v2) / ( v_mag(v1) * v_mag(v2) ) )

def v_AlmostEqual(v1,v2):
    v = v_sub(v1,v2)
    if round(v_mag(v),7) == 0:
       return True
    return False
