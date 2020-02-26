import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
%matplotlib inline

import math
# Constants
G = 6.67408e-11
g = 9.81

class NBody(object):
    '''
        This is a class for a Body given an initial position, velocity and mass
    '''
    def __init__(self, x0, y0, z0, vx0, vy0, vz0, mass):
        # Initial position
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        # Initial velocity 
        self.vx0 = vx0
        self.vy0 = vy0
        self.vz0 = vz0
        # Mass of the body
        self.mass = mass
        
def calculateMod(a, b, c):
    """ 
    Args:
       a : first number
       b : second number
       c : third number

    Returns:
        Modulus of three numbers
    """
    return math.sqrt(a**2 + b**2 + c**2)

def multiplyTuple(scalar, tupleInput):
    """ 
    Args:
       scalar: scalar value
       tupleInput: tuple

    Returns:
        Tuple 
     """
    return (scalar * tupleInput[0], scalar * tupleInput[1], scalar * tupleInput[2])

def f(x, t):
    return x**2 + 1
    
def solveEquation(arrayBodies, a, b, h, method):
    """ 
    Args:
       arrayBodies : array of bodies with first element being focused on
       a : lower limit
       b : upper limit
       h : timestep
       method: 'euler', 'rk2' or 'rk4'
       f : function, right hand side of the ODE

    Returns:
        A tuple with first element as time points and second element as position points corresponding to the time
    """
    # Initialize positions and velocities at their input values as a tuple
    initialPosition = (arrayBodies[0].x0, arrayBodies[0].y0, arrayBodies[0].z0)
    initialVelocity = (arrayBodies[0].vx0, arrayBodies[0].vy0, arrayBodies[0].vz0)
    
    
    # Save outputs in an array of tuples
    positionSave = [(arrayBodies[0].x0, arrayBodies[0].y0, arrayBodies[0].z0)]
    velocitySave = [(arrayBodies[0].vx0, arrayBodies[0].vy0, arrayBodies[0].vz0)]
    
    # Calculate acceleration as a tuple on the body being focused at
    numBodies = len(arrayBodies)
    accel = (0, 0, 0)
    for i in range(numBodies - 1):
        modVal1 = arrayBodies[i+1].x0
        modVal2 = arrayBodies[i+1].y0
        modVal3 = arrayBodies[i+1].z0
        firstVal = (G*arrayBodies[i+1].mass)/((calculateMod(modVal1 - initialPosition[0], modVal2 - initialPosition[1], modVal3 - initialPosition[2]))**3)
        aTuple = multiplyTuple(firstVal, (modVal1 - initialPosition[0], modVal2 - initialPosition[1], modVal3 - initialPosition[2]))
        accel = (accel[0] + aTuple[0], accel[1] + aTuple[1], accel[2] + aTuple[2])
        
    tpoints = np.arange(a, b, h)
    posPoints = []
    r0 = (0, 0, 0)
    for t in tpoints:
        posPoints.append(r0)
        # Euler's Method
        if (method == 'euler'):
            r0 = (r0[0] + h*f(r0[0], t), r0[1] + h*f(r0[1], t), r0[2] + h*f(r0[2], t))
        # Second Order Runge-Kutta Method 
        elif (method == 'rk2'):
            k1 = (h*f(r0[0], t), h*f(r0[1], t), h*f(r0[2], t))
            k2 = (h*f(r0[0] + 0.5*k1[0], t + 0.5*h), h*f(r0[1] + 0.5*k1[1], t + 0.5*h), h*f(r0[2] + 0.5*k1[2], t + 0.5*h))
            r0 = (r0[0] + k2[0], r0[1] + k2[1], r0[2] + k2[2])
        # Fourth Order Runge-Kutta Method
        elif (method == 'rk4'):
            k1 = (h*f(r0[0], t), h*f(r0[1], t), h*f(r0[2], t))
            k2 = (h*f(r0[0] + 0.5*k1[0], t + 0.5*h), h*f(r0[1] + 0.5*k1[1], t + 0.5*h), h*f(r0[2] + 0.5*k1[2], t + 0.5*h))
            k3 = (h*f(r0[0] + 0.5*k2[0], t + 0.5*h), h*f(r0[1] + 0.5*k2[1], t + 0.5*h), h*f(r0[2] + 0.5*k2[2], t + 0.5*h))
            k4 = (h*f(r0[0] + k3[0], t + h), h*f(r0[1] + k3[1], t + h), h*f(r0[2] + k3[2], t + h))
            r0 = (r0[0] + (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])/6, r0[1] + (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])/6, r0[2] + (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])/6)
    return (tpoints, posPoints)

# t = time calculated above
t = 0
Earth = NBody(0, 0, 0, 0, 0, 0, 5.972e24)
Projectile = NBody(0, 0, 10, 0, 0, 12000, 0.1)
bodyArray = [Projectile, Earth]
PARTC = solveEquation(bodyArray, 0, t, 0.01, 'euler')
plt.plot(PARTC[0], PARTC[1])

solveEquation(bodyArray, 0, 30000, 100, 'euler')
solveEquation(bodyArray, 0, 30000, 500, 'euler')

# To get the error, calculate using a different timestep and compare the results
solveEquation(bodyArray, 0, 30000, 500, 'rk4')
solveEquation(bodyArray, 0, 30000, 400, 'rk4')

%%timeit  
solveEquation(bodyArray, 0, 30000, 100, 'euler')

%%timeit
solveEquation(bodyArray, 0, 30000, 100, 'rk4')

Moon = NBody(0, 0, 384403, 0, 0, 1.022, 7.348e22)
bodyArray2 = [Moon, Earth]
PARTG1 = solveEquation(bodyArray2, 0, 30000, 100, 'euler')
plt.plot(PARTG1[0], PARTG1[1])

bodyArray2 = [Earth, Moon]
PARTG2 = solveEquation(bodyArray2, 0, 30000, 100, 'euler')
plt.plot(PARTG2[0], PARTG2[1])

STAR1 = (-0.5 * (1.49597871e11), -0.28867 * (1.49597871e11), 0, 14894.68, -25798.34, 0, 'euler')
STAR2 = (0, -0.57735 * (1.49597871e11), 0, -29789.36, 0, 0, 'euler')
STAR3 = STAR1 = (0.5 * (1.49597871e11), -0.28867 * (1.49597871e11), 0, 14894.68, 25798.34, 0, 'euler')

BodyArray4a = [STAR1, STAR2, STAR3]
BodyArray4b = [STAR2, STAR1, STAR3]
BodyArray4c = [STAR3, STAR1, STAR2]

%%timeit
PARTHA = solveEquation(BodyArray4a, 0, 30000, 10000, 'euler')

%%timeit
PARTHB = solveEquation(BodyArray4b, 0, 30000, 10000, 'euler')

%%timeit
PARTHC = solveEquation(BodyArray4c, 0, 30000, 10000, 'euler')


plt.plot(PARTHA[0], PARTHA[1])
plt.plot(PARTHB[0], PARTHB[1])
plt.plot(PARTHC[0], PARTHC[1])

%%timeit
PARTHA = solveEquation(BodyArray4a, 0, 30000, 10000, 'euler')

%%timeit
PARTHB = solveEquation(BodyArray4b, 0, 30000, 10000, 'euler')

%%timeit
PARTHC = solveEquation(BodyArray4c, 0, 30000, 10000, 'euler')


plt.plot(PARTHA[0], PARTHA[1])
plt.plot(PARTHB[0], PARTHB[1])
plt.plot(PARTHC[0], PARTHC[1])

%%timeit
PARTHA = solveEquation(BodyArray4a, 0, 30000, 10000, 'rk4')

%%timeit
PARTHB = solveEquation(BodyArray4b, 0, 30000, 10000, 'rk4')

%%timeit
PARTHC = solveEquation(BodyArray4c, 0, 30000, 10000, 'rk4')


plt.plot(PARTHA[0], PARTHA[1])
plt.plot(PARTHB[0], PARTHB[1])
plt.plot(PARTHC[0], PARTHC[1])