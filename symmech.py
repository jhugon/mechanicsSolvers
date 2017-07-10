#!/usr/bin/env python

import sympy
from sympy import *

def difftotal(expr, diffby, diffmap):
  """Take the total derivative with respect to a variable.

  Example:

      theta, t, theta_dot = symbols("theta t theta_dot")
      difftotal(cos(theta), t, {theta: theta_dot})

  returns

      -theta_dot*sin(theta)
  
  Taken from http://robotfantastic.org/total-derivatives-in-sympy.html
  """
  # Replace all symbols in the diffmap by a functional form
  fnexpr = expr.subs({s:s(diffby) for s in diffmap})
  # Do the differentiation
  diffexpr = diff(fnexpr, diffby)
  # Replace the Derivatives with the variables in diffmap
  derivmap = {Derivative(v(diffby), diffby):dv 
              for v,dv in diffmap.iteritems()}
  finaldiff = diffexpr.subs(derivmap)
  # Replace the functional forms with their original form
  return finaldiff.subs({s(diffby):s for s in diffmap})

def generateTimeDerivMapFromCoords(*args):
  """
  Creates a dictionary containing the first and second 
  time derivatives for each coordinate. This is needed 
  for the makeMomentaHamiltonian and eulerLegrangeEOM 
  functions.
  """
  result = {}
  for arg in args:
    if type(arg)==list or type(arg)==tuple:
      tmpresult = generateTimeDerivMapFromCoords(*arg)
      for key in tmpresult:
        result[key] = tmpresult[key]
    else:
      coordName = str(arg)
      velocity = symbols(coordName+"_dot")
      acceleration = symbols(coordName+"_dot_dot")
      result[arg] = velocity
      result[velocity] = acceleration
  return result

def eulerLegrangeEOM(lagrangian,coordinates,velocities,timeVar,timeDerivMap):
  """
    Finds the Euler-Legrange equations of motion given 
    a lagrangian, the coordinates, the velocities, 
    and the dictionary of time derivatives.
    Ensure that for each index, i, velocities[i] 
    corresponds to the time derivitive of 
    coordinates[i]. For the total derivitive, the 
    timeDerivMap is used to substitute the proper 
    derivative variable.
  """

  if len(coordinates) != len(velocities):
    raise Exception("len(coordinates) != len(velocities)")

  result = []
  for coord, vel in zip(coordinates,velocities):
    dLdcoord = diff(lagrangian,coord)
    dLdvel = diff(lagrangian,vel)
    ddtdLdvel = difftotal(dLdvel,timeVar,timeDerivMap)
    eqn = Eq(dLdcoord,ddtdLdvel)
    result.append(eqn)
  return result

def makeMomentaHamiltonian(lagrangian,coordinates,velocities,timeDerivMap):
  """
    Finds the conjugate momenta for each of the 
    generalized coordiantes, and also the 
    hamiltonian in terms of the canonical 
    coordiantes. Returns a list of the momenta
    coordinate variables, a list of equations
    defining the momenta, and the hamiltonian.
  """
  if len(coordinates) != len(velocities):
    raise Exception("len(coordinates) != len(velocities)")

  momenta = [symbols("p_"+str(coord)) for coord in coordinates]
  for mom, coord in zip(momenta,coordinates):
    if not timeDerivMap.has_key(mom):
      timeDerivMap[mom] = symbols("p_dot_"+str(coord))
  momentaInTermsOfVel = [diff(lagrangian,vel) for vel in velocities]
  momentaDefs = [Eq(mom,momInTermsOfVel) for mom,momInTermsOfVel in zip(momenta,momentaInTermsOfVel)]
  velocitiesInTermsOfMom = []
  for mom, vel, momInTermsOfVel in zip(momenta,velocities,momentaInTermsOfVel):
    slnList = solve(Eq(mom,momInTermsOfVel),vel)
    if len(slnList) != 1:
      raise Exception("Wrong number of solutions inverting for velocity in terms of momentum, {0} solutions".format(len(slnList)))
    velocitiesInTermsOfMom.append(slnList[0])

  hamiltonian = -lagrangian
  for mom, vel in zip(momenta,velocities):
    hamiltonian += mom*vel
  hamiltonian = hamiltonian.subs({vel:velInTermsOfMom for vel,velInTermsOfMom in zip(velocities,velocitiesInTermsOfMom)})

  return momenta, momentaDefs, hamiltonian

def hamiltonianEOM(hamiltonian,coordinates,momenta,timeDerivMap):
  """
    Finds the Hamiltonian equations of motion given the Hamiltonian, coordinates,
    momenta, and dictionary of time derivatives.  Ensure that for each index, i,
    momenta[i] is the conjuegate momentum of coordinates[i].  For the total 
    derivitive, the timeDerivMap is used to substitute the proper derivative
    variable.
  """

  if len(coordinates) != len(momenta):
    raise Exception("len(coordinates) != len(momenta)")

  eom = []
  for coord, mom in zip(coordinates, momenta):
    coord_dot = timeDerivMap[coord]
    mom_dot = timeDerivMap[mom]
    dpdt = - diff(hamiltonian,coord)
    dqdt =   diff(hamiltonian,mom)
    eom.append(Eq(dpdt,mom_dot))
    eom.append(Eq(dqdt,coord_dot))

  return eom

class SymComputer(object):

  def __init__(self,lagrangian,coordinates,velocities,t=Symbol('t')):
    for coord in coordinates:
      if str(coord) == str(t):
        raise Exception("Coordinate '{0}' and time variable '{1}' have the same name".format(str(coord),str(t)))

    self.lagrangian = lagrangian
    self.coordinates = coordinates
    self.velocities = velocities
    self.t = t

    self.tdm = None

    self.momenta = None
    self.momentaDefs = None
    self.hamiltonian = None

    self.eomsL = None
    self.eomsH = None

    self.haveEOMsL = False
    self.haveEOMsH = False

    self.tdm = generateTimeDerivMapFromCoords(self.coordinates)


  def solveEulerLegrange(self):
    self.eomsL =  eulerLegrangeEOM(self.lagrangian,self.coordinates,self.velocities,self.t,self.tdm)
    self.haveEOMsL = True

  def solveHamiltonian(self):
    self.momenta, self.momentaDefs, self.hamiltonian = makeMomentaHamiltonian(self.lagrangian,
                                                        self.coordinates,self.velocities,self.tdm)
    self.eomsH = hamiltonianEOM(self.hamiltonian,self.coordinates,self.momenta,self.tdm)
    self.haveEOMsH = True

  def getEulerLegrangeEOMs(self):
    if not self.haveEOMsL:
      self.solveEulerLegrange()
    return self.eomsL

  def getHamiltonianEOMs(self):
    if not self.haveEOMsH:
      self.solveHamiltonian()
    return self.eomsH

  def __str__(self):
    result = ""

    result += "***SymComputer***"+"\n"
    result += "L = "+str(self.lagrangian)+"\n"
    if self.hamiltonian:
      result += "H = "+str(simplify(self.hamiltonian))+"\n"
    if self.momenta:
      result += "momenta:"+"\n"
      for mdef in self.momentaDefs:
        result += "  "+str(mdef)+"\n"
    if self.eomsL:
      result += "Euler-Legrange EOM:"+"\n"
      for eom in self.eomsL:
        result += "  "+str(simplify(eom))+"\n"
    if self.eomsH:
      result += "Hamiltonian EOM:"+"\n"
      for eom in self.eomsH:
        result += "  "+str(simplify(eom))+"\n"
    result += "*****************"+"\n"
    return result
    
class MatrixLagrangian(object):
    """
    Use linear algebra to manipulate a lagrangian using sympy
    """

    def __init__(self,kineticMatrix,velocityColumn,noVelocityTerms,coordinates,velocities):
        """
        You can use this class for Lagrangian of the form:
            L = 0.5 * (q_dot row vector) * (kinetic matrix) * (q_dot column vector) 
                + (q_dot row vector) * (velocity dependence column vector)
                + (terms not dependant on velocity)

        As described in Goldstein, Poole, and Safko, Classical Mechanics 3rd Ed. p339

        The kinetic matrix is the first argument, velocity dependence is the second,
            and all terms which don't depend on velocity should be put third. The final two
            arguments are a list of coordinates and a list of cooresponding velocities.

        velocityColumn may be None.
        """

        if not isinstance(kineticMatrix, Matrix):
            raise Exception("KineticMatrix {} is not of type Matrix".format(kineticMatrix))
        if not isinstance(velocityColumn, Matrix):
            if not (velocityColumn is None):
                raise Exception("velocityColumn {} is not of type Matrix".format(velocityColumn))

        self.kineticMatrix = kineticMatrix
        self.velocityColumn = velocityColumn
        self.noVelocityTerms = noVelocityTerms
        self.coordinates = coordinates
        self.velocities = velocities
        self.velocitiesColumn = Matrix(self.velocities)
        self.velocitiesRow = self.velocitiesColumn.T

    def getLagrangian(self):
        """
        Returns a sympy expression for the Lagrangian represented by this object
        """
        velocityDependence = 0
        if not (self.velocityColumn is None):
            velocityDependence = self.velocitiesRow * self.velocityColumn
        return self.velocitiesRow*self.kineticMatrix*self.velocitiesColumn/2 + velocityDependence + self.noVelocityTerms

    def getDiagonalizedKineticLagrangian(self):
        """
        Returns a sympy expression for the Lagrangian represented by this object
        """
        velocityDependence = 0
        eigen = self.kineticMatrix.eigenvects()

    def getInvertedKineticHamiltonian(self):
        """
        Returns a sympy expression for the Hamiltonian represented by this object, a list of the momenta, 
            and a list of equations defining the momenta.
        """

        momenta = [symbols("p_"+str(coord)) for coord in self.coordinates]
        momentaInTermsOfVel = self.kineticMatrix * self.velocitiesColumn + self.velocityColumn
        momentaDefs = [Eq(mom,momInTermsOfVel) for mom,momInTermsOfVel in zip(momenta,momentaInTermsOfVel)]

        kineticMomentaColumn = Matrix(momenta)
        if not (self.velocityColumn is None):
            kineticMomentaColumn -= self.velocityColumn
        invKinetic = self.kineticMatrix**(-1)
        hamiltonian = (kineticMomentaColumn.T)*invKinetic*(kineticMomentaColumn)/2 - self.noVelocityTerms
        return hamiltonian, momenta, momentaDefs

if __name__ == "__main__":

  g = symbols('g')
  m = symbols('m')
  M = symbols('M')
  l = symbols('l')
  R = symbols('R')
  I = symbols('I')
  K = symbols('K')
  F = symbols('F')
  tau = symbols('tau')
  lam = symbols('lambda')

  t = symbols('t')

  x = symbols('x')
  y = symbols('y')
  z = symbols('z')

  x_dot = symbols('x_dot')
  y_dot = symbols('y_dot')
  z_dot = symbols('z_dot')

  theta = symbols('theta')
  theta_dot = symbols('theta_dot')

  thetaB = symbols('thetaB')
  thetaB_dot = symbols('thetaB_dot')

  Ax = symbols('Ax')
  Ay = symbols('Ay')
  Az = symbols('Az')

  L = m * x_dot**2 / 2 + m * y_dot**2 / 2 + m * z_dot**2 / 2
  sc = SymComputer(L,[x,y,z],[x_dot,y_dot,z_dot])
  sc.getEulerLegrangeEOMs()
  sc.getHamiltonianEOMs()
  print sc
  
  L = m * x_dot**2 / 2 + F*x
  sc = SymComputer(L,[x],[x_dot])
  sc.getEulerLegrangeEOMs()
  sc.getHamiltonianEOMs()
  print sc
  
  L = m * y_dot**2 / 2 - m*g*y
  sc = SymComputer(L,[y],[y_dot])
  sc.getEulerLegrangeEOMs()
  sc.getHamiltonianEOMs()
  print sc

  L = m * x_dot**2 / 2 + m * y_dot**2 / 2 + m * z_dot**2 / 2 + x_dot*Ax + y_dot*Ay + z_dot*Az - theta
  sc = SymComputer(L,[x,y,z],[x_dot,y_dot,z_dot])
  sc.getEulerLegrangeEOMs()
  sc.getHamiltonianEOMs()
  print sc
  
  L = m*l**2*theta_dot**2/2 + m*g*l*(1-sympy.cos(theta))
  sc = SymComputer(L,[theta],[theta_dot])
  sc.getEulerLegrangeEOMs()
  sc.getHamiltonianEOMs()
  print sc

  L = (M+m+I/R)*x_dot**2/2 + m*l**2*theta_dot**2/2 + m*l*x_dot*theta_dot*sympy.cos(theta) - m*g*l*(1-sympy.cos(theta)) + tau*(theta-x/R)
  sc = SymComputer(L,[x,theta],[x_dot,theta_dot])
  sc.getEulerLegrangeEOMs()
  sc.getHamiltonianEOMs()
  print sc

  L = (M+m)*x_dot**2/2 + I*thetaB_dot**2/2 + m*l**2*theta_dot**2/2 + m*l*x_dot*theta_dot*sympy.cos(theta) - m*g*l*(1-sympy.cos(theta)) + tau*(theta-thetaB) + lam*(R*thetaB+x)
  sc = SymComputer(L,[x,theta,thetaB],[x_dot,theta_dot,thetaB_dot])
  sc.getEulerLegrangeEOMs()
  sc.getHamiltonianEOMs()
  print sc
