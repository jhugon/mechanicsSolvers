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

  def __init__(self,lagrangian,coordinates,t=Symbol('t'),velocity_suffix="_dot"):
    for coord in coordinates:
      name = str(coord)
      if "_dot" in name:
        raise Exception("Coordinate '{0}' contains the velocity_suffix '{1}'".format(name,velocity_suffix))
      if name == str(t):
        raise Exception("Coordinate '{0}' and time variable '{1}' have the same name".format(name,str(t)))

    self.coordinates = coordinates
    self.lagrangian = lagrangian
    self.t = t
    self.velocity_suffix = velocity_suffix

    self.velocities = None
    self.tdm = None

    self.momenta = None
    self.momentaDefs = None
    self.hamiltonian = None

    self.eomsE = None
    self.eomsL = None

    self.setup()

  def setup(self):
    self.velocities = [str(x)+velocity_suffix]
    self.tdm = symmech.generateTimeDerivMapFromCoords(self.coordinates)
    

if __name__ == "__main__":

  def printSteps(L,coords,vels,t,tdm):
    print "########################################"
    print "########################################"
    print "########################################"
    momenta, momentaDefs, H = makeMomentaHamiltonian(L,coords,vels,tdm)
    elEOM =  eulerLegrangeEOM(L,coords,vels,t,tdm)
    hEOM =  hamiltonianEOM(H,coords,momenta,tdm)
    print "L = "+str(L)
    print "H = "+str(simplify(H))
    print "momenta:"
    for mdef in momentaDefs:
      print "  "+str(mdef)
    print "Euler-Legrange EOM:"
    for eom in elEOM:
      print "  "+str(simplify(eom))
    print "Hamiltonian EOM:"
    for eom in hEOM:
      print "  "+str(simplify(eom))
    print "########################################"
    print "########################################"
    print "########################################"
    

  g = symbols('g')
  m = symbols('m')
  M = symbols('M')
  l = symbols('l')
  R = symbols('R')
  I = symbols('I')
  K = symbols('K')
  F = symbols('F')
  tau = symbols('tau')

  t = symbols('t')

  x = symbols('x')
  x_dot = symbols('x_dot')
  y = symbols('y')
  y_dot = symbols('y_dot')
  z = symbols('z')
  z_dot = symbols('z_dot')

  theta = symbols('theta')
  theta_dot = symbols('theta_dot')

  Ax = symbols('Ax')
  Ay = symbols('Ay')
  Az = symbols('Az')

  tdm = generateTimeDerivMapFromCoords([x,y,z,Ax,Ay,Az,theta])
  
  #L = m * x_dot**2 / 2 + m * y_dot**2 / 2 + m * z_dot**2 / 2
  #printSteps(L,[x,y,z],[x_dot,y_dot,z_dot],t,tdm)
  
  #L = m * x_dot**2 / 2 + F*x
  #printSteps(L,[x],[x_dot],t,tdm)
  
  #
  #L = m * y_dot**2 / 2 - m*g*y
  #printSteps(L,[y],[y_dot],t,tdm)

  #L = m * x_dot**2 / 2 + m * y_dot**2 / 2 + m * z_dot**2 / 2 + x_dot*Ax + y_dot*Ay + z_dot*Az - theta
  #printSteps(L,[x,y,z],[x_dot,y_dot,z_dot],t,tdm)
  #
  #L = m*l**2*theta_dot**2/2 + m*g*l*sympy.cos(theta)
  #printSteps(L,[theta],[theta_dot],t,tdm)

  #L = (M+m)*x_dot**2/2 + m*x_dot*l*theta_dot*sympy.cos(theta) + m*l**2*theta_dot**2/2 + m*g*l*sympy.cos(theta)
  #printSteps(L,[x,theta],[x_dot,theta_dot],t,tdm)

  L = M*x_dot**2/2 + I*theta_dot**2/2 + K*x_dot*theta_dot*sympy.cos(theta)+g*K*sympy.cos(theta) + tau*theta
  printSteps(L,[x,theta],[x_dot,theta_dot],t,tdm)


