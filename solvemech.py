#!/usr/bin/env python3

from symmech import SymComputer
import sympy
from sympy import Symbol
import scipy
import scipy.integrate
import types

class SolveMech(object):

  def __init__(self,lagrangian,coordinates,velocities,t=Symbol("t"),derivative_suffix="_dot",velocity_suffix="_v"):
    """
    Takes a lagrangian, made from sympy variables. Coordinates and velocities
    are lists of sympy variables used in the lagrangian. Make sure that coorinates[i]
    matches velocities[i].
    """

    self.symComputer = SymComputer(lagrangian,coordinates,velocities,t=t)
    if derivative_suffix in velocity_suffix or velocity_suffix in derivative_suffix:
        raise Exception("derivative_suffix and velocity_suffix can't be substrings of eachother")
    self.derivative_suffix = derivative_suffix
    self.velocity_suffix = velocity_suffix
    self.diffeq = None
    self.constValsList = None
    self.allVarList = None

  def __str__(self):
    result = ""
    result += "****SolveMech**** "
    result += str(self.symComputer)
    return result

  def __call__(self,yarr,t):
    #print "call: ",yarr,t
    constVals = []
    for p in self.constValsList:
      try:
        constVals.append(float(p))
      except TypeError:
        constVals.append(self.callConstFunc(p,yarr,t))
    allvars = [t]+list(yarr)+constVals
    result = scipy.zeros(len(yarr))
    for i in range(len(yarr)):
        result[i] = self.diffeq[i](*allvars)
    return result

  def solveEulerLegrange(self,times,initialVals,constantValueDict):
    """
    Requires a list of time values to solve for, a list of initial values for
    the coordinates concatenated with one for the velocities, and a dictionary of
    values to assign to the constant variables (where the keys are sympy
    variables).

    returns a list of lists for each time point, where for each time
    point is a list of values for each coord and each velocity in the
    same order as in the initial values.

    Order of initial values and results:

    [coord1,coord2,...,vel_coord1,vel_coord2,...]

    Make sure initial values are in the same order as the coords and velocities
    fed to the contructor!
    """
    eoms = self.symComputer.getEulerLegrangeEOMs()
    print(eoms)
    eomDict = self.rearangeForDerivs(eoms,secondOrder=True)
    print(self.symComputer.tdm)
    print(eomDict)

    varList = []
    for coord in self.symComputer.coordinates:
      varList.append(coord)
    for vel in self.symComputer.velocities:
      varList.append(vel)

    constList = list(constantValueDict.keys())
    self.constValsList = [constantValueDict[x] for x in constList]
    allVarList = [self.symComputer.t]+varList+constList
    self.allVarList = allVarList

    ## create functions to be called with list [t]+[coord1,coord2,...]+[vel1,vel2,...]+[const1,const2,...]
    exprList = []
    for var in varList:
      var_dot = self.symComputer.tdm[var]
      if var_dot in eomDict:
        exprList.append(eomDict[var_dot])
      else:
        exprList.append(var_dot)
    print("vars:", varList)
    print("diffexps:", exprList)
    print("allVarList: ",allVarList)
    self.diffeq = []
    for expr in exprList:
      self.diffeq.append(sympy.lambdify(allVarList,expr))

    #print "about to call odeint!"
    #print "with ", self, initialVals, times
    result = scipy.integrate.odeint(self,initialVals,times)
    #print result

    self.diffeq = None
    self.constValsList = None
    self.allVarList = None

    return result


  def solveHamiltonian(self,times,initialVals,constantValueDict):
    """
    Requires a list of time values to solve for, a list of initial values for
    the coordinates concatenated with one for the momenta, and a dictionary of
    values to assign to the constant variables (where the keys are sympy
    variables).

    returns a list of lists for each time point, where for each time
    point is a list of values for each coord and each momenta in the
    same order as in the initial values.

    Order of initial values and results:

    [coord1,coord2,...,momentum_coord1,momentum_coord2,...]

    Make sure initial values are in the same order as the coords and velocities
    fed to the contructor!
    """

    eoms = self.symComputer.getHamiltonianEOMs()
    eomDict = self.rearangeForDerivs(eoms)

    varList = []
    for coord in self.symComputer.coordinates:
      varList.append(coord)
    for mom in self.symComputer.momenta:
      varList.append(mom)

    constList = list(constantValueDict.keys())
    self.constValsList = [constantValueDict[x] for x in constList]
    allVarList = [self.symComputer.t]+varList+constList
    self.allVarList = allVarList

    ## create functions to be called with list [t]+[coord1,coord2,...]+[mom1,mom2,...]+[const1,const2,...]
    exprList = []
    for var in varList:
      var_dot = self.symComputer.tdm[var]
      exprList.append(eomDict[var_dot])
    self.diffeq = []
    for expr in exprList:
      self.diffeq.append(sympy.lambdify(allVarList,expr))

    #print "about to call odeint!"
    #print "with ", self, initialVals, times
    result = scipy.integrate.odeint(self,initialVals,times)
    #print result

    self.diffeq = None
    self.constValsList = None
    self.allVarList = None

    return result

  def rearangeForDerivs(self,eoms,secondOrder=False):
    suffix = self.derivative_suffix
    if secondOrder:
        suffix += suffix
    derivSet = set()
    for eom in eoms:
      for x in sympy.preorder_traversal(eom):
        if type(x) == sympy.Symbol and suffix in str(x):
          derivSet.add(x)
    derivList = list(derivSet)
    if len(derivList) == 0:
      raise Exception("No time derivative was found in eom: '{0}', where derivative has suffix of '{1}'".format(eom,suffix))
    solns = sympy.solve(eoms,derivList)
    if len(solns) < len(derivList):
      raise Exception("Too few solutions, {2}, were found for '{1}'".format(eom,derivList,solns))
    if len(solns) > len(derivList):
      raise Exception("Too many solutions, {2}, were found for '{1}'".format(eom,derivList,solns))
    for deriv in derivList:
      try:
        solns[deriv]
      except KeyError:
        raise Exception("Couldn't find solution for {1}, in solutions {2}, in EOMs '{0}'".format(eom,deriv,solns))
    return solns

  def callConstFunc(self,func,yarr,t):
    allVarNameList = [str(i) for i in self.allVarList]
    params = {}
    for name,val in zip(allVarNameList,[t]+list(yarr)+self.constValsList):
      try:
        params[name] = float(val)
      except TypeError:
        pass
    return func(params)

if __name__ == "__main__":
  from matplotlib import pyplot as mpl

  t = Symbol("t")
  x = Symbol("x")
  y = Symbol("y")
  z = Symbol("z")
  x_dot = Symbol("x_dot")
  y_dot = Symbol("y_dot")
  z_dot = Symbol("z_dot")
  m = Symbol("m")
  M = Symbol("M")
  g = Symbol("g")
  constValsDict = {
    m:20.,
    #M:0.,
    g:10.,
    #g: lambda d: 10.-10*d['t'],
  }

  times = scipy.linspace(0,1,20)
  initialValsL = [0.,10.0]
  initialValsH = [0.,10.0*constValsDict[m]]
  

  L = m*x_dot**2/2 - m*g*x

  sm = SolveMech(L,[x],[x_dot])
  timeSeriesH = sm.solveHamiltonian(times,initialValsH,constValsDict)
  timeSeriesL = sm.solveEulerLegrange(times,initialValsL,constValsDict)
  print(sm)

  mpl.subplot(2, 1, 1)
  mpl.plot(times,timeSeriesH[:,0], 'b-')
  mpl.plot(times,timeSeriesL[:,0], 'ro')
  mpl.ylabel('Position')
  mpl
  mpl.subplot(2, 1, 2)
  mpl.plot(times,timeSeriesH[:,1], 'b-')
  mpl.plot(times,timeSeriesL[:,1], 'ro')
  mpl.xlabel('time (s)')
  mpl.ylabel('Velocity/Momentum')
  mpl
  mpl.show()
