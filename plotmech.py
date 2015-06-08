#!/usr/bin/env python

import scipy
from sympy import Symbol
import solvemech
from matplotlib import pyplot as mpl

def plotmech(times,results,coordinates,velmoms,fn):
  if len(coordinates) != len(velmoms):
    raise Exception("coordinates and velmoms must have the same length")

  allvars = coordinates+velmoms
  nvars = len(allvars)
  if nvars != results.shape[1]:
    raise Exception("The length of coordinates + velmoms must equal the number of result varialbes")
  if len(times) != results.shape[0]:
    raise Exception("The length of times must equal the result length")

  mpl.clf()
  for ivar in range(nvars):
    mpl.subplot(nvars, 1, 1+ivar)
    mpl.plot(times,results[:,ivar], 'b-')
    mpl.xlabel("t")
    mpl.ylabel(str(allvars[ivar]))
  mpl.savefig(fn)

if __name__ == "__main__":

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
    M:0.,
    g:10.,
  }

  times = scipy.linspace(0,1,20)
  initialValsL = [0.,10.0]
  initialValsH = [0.,10.0*constValsDict[m]]
  

  L = m*x_dot**2/2 - m*g*x

  sm = solvemech.SolveMech(L,[x],[x_dot])
  timeSeriesH = sm.solveHamiltonian(times,initialValsH,constValsDict)
  plotmech(times,timeSeriesH,["x"],["$p_x$"],"mgH.png")
  timeSeriesL = sm.solveEulerLegrange(times,initialValsL,constValsDict)
  plotmech(times,timeSeriesL,[x],[x_dot],"mgL.png")
