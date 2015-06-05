#!/usr/bin/env python

import symmech
import sympy
from sympy import Symbol

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

  L = m*x_dot**2/2

  tdm = symmech.generateTimeDerivMapFromCoords(x,y,z)

  #eomsL = symmech.eulerLegrangeEOM(L,[x],[x_dot],t,tdm)
  
  momenta, momentaDefs, H = symmech.makeMomentaHamiltonian(L,[x],[x_dot],tdm)
  print momenta
  print momentaDefs
  print H

  eomsH = symmech.hamiltonianEOM(H,[x],momenta,tdm)
  print eomsH
