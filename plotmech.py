#!/usr/bin/env python3

import scipy
import solvemech
from matplotlib import pyplot as mpl


def plotmech(times, results, coordinates, velmoms, fn):
    if len(coordinates) != len(velmoms):
        raise Exception("coordinates and velmoms must have the same length")

    allvars = coordinates + velmoms
    nvars = len(allvars)
    if nvars != results.shape[1]:
        raise Exception(
            "The length of coordinates + velmoms must equal the number of result variables"
        )
    if len(times) != results.shape[0]:
        raise Exception("The length of times must equal the result length")

    mpl.clf()
    for ivar in range(nvars):
        mpl.subplot(nvars, 1, 1 + ivar)
        mpl.plot(times, results[:, ivar], "b-")
        mpl.xlabel("t")
        mpl.ylabel(str(allvars[ivar]))
    mpl.savefig(fn)


if __name__ == "__main__":

    import sympy
    from sympy import Symbol

    t = Symbol("t")
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    theta = Symbol("theta")
    thetaB = Symbol("thetaB")
    x_dot = Symbol("x_dot")
    y_dot = Symbol("y_dot")
    z_dot = Symbol("z_dot")
    theta_dot = Symbol("theta_dot")
    thetaB_dot = Symbol("thetaB_dot")
    m = Symbol("m")
    M = Symbol("M")
    g = Symbol("g")
    I = Symbol("I")
    R = Symbol("R")
    l = Symbol("l")
    tau = Symbol("tau")
    lam = Symbol("lam")
    hev = Symbol("hev")

    constValsDict = {
        m: 20.0,
        M: 20.0,
        g: 9.81,
        I: 1.0,
        tau: lambda d: 10 if d["t"] < 1 else -20,
        R: 1.0,
        l: 1.0,
        lam: 1e-3,
        hev: lambda d: 0 if d["t"] < 1 else 1,
    }

    times = scipy.linspace(0, 3, 5000)
    initialValsL = [0.0, 10.0]
    initialValsH = [0.0, 10.0 * constValsDict[m]]

    L = m * x_dot ** 2 / 2 - m * g * x * hev

    sm = solvemech.SolveMech(L, [x], [x_dot])
    timeSeriesH = sm.solveHamiltonian(times, initialValsH, constValsDict)
    plotmech(times, timeSeriesH, ["x"], ["$p_x$"], "mgH.png")
    timeSeriesL = sm.solveEulerLegrange(times, initialValsL, constValsDict)
    plotmech(times, timeSeriesL, [x], [x_dot], "mgL.png")

    L2 = (
        (M + m + I / R) * x_dot ** 2 / 2
        + m * l ** 2 * theta_dot ** 2 / 2
        + m * l * x_dot * theta_dot * sympy.cos(theta)
        - m * g * l * (1 - sympy.cos(theta))
        + tau * (theta - x / R)
    )
    sm2 = solvemech.SolveMech(L2, [x, theta], [x_dot, theta_dot])
    timeSeriesL2 = sm2.solveEulerLegrange(times, [0.0, 0.0, 0.0, 0.0], constValsDict)
    plotmech(times, timeSeriesL2, [x, theta], [x_dot, theta_dot], "L2.png")

    L3 = (
        (M + m) * x_dot ** 2 / 2
        + I * thetaB_dot ** 2 / 2
        + m * l ** 2 * theta_dot ** 2 / 2
        + m * l * x_dot * theta_dot * sympy.cos(theta)
        - m * g * l * (1 - sympy.cos(theta))
        + tau * (theta - thetaB)
        + lam * (R * thetaB + x)
    )
    sm3 = solvemech.SolveMech(L3, [x, theta, thetaB], [x_dot, theta_dot, thetaB_dot])
    timeSeriesL3 = sm3.solveEulerLegrange(
        times, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], constValsDict
    )
    plotmech(
        times,
        timeSeriesL3,
        [x, theta, thetaB],
        [x_dot, theta_dot, thetaB_dot],
        "L3.png",
    )
