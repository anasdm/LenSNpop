import glob
import numpy as np
import matplotlib.pyplot as plt
import pickle

class ColorImage:
    def __init__(self):
        self.rnorm = None
        self.gnorm = None
        self.I = None
        self.M = None
        self.nonlin = 5.0
        self.m = 0.5
        self.bMinusr = 0.8
        self.bMinusg = 0.4

    def clip(self, arr, nsig=3.5):
        a = arr.flatten()
        a.sort()
        a = a[int(a.size * 0.05):int(a.size * 0.8)]
        m, s, l = a.mean(), a.std(), a.size
        while True:
            a = a[np.abs(a - m) < s * nsig]
            if a.size == l:
                return m, s
            m, s, l = a.mean(), a.std(), a.size

    def createModel(self, b, g, r):
        bMinusr = self.bMinusr
        bMinusg = self.bMinusg
        b0 = b.copy()
        g0 = g.copy()
        r0 = r.copy()
        
        w = int(r.shape[0] / 2 - 5)
        rb = r0 / b0
        gb = g0 / b0
        rnorm = np.median(rb[w:-w, w:-w])
        gnorm = np.median(gb[w:-w, w:-w])
        r0 /= rnorm
        g0 /= gnorm
        r0 *= 10**(0.4 * bMinusr)
        g0 *= 10**(0.4 * bMinusg)

        r0 /= 620.0
        g0 /= 540.0
        b0 /= 460.0

        I = (r0 + g0 + b0) / 3.0
        self.I = I
        self.rnorm = rnorm
        self.gnorm = gnorm
        return self.colorize(b, g, r)

    def colorize(self, b, g, r, newI=False):
        bMinusr = self.bMinusr
        bMinusg = self.bMinusg
        rnorm = self.rnorm
        gnorm = self.gnorm
        m = self.m
        nonlin = self.nonlin
        I = self.I.copy()

        b = b.copy()
        g = g.copy()
        r = r.copy()

        w = int(r.shape[0] / 2 - 5)
        r /= rnorm
        g /= gnorm
        r *= 10**(0.4 * bMinusr)
        g *= 10**(0.4 * bMinusg)

        r /= 620.0
        g /= 540.0
        b /= 460.0

        sdev = self.clip(I)[1]
        m *= sdev
        if self.M is None:
            M = I[w:-w, w:-w].max()
        else:
            M = self.M
        nonlin *= sdev

        if newI:
            I = (b + g + r) / 3.0
        f = np.arcsinh((I - m) / nonlin) / np.arcsinh((M - m) / nonlin)
        f[I < m] = 0.0
        f[I > M] = 1.0
        R = r * f / I
        G = g * f / I
        B = b * f / I

        R[I <= 0] = 0.0
        G[I <= 0] = 0.0
        B[I <= 0] = 0.0

        R[R <= 0] = 0.0
        G[G <= 0] = 0.0
        B[B <= 0] = 0.0

        R[R > 1] = 1.0
        G[G > 1] = 1.0
        B[B > 1] = 1.0

        white = True
        if white:
            cond = (f == 1)
            R[cond] = 1.0
            G[cond] = 1.0
            B[cond] = 1.0

        arr = np.empty((R.shape[0], R.shape[1], 3))
        arr[:, :, 0] = R
        arr[:, :, 1] = G
        arr[:, :, 2] = B

        return arr
