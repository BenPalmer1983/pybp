#!/bin/python3
################################################################################

import sys
import os
from randnum import RandomLCG as RandomLCG
import numpy

class eosO:
  def __init__(self, B0, B0P, V0, E0):
    self.B0 = B0
    self.B0P = B0P
    self.V0 = V0
    self.E0 = E0
    self.b0Gpa()

  def b0Gpa(self):
    self.B0_GPa = self.B0 * 1.47105e4

class eos:
  def __init__(self):
    self.eosO = eosO(0,0,0,0)
    self.rss = 0

  def loadData(self, volume, energy):
    self.volume = volume
    self.energy = energy

  def fitEoS(self):
    # 2nd order polynomial fit
    self.pFit = numpy.polyfit(self.volume, self.energy, 2)
    # Starting points
    self.eosO.V0 = (-1 * self.pFit[1]) / (2 * self.pFit[0])
    self.eosO.E0 = (self.pFit[0] * self.eosO.V0 * self.eosO.V0) + (self.pFit[1] * self.eosO.V0) + self.pFit[2]
    self.eosO.B0 = 2.0e0 * self.pFit[0] * self.eosO.V0
    self.eosO.B0P = 2.0e0
    self.eosO.b0Gpa()
    # Calculate error between fit and data
    optRss = self.EoS_rss(self.eosO, self.volume, self.energy)
    # Randomise and pick best options
    rg=RandomLCG()      # Make new rand
    trialEoS = self.eosO
    trialRss = self.EoS_rss(trialEoS, self.volume, self.energy)
    for n in range (0,1000):
      if(n<=100):
        self.perturbB0P(trialEoS, rg)
      else:
        self.perturb(trialEoS, rg)
      trialRss = self.EoS_rss(trialEoS, self.volume, self.energy)
      if(trialRss<optRss):
        optRss = trialRss
        self.eosO = trialEoS
      trialEoS = self.eosO
    self.eosO.rss = self.EoS_rss(self.eosO, self.volume, self.energy)

  def EoS_rss(self, eosO, volumeIn, energyIn):
    rss = 0.0e0
    i = 0
    for V in volumeIn:
      E = self.calc(eosO, V)
      rss = rss + (E - energyIn[i])**2
      i = i + 1
    return rss

  def perturb(self, eosO, rng):
    # Perturb each value by a small amount (fraction of 0.1%)
    rVal = rng.rng()
    eosO.B0 = eosO.B0 * (1.0 + 0.001 * (rVal - 0.5))
    rVal = rng.rng()
    eosO.B0P = eosO.B0P * (1.0 + 0.001 * (rVal - 0.5))
    rVal = rng.rng()
    eosO.E0 = eosO.E0 * (1.0 + 0.001 * (rVal - 0.5))
    rVal = rng.rng()
    eosO.V0 = eosO.V0 * (1.0 + 0.001 * (rVal - 0.5))

  def perturbB0P(self, eosO, rng):
    # Random value between 2 and 9
    rVal = rng.rng()
    eosO.B0P = 2.0 + 7.0 * rVal

  def calc(self, eosO, V):
    eta = (V/eosO.V0)**(1/3.0)
    E = eosO.E0 + (9/16.0) * (eosO.B0 * eosO.V0) * ((eta*eta - 1)*(eta*eta - 1)) * (6.0 + eosO.B0P * (eta * eta - 1) - 4 * eta * eta )
    return E

  def display(self, frontPadding=""):
    self.eosO.b0Gpa()
    print(frontPadding,"B0   ",self.eosO.B0,"   B0 (GPa)  ",self.eosO.B0_GPa)
    print(frontPadding,"B0P  ",self.eosO.B0P)
    print(frontPadding,"E0   ",self.eosO.E0)
    print(frontPadding,"V0   ",self.eosO.V0)
    print(frontPadding,"rss  ",self.eosO.rss)











################################################################################
