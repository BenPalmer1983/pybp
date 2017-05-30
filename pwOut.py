#!/bin/python3
################################################################################

import math
import os

from bpstandard import oFileData as oFileData
from bpstandard import oStrings as oStrings



###########################
# pwscf output file
###########################
class pwOut:

  def __init__(self, fileName=None):
    self.successful = False
    if(fileName is not None):
      self.loadDataFile(fileName)

  def loadDataFile(self, fileName):
    if os.path.isfile(fileName):
      self.pwFile = oFileData()               # make new oFileData object
      self.pwFile.loadFile(fileName)
      self.successful = self.checkSuccessful()
      if(self.successful):
        self.readData()

  def readData(self):
    self.atomCount = self.getAtomCount()
    self.totalEnergy = self.getEnergy()
    self.totalForce = self.getForce()
    self.stress = self.getStress()
    self.stressAvg = self.getStressAvg()
    # Per atom
    self.energyPerAtom = self.totalEnergy / self.atomCount
    self.forcePerAtom = self.totalForce / self.atomCount
    # for plotting
    self.energyPerAtomPlot = self.energyPerAtom - math.floor(float(1000 * self.energyPerAtom)) / 1000
    self.forcePerAtomPlot = self.forcePerAtom - math.floor(float(1000 * self.forcePerAtom)) / 1000
    self.stressPlot = self.stressAvg - math.floor(float(1000 * self.stressAvg)) / 1000

  def checkSuccessful(self):
    result = False
    for i in range(0,self.pwFile.lineCount):
      if '   JOB DONE.' in self.pwFile.fileData[i]:
        result = True
    return result

  def getEnergy(self):
    energy = 0.0
    for i in range(0,self.pwFile.lineCount):
      if '!    total energy' in self.pwFile.fileData[i]:
        tString = self.pwFile.fileData[i].upper()
        tString = tString.split("=")
        tString = tString[1]
        tString = tString.split("RY")
        energy = float(tString[0])
    return energy

  def getEnergyPerAtom(self):
    energy = 0.0
    for i in range(0,self.pwFile.lineCount):
      if 'number of atoms/cell      =' in self.pwFile.fileData[i]:
        tString = self.pwFile.fileData[i].upper()
        tString = tString.split("=")
        nAtoms = int(tString[1])
      if '!    total energy' in self.pwFile.fileData[i]:
        tString = self.pwFile.fileData[i].upper()
        tString = tString.split("=")
        tString = tString[1]
        tString = tString.split("RY")
        energy = float(tString[0])
    energyPerAtom = energy / nAtoms
    return energyPerAtom

  def getForce(self):
    force = 0.0
    for i in range(0,self.pwFile.lineCount):
      if '     Total force =' in self.pwFile.fileData[i]:
        tString = self.pwFile.fileData[i][20:34]
        force = float(tString)
    return force

  def getStress(self):
    stress = [None]*9
    for i in range(0,self.pwFile.lineCount):
      if 'total   stress ' in self.pwFile.fileData[i]:
        i = i + 1
        fileRow = oStrings.removeDouble(self.pwFile.fileData[i]," ")
        fileRow = oStrings.trimEnds(fileRow)
        stressArr = fileRow.split(" ")
        stress[0] = float(stressArr[0])
        stress[1] = float(stressArr[1])
        stress[2] = float(stressArr[2])
        i = i + 1
        fileRow = oStrings.removeDouble(self.pwFile.fileData[i]," ")
        fileRow = oStrings.trimEnds(fileRow)
        stressArr = fileRow.split(" ")
        stress[3] = float(stressArr[0])
        stress[4] = float(stressArr[1])
        stress[5] = float(stressArr[2])
        i = i + 1
        fileRow = oStrings.removeDouble(self.pwFile.fileData[i]," ")
        fileRow = oStrings.trimEnds(fileRow)
        stressArr = fileRow.split(" ")
        stress[6] = float(stressArr[0])
        stress[7] = float(stressArr[1])
        stress[8] = float(stressArr[2])
    return stress

  def getStressAvg(self):
    self.stressAvg = 0.0
    for i in range(0,8):
      self.stressAvg = self.stressAvg + abs(self.stress[i])/9.0
    return self.stressAvg

  def getAlat(self):
    for i in range(0,self.pwFile.lineCount):
      if 'lattice parameter (alat)  =' in self.pwFile.fileData[i]:
        tString = self.pwFile.fileData[i].upper()
        tString = tString.split("=")
        tString = tString[1]
        tString = tString.split("A.U.")
        aLat = float(tString[0])
    return aLat

  def getAlat_Relaxed(self):
    inRelaxed = False
    for i in range(0,self.pwFile.lineCount):
      if 'A final scf calculation' in self.pwFile.fileData[i]:
        inRelaxed = True
      if(inRelaxed):
        if 'lattice parameter (alat)  =' in self.pwFile.fileData[i]:
          tString = self.pwFile.fileData[i].upper()
          tString = tString.split("=")
          tString = tString[1]
          tString = tString.split("A.U.")
          aLat = float(tString[0])
    return aLat

  def getAlat_RelaxedNormalised(self):
    inRelaxed = False
    for i in range(0,self.pwFile.lineCount):
      if 'A final scf calculation' in self.pwFile.fileData[i]:
        inRelaxed = True
      if(inRelaxed):
        if 'lattice parameter (alat)  =' in self.pwFile.fileData[i]:
          tString = self.pwFile.fileData[i].upper()
          tString = tString.split("=")
          tString = tString[1]
          tString = tString.split("A.U.")
          aLat = float(tString[0])
        if 'crystal axes:' in self.pwFile.fileData[i]:
          i = i + 1
          dataRow = oStrings.removeDouble(self.pwFile.fileData[i]," ")
          dataRow = oStrings.trimEnds(dataRow)
          dataArr = dataRow.split(" ")
          aLat = aLat * float(dataArr[3])
    return aLat

  def getCell_Relaxed(self):
    cell = [[0 for x in range(3)] for y in range(3)]
    inRelaxed = False
    for i in range(0,self.pwFile.lineCount):
      if 'A final scf calculation' in self.pwFile.fileData[i]:
        inRelaxed = True
      if(inRelaxed):
        if 'crystal axes:' in self.pwFile.fileData[i]:
          i = i + 1
          dataRow = oStrings.removeDouble(self.pwFile.fileData[i]," ")
          dataRow = oStrings.trimEnds(dataRow)
          dataArr = dataRow.split(" ")
          cell[0][0] = float(dataArr[3])
          cell[0][1] = float(dataArr[4])
          cell[0][2] = float(dataArr[5])
          i = i + 1
          dataRow = oStrings.removeDouble(self.pwFile.fileData[i]," ")
          dataRow = oStrings.trimEnds(dataRow)
          dataArr = dataRow.split(" ")
          cell[1][0] = float(dataArr[3])
          cell[1][1] = float(dataArr[4])
          cell[1][2] = float(dataArr[5])
          i = i + 1
          dataRow = oStrings.removeDouble(self.pwFile.fileData[i]," ")
          dataRow = oStrings.trimEnds(dataRow)
          dataArr = dataRow.split(" ")
          cell[2][0] = float(dataArr[3])
          cell[2][1] = float(dataArr[4])
          cell[2][2] = float(dataArr[5])
    return cell

  def getVol(self):
    vol = 0
    for i in range(0,self.pwFile.lineCount):
      if 'unit-cell volume          =' in self.pwFile.fileData[i]:
        tString = self.pwFile.fileData[i].upper()
        tString = tString.split("=")
        tString = tString[1]
        tString = tString.split("(A.U.)^3")
        vol = float(tString[0])
    return vol

  def getVolPerAtom(self):
    vol = 0
    for i in range(0,self.pwFile.lineCount):
      if 'unit-cell volume          =' in self.pwFile.fileData[i]:
        tString = self.pwFile.fileData[i].upper()
        tString = tString.split("=")
        tString = tString[1]
        tString = tString.split("(A.U.)^3")
        vol = float(tString[0])
      if 'number of atoms/cell      =' in self.pwFile.fileData[i]:
        tString = self.pwFile.fileData[i].upper()
        tString = tString.split("=")
        nAtoms = int(tString[1])
    volPerAtom = vol / nAtoms
    return volPerAtom

  def getAtomCount(self):
    for i in range(0,self.pwFile.lineCount):
      if 'number of atoms/cell      =' in self.pwFile.fileData[i]:
        tString = self.pwFile.fileData[i].upper()
        tString = tString.split("=")
        nAtoms = int(tString[1])
    return nAtoms




###########################
# Batch
###########################
class pwOutBatch:
  def __init__(self):
    self.count = 0

  def loadFile(self, fileName):
    self.count = self.count + 1
    pwData = pwOut(outFile)
    self.files.append(pwData)


###########################
# Batch
###########################
class pwCompare:
  def __init__(self, pwOutA, pwOutB):
    self.pwOutA = pwOutA
    self.pwOutB = pwOutB
    self.pwOutA.readData()
    self.pwOutB.readData()

  def compareEnergy(self):
    eA = self.pwOutA.totalEnergy
    eB = self.pwOutB.totalEnergy
    eDiff = abs(eB-eA)
    return eDiff

  def compareEnergyPA(self):
    eA = self.pwOutA.energyPerAtom
    eB = self.pwOutB.energyPerAtom
    eDiff = abs(eB-eA)
    return eDiff

  def compareForce(self):
    eA = self.pwOutA.totalForce
    eB = self.pwOutB.totalForce
    eDiff = abs(eB-eA)
    return eDiff

  def compareForcePA(self):
    eA = self.pwOutA.forcePerAtom
    eB = self.pwOutB.forcePerAtom
    eDiff = abs(eB-eA)
    return eDiff

  def compareStress(self):
    sDiff = 0.0
    for i in range(0,9):
      sDiff = sDiff + abs(self.pwOutB.stress[i]-self.pwOutA.stress[i])
    sDiff = sDiff / 9.0
    return sDiff



################################################################################
