#!/bin/python3
################################################################################
#
#  Used to read input file and make a PW input file
#
#
#

import os
import sys
from bpstandard import oFileData as oFileData
from bpstandard import oStrings as oStrings
from pwIn import pwIn as pwIn

class prepPwIn:

  def __init__(self, fileName=None):
    self.loadDefaults()
    if(fileName is not None):
      self.cfName = fileName
      self.loadFile(self.cfName)


  def loadDefaults(self):
    # directories
    self.outdir = None
    self.ppdir = None


    # template: from file
    self.pwtemplate = None

    # template: build from input
    self.atomList = []
    self.massList = []
    self.ppList = []
    self.structure = None
    self.inputAlat = None
    self.perturb = None

  # template: from file or build from input
    self.cell_xx = 1.0
    self.cell_xy = 0.0
    self.cell_xz = 0.0
    self.cell_yx = 0.0
    self.cell_yy = 1.0
    self.cell_yz = 0.0
    self.cell_zx = 0.0
    self.cell_zy = 0.0
    self.cell_zz = 1.0
    self.copies = 1
    self.nosmear = False
    self.degauss = None
    self.ecutwfc = 50
    self.ecutrho = None
    self.kpoints = "5 5 5   1 1 1"
    self.nspin = None
    self.mixing_mode = None
    self.heat = None
    self.magstructure = "FM"

  def make(self):
    self.pwIn = pwIn()

    # Template: From File
    if self.pwtemplate is not None:
      self.pwIn.loadFile(self.pwtemplate)
      self.pwIn.increaseSize(self.copies)
    # Template: Build
    else:
      self.pwIn.makeTemplate()
      self.pwIn.changeAtomSpecies(self.atomList,self.massList,self.ppList)
      self.pwIn.setMixingLocalTF()
      self.pwIn.makeStructure(self.structure,self.copies,self.copies,self.copies,self.perturb)
      self.pwIn.changeAlat(self.copies * self.inputAlat,True)
      self.nat = self.pwIn.nat
      self.pwIn.changeCell(self.cell_xx,self.cell_xy,self.cell_xz,self.cell_yx,self.cell_yy,self.cell_yz,self.cell_zx,self.cell_zy,self.cell_zz)
      self.pwIn.changeEcutwfc(self.ecutwfc)
      if(self.ecutrho is None):
        self.ecutrho = 4.0 * self.ecutwfc
      self.pwIn.changeEcutrho(self.ecutrho)
      self.pwIn.changeKpoints(self.kpoints)
      if(not self.nosmear):
        self.pwIn.changeDegauss(self.degauss)  # if it exists

    # Template: Both
    self.pwIn.changePrefix("template")
    if(self.outdir is not None):
      if(self.outdir.upper() == "ENV"):
        pwscratch = str(os.environ['PWSCRATCH'])
        if(pwscratch != ""):
          self.pwIn.changeOutdir(pwscratch)
      else:
        self.pwIn.changeOutdir(self.outdir)
    if(self.ppdir is not None):
      if(self.ppdir.upper() == "ENV"):
        pwpp = str(os.environ['PWPP'])
        if(pwpp != ""):
          self.pwIn.changePPdir(pwpp)
      else:
        self.pwIn.changePPdir(self.ppdir)
    # Prefix etc
    if(self.nosmear):
      self.pwIn.removeSmearing()
    # Magnetism
    self.pwIn.setMagStructure(self.magstructure)
    if(self.nspin==0):
      self.pwIn.setMagnetic0()
    if(self.nspin==2):
      self.pwIn.setMagnetic2()
    # Mixing Mode
    if(self.mixing_mode=="TF"):
      self.pwIn.setMixingTF()
    if(self.mixing_mode=="local-TF"):
      self.pwIn.setMixingLocalTF()

    # Heat coords
    #if(self.heat is not None):
    #  self.pwIn.heat(self.heat)

  def getPW(self):
    return self.pwIn


  # Static methods

  @staticmethod
  def checkKeyword(lineIn, keyword, verbose=False):
    result = None
    if(lineIn != ""):
      lineInArr = lineIn.split("#")
      lineIn = lineInArr[0]
      if(lineIn != ""):
        lineInUC = lineIn.upper()
        keywordUC = keyword.upper()
        keywordLen = len(keywordUC)
        if(lineInUC[0:keywordLen]==keywordUC):
          result = oStrings.removeDouble(lineIn," ")
          if(verbose):
            print(result)
    return result




################################################################################

  def loadFile(self, fileName):
    print("name "+fileName)
    self.cFile = oFileData()
    self.cFile.loadFile(fileName)

    for i in range(0,self.cFile.lineCount):
      # Directories
      #######################
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$OUTDIR", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.outdir = fileRowArr[1]
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$PPDIR", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.ppdir = fileRowArr[1]

      # Template: from file
      #######################
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$PWTEMPLATE", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.pwtemplate = fileRowArr[1]

      # Template: build from input
      #######################
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ATOMLIST", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        for i in range(1,len(fileRowArr)):
          self.atomList.append(fileRowArr[i])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$MASSLIST", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        for i in range(1,len(fileRowArr)):
          self.massList.append(fileRowArr[i])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$PPLIST", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        for i in range(1,len(fileRowArr)):
          self.ppList.append(fileRowArr[i])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$STRUCTURE", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.structure = fileRowArr[1]
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$MAGSTRUCTURE", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.magstructure = fileRowArr[1]
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ALAT", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.inputAlat = float(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$DEGAUSS", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.degauss = float(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ECUTWFC", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.ecutwfc = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ECUTRHO", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.ecutrho = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$KPOINTS", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.kpoints = str(fileRowArr[1])+" "+str(fileRowArr[2])+" "
        self.kpoints = self.kpoints + str(fileRowArr[3])+" "+str(fileRowArr[4])+" "
        self.kpoints = self.kpoints + str(fileRowArr[5])+" "+str(fileRowArr[6])


      # template: from file or build from input
      #######################
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$CELL", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.cell_xx = float(fileRowArr[1])
        self.cell_xy = float(fileRowArr[2])
        self.cell_xz = float(fileRowArr[3])
        self.cell_yx = float(fileRowArr[4])
        self.cell_yy = float(fileRowArr[5])
        self.cell_yz = float(fileRowArr[6])
        self.cell_zx = float(fileRowArr[7])
        self.cell_zy = float(fileRowArr[8])
        self.cell_zz = float(fileRowArr[9])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$COPIES", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.copies = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$HEAT", False)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.heat = float(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$NOSMEAR", False)
      if(keywordResult is not None):
        self.nosmear = True
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$NSPIN", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.nspin = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$MIXING_MODE", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.mixing_mode = fileRowArr[1]






################################################################################
