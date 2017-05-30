#!/bin/python3
################################################################################


###########################
# pwscfUnitVector
###########################
class strains:
  @staticmethod
  def identityMatrix():
    iMatrix = [[0 for x in range(3)] for y in range(3)]
    for i in range (0,3):
      for j in range (0,3):
        if(i==j):
          iMatrix[i][j] = 1.0e0
        else:
          iMatrix[i][j] = 0.0e0
    return iMatrix

#################
# Cubic Strains
#################

  @staticmethod
  def orthorhombic(strain):
    sVec = strains.identityMatrix()
    sVec[0][0] = sVec[0][0] + strain
    sVec[1][1] = sVec[1][1] - strain
    sVec[2][2] = sVec[2][2] + (strain * strain) / (1.0 - strain * strain)
    return sVec

  @staticmethod
  def tetragonal(strain):
    sVec = strains.identityMatrix()
    sVec[1][0] = sVec[1][0] + strain/2.0  # e6/2
    sVec[0][1] = sVec[0][1] + strain/2.0  # e6/2
    sVec[2][2] = sVec[2][2] + (strain * strain) / (4.0 - strain * strain)   # e3
    return sVec

class structures:
  @staticmethod
  def primitiveStructure(structure, symbol):
    ####################
    # ISOLATED
    ####################
    if(structure=="ISO"):
      primitive = [[0 for x in range(4)] for y in range(1)]
      primitive[0][0] = str(symbol)
      primitive[0][1] = 0.0
      primitive[0][2] = 0.0
      primitive[0][3] = 0.0
    ####################
    # FCC
    ####################
    if(structure=="FCC"):
      primitive = [[0 for x in range(4)] for y in range(4)]
      primitive[0][0] = str(symbol)
      primitive[0][1] = 0.0
      primitive[0][2] = 0.0
      primitive[0][3] = 0.0
      primitive[1][0] = str(symbol)
      primitive[1][1] = 0.5
      primitive[1][2] = 0.5
      primitive[1][3] = 0.0
      primitive[2][0] = str(symbol)
      primitive[2][1] = 0.5
      primitive[2][2] = 0.0
      primitive[2][3] = 0.5
      primitive[3][0] = str(symbol)
      primitive[3][1] = 0.0
      primitive[3][2] = 0.5
      primitive[3][3] = 0.5
    ####################
    # BCC
    ####################
    if(structure=="BCC"):
      primitive = [[0 for x in range(4)] for y in range(2)]
      primitive[0][0] = str(symbol)
      primitive[0][1] = 0.0
      primitive[0][2] = 0.0
      primitive[0][3] = 0.0
      primitive[1][0] = str(symbol)
      primitive[1][1] = 0.5
      primitive[1][2] = 0.5
      primitive[1][3] = 0.5
    ####################
    # Zincblende
    ####################
    if(structure=="ZB"):
      primitive = [[0 for x in range(4)] for y in range(8)]
      primitive[0][0] = str(symbol)
      primitive[0][1] = 0.0
      primitive[0][2] = 0.0
      primitive[0][3] = 0.0
      primitive[1][0] = str(symbol)
      primitive[1][1] = 0.5
      primitive[1][2] = 0.5
      primitive[1][3] = 0.0
      primitive[2][0] = str(symbol)
      primitive[2][1] = 0.5
      primitive[2][2] = 0.0
      primitive[2][3] = 0.5
      primitive[3][0] = str(symbol)
      primitive[3][1] = 0.0
      primitive[3][2] = 0.5
      primitive[3][3] = 0.5
      primitive[4][0] = str(symbol)
      primitive[4][1] = 0.25
      primitive[4][2] = 0.25
      primitive[4][3] = 0.25
      primitive[5][0] = str(symbol)
      primitive[5][1] = 0.75
      primitive[5][2] = 0.75
      primitive[5][3] = 0.25
      primitive[6][0] = str(symbol)
      primitive[6][1] = 0.25
      primitive[6][2] = 0.75
      primitive[6][3] = 0.75
      primitive[7][0] = str(symbol)
      primitive[7][1] = 0.75
      primitive[7][2] = 0.25
      primitive[7][3] = 0.75
####return
    return primitive

  @staticmethod
  def buildStructure(primitive, copy, copyB=None, copyC=None, perturb=None):
    # Handle arguments
    copyX = copy
    if copyB is None:
      copyY = copy
    else:
      copyY = copyB
    if copyC is None:
      copyZ = copy
    else:
      copyZ = copyC
    if perturb is None:
      perturb = 0.0
    # Random number
    rNum = RandomLCG()
    newSize = len(primitive)*copyX*copyY*copyZ
    # Define
    structure = [[0 for x in range(4)] for y in range(newSize)]
    n = 0
    rNumArr = [0]*3
    for i in range (0,copyX):
      for j in range (0,copyY):
        for k in range (0,copyZ):
          for row in primitive:
            for m in range(0,3):
              rNumArr[m] = rNum.rng()
            structure[n][0] = row[0]
            structure[n][1] = (i + row[1]+(perturb*(0.5-rNumArr[0]))) / copyX
            structure[n][2] = (j + row[2]+(perturb*(0.5-rNumArr[1]))) / copyY
            structure[n][3] = (k + row[3]+(perturb*(0.5-rNumArr[2]))) / copyZ
            n = n + 1
    return structure

  @staticmethod
  def printStructure(structure):
    for row in structure:
      for field in row:
        print(field,sep="",end="")
        print("   ",sep="",end="")
      print()












################################################################################
