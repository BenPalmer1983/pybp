#!/bin/python3
################################################################################
#  Requires:
#    bpstandard.py
#

from bpstandard import RandomLCG as RandomLCG

class structures:
  @staticmethod
  def primitiveStructure(structure, symbol):

    structure = structure.upper()

    ####################
    # ISOLATED
    ####################
    if(structure=="ISO"):
      primitive = [[0 for x in range(4)] for y in range(1)]
      primitive[0][0] = str(symbol[0])
      primitive[0][1] = 0.5
      primitive[0][2] = 0.5
      primitive[0][3] = 0.5
    ####################
    # FCC
    ####################
    if(structure=="FCC"):
      primitive = [[0 for x in range(4)] for y in range(4)]
      primitive[0][0] = str(symbol[(0%len(symbol))])
      primitive[0][1] = 0.0
      primitive[0][2] = 0.0
      primitive[0][3] = 0.0
      primitive[1][0] = str(symbol[(1%len(symbol))])
      primitive[1][1] = 0.5
      primitive[1][2] = 0.5
      primitive[1][3] = 0.0
      primitive[2][0] = str(symbol[(2%len(symbol))])
      primitive[2][1] = 0.5
      primitive[2][2] = 0.0
      primitive[2][3] = 0.5
      primitive[3][0] = str(symbol[(3%len(symbol))])
      primitive[3][1] = 0.0
      primitive[3][2] = 0.5
      primitive[3][3] = 0.5
    ####################
    # BCC
    ####################
    if(structure=="BCC"):
      primitive = [[0 for x in range(4)] for y in range(2)]
      primitive[0][0] = str(symbol[(0%len(symbol))])
      primitive[0][1] = 0.0
      primitive[0][2] = 0.0
      primitive[0][3] = 0.0
      primitive[1][0] = str(symbol[(1%len(symbol))])
      primitive[1][1] = 0.5
      primitive[1][2] = 0.5
      primitive[1][3] = 0.5
    ####################
    # Zincblende
    ####################
    if(structure=="ZB"):
      primitive = [[0 for x in range(4)] for y in range(8)]
      primitive[0][0] = str(symbol[(0%len(symbol))])
      primitive[0][1] = 0.0
      primitive[0][2] = 0.0
      primitive[0][3] = 0.0
      primitive[1][0] = str(symbol[(1%len(symbol))])
      primitive[1][1] = 0.5
      primitive[1][2] = 0.5
      primitive[1][3] = 0.0
      primitive[2][0] = str(symbol[(2%len(symbol))])
      primitive[2][1] = 0.5
      primitive[2][2] = 0.0
      primitive[2][3] = 0.5
      primitive[3][0] = str(symbol[(3%len(symbol))])
      primitive[3][1] = 0.0
      primitive[3][2] = 0.5
      primitive[3][3] = 0.5
      primitive[4][0] = str(symbol[(4%len(symbol))])
      primitive[4][1] = 0.25
      primitive[4][2] = 0.25
      primitive[4][3] = 0.25
      primitive[5][0] = str(symbol[(5%len(symbol))])
      primitive[5][1] = 0.75
      primitive[5][2] = 0.75
      primitive[5][3] = 0.25
      primitive[6][0] = str(symbol[(6%len(symbol))])
      primitive[6][1] = 0.25
      primitive[6][2] = 0.75
      primitive[6][3] = 0.75
      primitive[7][0] = str(symbol[(7%len(symbol))])
      primitive[7][1] = 0.75
      primitive[7][2] = 0.25
      primitive[7][3] = 0.75
####return
    return primitive

  @staticmethod
  def buildStructure(primitive, copy, copyB=None, copyC=None, perturb=None):
    # Handle arguments
    copyX = int(copy)
    if copyB is None:
      copyY = copy
    else:
      copyY = int(copyB)
    if copyC is None:
      copyZ = copy
    else:
      copyZ = int(copyC)
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
            structure[n][1] = ((i + row[1]+(perturb*(0.5-rNumArr[0]))) / copyX) % 1.0
            structure[n][2] = ((j + row[2]+(perturb*(0.5-rNumArr[1]))) / copyY) % 1.0
            structure[n][3] = ((k + row[3]+(perturb*(0.5-rNumArr[2]))) / copyZ) % 1.0
            n = n + 1
    return structure

  @staticmethod
  def printStructure(structure):
    for row in structure:
      for field in row:
        print(field,sep="",end="")
        print("   ",sep="",end="")
      print()
