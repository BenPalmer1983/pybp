#!/bin/python3
################################################################################

import sys
import os
import hashlib
from randnum import RandomLCG as RandomLCG


###########################
# File Data
###########################
class oFileData:

  def __init__(self):
    self.fileName = None
    self.lineCount = 0
    self.fileData = []
    self.onDisk = False

  def loadFile(self, fileName):
    self.fileName = fileName
    if os.path.isfile(self.fileName):
      self.onDisk = True
      rFile = open(self.fileName, 'r')
      for line in rFile:
        self.fileData.append(line.rstrip('\n'))
        self.lineCount = self.lineCount + 1   # increment counter

  def writeFile(self, fileName=None):
    if(fileName is not None):
      self.fileName = fileName
      self.onDisk = True
    wFile = open(self.fileName, 'w')
    for fileLine in self.fileData:
      wFile.write(str(fileLine)+'\n')

  def printFile(self):
    for i in range(0,self.lineCount):
      print(self.fileData[i])

  def md5hash(self):
    myhash = hashlib.md5()
    message = ""
    for fileLine in self.fileData:
      message = message + fileLine
    myhash.update(message.encode())
    self.data_md5hash = myhash.hexdigest()

###########################
# Strings
###########################
class oStrings:

  @staticmethod
  def removeDouble(inputStr, dChar):
    outputStr = ""
    lastChar = None
    for char in inputStr:
      if(lastChar is not None):
        if(char == dChar):
          if(char != lastChar):
            outputStr = outputStr+char
        else:
          outputStr = outputStr+char
      else:
        outputStr = outputStr+char
      lastChar = char
    return outputStr

  @staticmethod
  def removeSpaces(inputStr):
    outputStr = ""
    for char in inputStr:
      if(char != " "):
        outputStr = outputStr+char
    return outputStr

  @staticmethod
  def trimEnds(inputStr):
    outputStr = ""
    start = 0
    end = len(inputStr)
    for i in range (0,len(inputStr)):
      if(inputStr[i:i+1]!=" "):
        start = i
        break
    for i in range (0,len(inputStr)):
      j = len(inputStr) - (i+1)
      if(inputStr[j:j+1]!=" "):
        end = j+1
        break
    output = inputStr[start:end]
    return output

###########################
# General functions
###########################
class general:

  @staticmethod
  def mkDir(dirIn):
    cmdIn = "mkdir -p "+dirIn
    os.system(cmdIn)










################################################################################
