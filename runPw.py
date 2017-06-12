#!/bin/python3
################################################################################

import os
from pwOut import pwOut as pwOut


class runPw:

  @staticmethod
  def run(fileName, dftDir=None, runCode=1, procs=0, verbose=False):
    if(dftDir is None):
      inFile = fileName+".in"
      outFile = fileName+".out"
    else:
      inFile = dftDir+"/"+fileName+".in"
      outFile = dftDir+"/"+fileName+".out"

    if(procs==0 or procs is None):
      try:
        procs = os.environ['procCount']
      except:
        procs = 4

    cmdIn = "mpirun "
    cmdIn = cmdIn + "-np " + str(procs) + " "
    cmdIn = cmdIn + "pw.x "
    cmdIn = cmdIn + "< " + inFile + " "
    cmdIn = cmdIn + "> " + outFile

    if(verbose):
      print (runCode, cmdIn)

    if(runCode==1):
      # Run anyway
      os.system(cmdIn)
    if(runCode==2):
      # Check if successful
      pwDat = pwOut(outFile)
      if(not pwDat.successful):
        os.system(cmdIn)
    return cmdIn






################################################################################
