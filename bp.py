#!/bin/python3
################################################################################
#
# Calculate Bulk Property
#
# Elastic constants: method from First principles calculations of elastic
# properties of metals 1994  Mehl, Klien, Papaconstantopoulos
#

#
import getopt
import math
import numpy
import copy
import time
#
import os
import sys
#include = os.environ['PYLIB']
#sys.path.append(include)

from bpstandard import oFileData as oFileData
from bpstandard import oStrings as oStrings
from bpstandard import general as general
from structures import structures as structures
from pwIn import pwIn as pwIn
from pwOut import pwOut as pwOut
from pwOut import pwCompare as pwCompare
from runPw import runPw as runPw
from gnuplot import gnuplot as gnuplot
from strains import strains as strains
from eos import eos as eos
from randnum import RandomLCG as RandomLCG


###########################
# Input/Output
###########################
class oInOut:
  @staticmethod
  def readInputFile():
    print("Read")




###########################
# bpCalc
###########################
class bpCalc:
  iCount = 0  # Instance counter

  def __init__(self):
    bpCalc.iCount = bpCalc.iCount + 1   # Increment instance counter
    self.cmdLog = []
    self.resultsArr = []
    self.startTime = time.time()
    #self.resultsArr.append(self.startTime)
    # Init files
    self.tFile = oFileData()
    self.wFile = oFileData()
    self.dFile = oFileData()
    # Load control file
    self.controlFile()
    # Load pwscf input file
    #self.tFile.loadFile(self.tFileName)
    # Store input data
    #self.pwIn = pwIn(self.tFileName)
    #self.pwIsoIn = pwIn(self.tFileName)
    # Change pwIn file
    #self.pwIn.changeEcutwfc(self.ecutwfc)
    #self.pwIn.changeEcutrho(self.ecutrho)
    #self.pwIn.changeKpoints(self.kpoints)
    # Change pwIsoIn file
    #self.pwIsoIn.changeEcutwfc(self.ecutwfc)
    #self.pwIsoIn.changeEcutrho(self.ecutrho)
    #self.pwIsoIn.changeKpoints(self.kpoints)

  def controlFile(self):
    print ("Reading control file")
    print ("==================================================")
    self.args = str(sys.argv)   # Store input argument as the templateFile
    self.cFileName = sys.argv[1]
    self.cFile = oFileData()
    self.cFile.loadFile(self.cFileName)
    self.tmpDir = "tmp"
    self.runType = "full"
    self.copies = 2
    self.bmstrain = 0.05
    self.orstrain = 0.02
    self.testrain = 0.02
    self.bmsteps = 10
    self.orsteps = 10
    self.testeps = 10
    self.outdir = ""
    self.ppdir = ""
    self.pwtemplate = None

    self.atomList = []
    self.massList = []
    self.ppList = []

    for i in range(0,self.cFile.lineCount):
      # Directories
      #######################
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$TMPDIR", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.tmpDir = fileRowArr[1]
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$OUTDIR", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.outdir = fileRowArr[1]
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$PPDIR", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.ppdir = fileRowArr[1]

      # PWscf Settings
      #######################
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$PWTEMPLATE", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.pwtemplate = fileRowArr[1]
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ATOMLIST", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        for i in range(1,len(fileRowArr)):
          self.atomList.append(fileRowArr[i])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$MASSLIST", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        for i in range(1,len(fileRowArr)):
          self.massList.append(fileRowArr[i])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$PPLIST", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        for i in range(1,len(fileRowArr)):
          self.ppList.append(fileRowArr[i])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$STRUCTURE", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.structure = fileRowArr[1]
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ALAT", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.inputAlat = float(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$COPIES", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.copies = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$NSPIN", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.nspin = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$MIXING_MODE", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.mixing_mode = fileRowArr[1]


      # Run Settings
      #######################
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$PROCS", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.procCount = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$CONVERGE", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.convergeChoice = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$RUN", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.runCode = int(fileRowArr[1])


      # Convergence Settings
      #######################
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ETHRESHOLD", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.eThreshold = float(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$FTHRESHOLD", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.fThreshold = float(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$STHRESHOLD", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.sThreshold = float(fileRowArr[1]) / 1.471050658e7
      ####
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$MINECUTWFC", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.min_ecutwfc = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ECUTWFC", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.ecutwfc = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ECUTRHO", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.ecutrho = int(fileRowArr[1])
      ####
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$KPOINTS", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.kpoints = str(fileRowArr[1])+" "+str(fileRowArr[2])+" "
        self.kpoints = self.kpoints + str(fileRowArr[3])+" "+str(fileRowArr[4])+" "
        self.kpoints = self.kpoints + str(fileRowArr[5])+" "+str(fileRowArr[6])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$STARTKPOINTS", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.kPointStart = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ENDKPOINTS", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.kPointEnd = int(fileRowArr[1])


      # Convergence Settings
      #######################
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$BMSTRAIN", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.bmstrain = float(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ORSTRAIN", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.orstrain = float(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$TESTRAIN", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.testrain = float(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$BMSTEPS", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.bmsteps = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$ORSTEPS", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.orsteps = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$TESTEPS", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.testeps = int(fileRowArr[1])







    print ("==================================================")
    print()

    print(self.atomList)

    ## Mk dir
    general.mkDir(self.tmpDir)
    ## Type of run
    self.runCode = 0 # Default
    if(self.runType.upper()=="FULL" or self.runType[0:1]=="1"):       # Full run from fresh
      self.runCode = 1
      self.runType = "FULL"
    if(self.runType.upper()=="CONTINUE" or self.runType[0:1]=="2"):   # Completes where there's no output file
      self.runCode = 2
      self.runType = "CONTINUE"
    if(self.runType.upper()=="NONE" or self.runType[0:1]=="3"):       # Don't run DFT
      self.runCode = 0

    ## Log
    self.resultsArr.append("Run Options     ")
    self.resultsArr.append("==========================================")
    self.resultsArr.append("Dir:       "+self.tmpDir)
    self.resultsArr.append("Code:      "+str(self.runCode)+"  ("+str(self.runType)+")")
    self.resultsArr.append("Structure: "+str(self.structure))
    self.resultsArr.append("  ")


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


  def replaceInput(self, fileName):
    self.pwIn.outputFile(fileName)

#############################
#############################

  def run(self):
    print ("Run")
    self.makeStructure()
    self.relaxed()
    self.isolatedAtom()
    #self.unperturbed()
    #self.bulkModulus()
    #self.cubicElasticConstants()
    #self.outputResults()

  def makeStructure(self):       # Use vc-relax to fine optimum settings
    print ("1. Make atom structure")
    self.pwIn = pwIn()
    if(self.pwtemplate is None):
      # Make structure
      print ("   build structure")
      self.primitive = structures.primitiveStructure(self.structure,self.atomList)
      self.expanded = structures.buildStructure(self.primitive,self.copies)
      # Make template file
      self.pwIn.makeTemplate()
      self.pwIn.changeAlat(self.copies * self.inputAlat, True)
      self.pwIn.changeAtoms(self.expanded)
      self.pwIn.changeAtomSpecies(self.atomList,self.massList,self.ppList)
      self.pwIn.changePrefix("template")
      if(self.nspin==0):
        self.pwIn.setMagnetic0()
      if(self.nspin==2):
        self.pwIn.setMagnetic2()
      if(self.mixing_mode=="TF"):
        self.pwIn.setMixingTF()
      if(self.mixing_mode=="local-TF"):
        self.pwIn.setMixingLocalTF()
    else:
      # Load all from a template file
      print ("   load from file")
      self.pwIn.loadFile(self.pwtemplate)
    # Change values (regardless of template file or built)
    self.pwIn.changeEcutwfc(self.ecutwfc)
    self.pwIn.changeEcutrho(self.ecutrho)
    self.pwIn.changeKpoints(self.kpoints)
    self.pwIn.changeOutdir(self.outdir)
    self.pwIn.changePPdir(self.ppdir)
    self.pwIn.outputFile(self.tmpDir+"/template.in")

    #structures.printStructure(self.expanded)

  def relaxed(self):       # Use vc-relax to fine optimum settings
    print ("2. Relaxed aLat")
    # Set up input file
    self.pwIn.changeCalculation("vc-relax")
    self.pwIn.outputFile(self.tmpDir+"/vc-relax.in")
    # Run pwscf
    runPw.run("vc-relax",self.tmpDir,self.runCode)
    self.pwRelaxed = pwOut(self.tmpDir+"/vc-relax.out")
    aLat = self.pwRelaxed.getAlat_Relaxed()   # Read aLat from relaxed output file
    self.aLat_relaxed = aLat
    rCell = self.pwRelaxed.getCell_Relaxed()   # Read cell from relaxed output file
    self.pwIn.changeAlat(aLat)                # Update aLat in input file
    self.pwIn.changeCellParameters(rCell)     # Update cell parameters in input file
    self.pwIn.normaliseParameters()           # normalise cell/aLat
    self.pwIn.changeCalculation("scf")
    self.replaceInput(self.tmpDir+"/opt.in")

  def isolatedAtom(self):       # Use vc-relax to fine optimum settings
    print ("3. Isolated Atom")
    atomLabels = []
    atomMasses = []
    atomPPs = []
    atomLabels.append("AL")
    atomMasses.append("26.982")
    atomPPs.append("Al.pbe-nl-kjpaw_psl.1.0.0.UPF")

    # Load relaxed file
    pwRelaxed = pwOut(self.tmpDir+"/vc-relax.out")
    aLat = pwRelaxed.getAlat_RelaxedNormalised()

    primitive = structures.primitiveStructure("ISO",atomLabels)
    pwIn_iso = pwIn()
    pwIn_iso.makeIsolated()
    pwIn_iso.changeAlat(aLat)
    pwIn_iso.changeOutdir(self.outdir)
    pwIn_iso.changePPdir(self.ppdir)
    pwIn_iso.changeAtomSpecies(atomLabels,atomMasses,atomPPs)
    pwIn_iso.changeAtoms(primitive)
    pwIn_iso.outputFile(self.tmpDir+"/isolated.in")
    runPw.run("isolated",self.tmpDir,self.runCode,0,True)
    pwOut_iso = pwOut(self.tmpDir+"/isolated.out")
    self.isolated_energy = pwOut_iso.getEnergyPerAtom()

  def unperturbed(self):       # Use vc-relax to fine optimum settings
    print ("4. Unperturbed")
    print ("     Input file: ",self.tmpDir+"/opt.in")
    # Load file
    pwIn_up = pwIn(self.tmpDir+"/opt.in")
    # Ensure scf
    pwIn_up.changeCalculation("scf")
    # Run pwscf
    pwIn_up.outputFile(self.tmpDir+"/unperturbed.in")
    runPw.run("unperturbed",self.tmpDir,self.runCode,0,True)
    pwOut_up = pwOut(self.tmpDir+"/unperturbed.out")
    self.unperturbed_energy = pwOut_up.getEnergyPerAtom()
    self.unperturbed_volume = pwOut_up.getVolPerAtom()
    self.cohesive_energy = self.unperturbed_energy - self.isolated_energy
    print ("     Energy per atom: ",self.unperturbed_energy)
    print ("     Volume per atom: ",self.unperturbed_volume)

  def bulkModulus(self):
    print("5. Bulk Modulus")
    bmVol = []
    bmEnergy = []
    pwIn_bm = pwIn(self.tmpDir+"/opt.in")
    pwIn_bm.changeCalculation("scf")
    iMatrix = strains.identityMatrix()
    # file name stub
    stub = "bm"
    # Loop through and make input files
    for i in range (0,2*self.bmsteps+1):
    # Transform
      strain = (i-self.bmsteps) * (self.bmstrain / self.bmsteps)
      pwIn_bm.resetCellParameters()
      pwIn_bm.scale(1.00 + strain)
      if((i-self.bmsteps)==0):
        bmVol.append(self.unperturbed_volume)
        bmEnergy.append(self.unperturbed_energy)
      else:
        pwIn_bm.outputFile(self.tmpDir+"/"+stub+"_"+str(i)+".in")
        runPw.run(stub+"_"+str(i),self.tmpDir,self.runCode,0,True)
        pwOut_bm = pwOut(self.tmpDir+"/"+stub+"_"+str(i)+".out")
        bmVol.append(pwOut_bm.getVolPerAtom())
        bmEnergy.append(pwOut_bm.getEnergyPerAtom())
    self.eos = eos()
    self.eos.loadData(bmVol, bmEnergy)
    self.eos.fitEoS()
    #self.eos.display()

    gp = gnuplot()
    gp.reset()
    gp.setDir(self.tmpDir+"/plots")
    gp.title("Equation of State")
    gp.axisLabel("x1","Volume/Atom (Bohr^3)")
    gp.axisLabel("y1","Energy/Atom (Ry)")
    gp.outputPlot("eos")
    gp.addPlot(bmVol,bmEnergy,"x1y1","EoS",1,1)
    gp.makePlot()



  def cubicElasticConstants(self):
    print("6. Cubic Elastic Constants")
    self.cubicElasticConstants_Orthorhombic()
    self.cubicElasticConstants_Tetragonal()
    self.c12 = self.eos.eosO.B0 - self.ortho / 3
    self.c11 = self.ortho + self.c12
    self.c44 = self.tetra

    self.c11_GPa = self.c11 * 1.47105e4
    self.c12_GPa = self.c12 * 1.47105e4
    self.c44_GPa = self.c44 * 1.47105e4

  def cubicElasticConstants_Orthorhombic(self):
    print("6a. Cubic Elastic Constants - Orthorhombic")
    energy = []
    strain = []
    stub = "ortho"
    for i in range (0,self.orsteps+1):
    # Load opt file
      pwIn_ec = pwIn(self.tmpDir+"/opt.in")
    # strain amount
      orthoStrain = i * (self.orstrain / self.orsteps)
    # Transform
      sVec = strains.orthorhombic(orthoStrain)
      pwIn_ec.transform(sVec)
      if(i==0):
        energy.append(self.unperturbed_energy)
        strain.append(0.0)
      else:
        strain.append(orthoStrain)
        pwIn_ec.outputFile(self.tmpDir+"/"+stub+"_"+str(i)+".in")
        runPw.run(stub+"_"+str(i),self.tmpDir,self.runCode,0,True)
        pwOut_ec = pwOut(self.tmpDir+"/"+stub+"_"+str(i)+".out")
        energy.append(pwOut_ec.getEnergyPerAtom())
    pFit = numpy.polyfit(strain, energy, 2)
    self.ortho = pFit[0] / self.unperturbed_volume
    self.ortho_GPa = self.ortho * 1.47105e4
    print(pFit[0],self.ortho_GPa)

    gp = gnuplot()
    gp.reset()
    gp.setDir(self.tmpDir+"/plots")
    gp.title("Orthogonal Strain")
    gp.axisLabel("x1","Strain")
    gp.axisLabel("y1","Energy/Atom (Ry)")
    gp.outputPlot("ortho")
    gp.addPlot(strain,energy,"x1y1","EoS",1,1)
    gp.makePlot()

  def cubicElasticConstants_Tetragonal(self):
    print("6a. Cubic Elastic Constants - Tetragonal")
    energy = []
    strain = []
    stub = "tetra"
    for i in range (0,self.orsteps+1):
    # Load opt file
      pwIn_ec = pwIn(self.tmpDir+"/opt.in")
    # strain amount
      tetraStrain = i * (self.testrain / self.testeps)
    # Transform
      sVec = strains.tetragonal(tetraStrain)
      pwIn_ec.transform(sVec)
      if(i==0):
        energy.append(self.unperturbed_energy)
        strain.append(0.0)
      else:
        strain.append(tetraStrain)
        pwIn_ec.outputFile(self.tmpDir+"/"+stub+"_"+str(i)+".in")
        runPw.run(stub+"_"+str(i),self.tmpDir,self.runCode,0,True)
        pwOut_ec = pwOut(self.tmpDir+"/"+stub+"_"+str(i)+".out")
        energy.append(pwOut_ec.getEnergyPerAtom())
    pFit = numpy.polyfit(strain, energy, 2)
    self.tetra = (2*pFit[0]) / self.unperturbed_volume
    self.tetra_GPa = self.tetra * 1.47105e4

    gp = gnuplot()
    gp.reset()
    gp.setDir(self.tmpDir+"/plots")
    gp.title("Tetragonal Strain")
    gp.axisLabel("x1","Strain")
    gp.axisLabel("y1","Energy/Atom (Ry)")
    gp.outputPlot("tetra")
    gp.addPlot(strain,energy,"x1y1","EoS",1,1)
    gp.makePlot()


  def outputResults(self):
    outputText =              "===================================================\n"
    outputText = outputText + "Results\n"
    outputText = outputText + "===================================================\n"
    outputText = outputText + " \n"
    outputText = outputText + "Unperturbed energy (Ry): "+str(self.unperturbed_energy)+"\n"
    outputText = outputText + "Isolated energy (Ry): "+str(self.isolated_energy)+"\n"
    outputText = outputText + "Cohesive energy (Ry): "+str(self.cohesive_energy)+"  (ev): "+str(self.cohesive_energy*13.6)+"\n"
    outputText = outputText + " \n"
    outputText = outputText + "B0: " + str(self.eos.eosO.B0) + "  B0 (GPA): " + str(self.eos.eosO.B0_GPa)+"\n"
    outputText = outputText + "B0P: " + str(self.eos.eosO.B0P)+"\n"
    outputText = outputText + "E0: " + str(self.eos.eosO.E0)+"\n"
    outputText = outputText + "V0: " + str(self.eos.eosO.V0)+"\n"
    outputText = outputText + "C11 (GPA): " + str(self.c11_GPa)+"\n"
    outputText = outputText + "C12 (GPA): " + str(self.c12_GPa)+"\n"
    outputText = outputText + "C44 (GPA): " + str(self.c44_GPa)+"\n"

    print()
    print(outputText)

    general.mkDir(self.tmpDir+"/summary")
    bpFile = self.tmpDir+"/summary"+"/summary.out"
    outputFile = open(bpFile, 'w')
    outputFile.write(outputText)
    outputFile.close()


########################################################################
## Make BP calc and run it

newBP = bpCalc()
newBP.run()





########################################################################
