#!/bin/python3
################################################################################


from bpstandard import oFileData as oFileData
from bpstandard import oStrings as oStrings
from strains import strains as strains
from structures import structures as structures


###########################
# pwscf input file
###########################
class pwIn:

  def __init__(self, fileIn=None):
    if(fileIn is not None):
      self.loadFile(fileIn)

  def loadFile(self, fileIn):
    self.filePath = fileIn
    self.pwFile = oFileData()               # make new oFileData object
    self.pwFile.loadFile(fileIn)            # load from file
    self.extractData()

  def getInput(self, inputVal):
    inputVal = inputVal.replace(",","")
    inputVal = inputVal.split("=")
    outputVal = float(inputVal[1])
    return outputVal

  def extractData(self):
    self.aLat = 0.0e0
    self.cellParameters = [[0 for x in range(3)] for y in range(3)]
    self.nat = 0
    self.ntype = 0
    self.atomSpecies_symbol = []
    self.atomSpecies_mass = []
    self.atomSpecies_pp = []
    for i in range(0,self.pwFile.lineCount):
      if 'CELLDM(1)' in self.pwFile.fileData[i].upper():
        self.aLat = self.getInput(self.pwFile.fileData[i])
      if 'CELL_PARAMETERS' in self.pwFile.fileData[i].upper():
        for j in range(0,3):
          i = i + 1
          fileRow = oStrings.removeDouble(self.pwFile.fileData[i]," ")
          fileRow = fileRow.split(" ")
          for k in range(0,3):
            self.cellParameters[j][k] = float(fileRow[k])
      if 'NTYP' in self.pwFile.fileData[i].upper():
        fileRow = self.pwFile.fileData[i].replace(",", "")
        fileRowArr = fileRow.split("=")
        self.ntype = int(fileRowArr[1])
      if 'NAT' in self.pwFile.fileData[i].upper():
        fileRow = self.pwFile.fileData[i].replace(",", "")
        fileRowArr = fileRow.split("=")
        self.nat = int(fileRowArr[1])
      if 'ATOMIC_SPECIES' in self.pwFile.fileData[i].upper():
        for j in range(0, self.ntype):
          i = i + 1
          fileRow = oStrings.removeDouble(self.pwFile.fileData[i]," ")
          fileRowArr = fileRow.split(" ")
          self.atomSpecies_symbol.append(fileRowArr[0])
          self.atomSpecies_mass.append(fileRowArr[1])
          self.atomSpecies_pp.append(fileRowArr[2])

  def readAtomSpecies(self):
    self.atomSpecies_symbol = []
    self.atomSpecies_mass = []
    self.atomSpecies_pp = []
    for i in range(0,self.pwFile.lineCount):
      if 'ATOMIC_SPECIES' in self.pwFile.fileData[i].upper():
        for j in range(0, self.ntype):
          i = i + 1
          fileRow = oStrings.removeDouble(self.pwFile.fileData[i]," ")
          fileRowArr = fileRow.split(" ")
          self.atomSpecies_symbol.append(fileRowArr[0])
          self.atomSpecies_mass.append(fileRowArr[1])
          self.atomSpecies_pp.append(fileRowArr[2])

  def reloadFile(self):
    fileIn = self.filePath
    self.loadFile(self, fileIn)

  def changeCalculation(self, calculation):
    for i in range(0,self.pwFile.lineCount):
      if 'calculation' in self.pwFile.fileData[i]:
        self.pwFile.fileData[i] = 'calculation = "'+calculation+'",'

  def changeOutdir(self, dirIn):
    for i in range(0,self.pwFile.lineCount):
      if 'outdir' in self.pwFile.fileData[i]:
        self.pwFile.fileData[i] = 'outdir = "'+dirIn+'",'

  def changePPdir(self, dirIn):
    for i in range(0,self.pwFile.lineCount):
      if 'pseudo_dir' in self.pwFile.fileData[i]:
        self.pwFile.fileData[i] = 'pseudo_dir = "'+dirIn+'",'

  def changeEcutwfc(self, ecutwfc):
    for i in range(0,self.pwFile.lineCount):
      if 'ECUTWFC' in self.pwFile.fileData[i].upper():
        self.pwFile.fileData[i] = 'ecutwfc = '+str(ecutwfc)+','

  def changeEcutrho(self, ecutrho):
    for i in range(0,self.pwFile.lineCount):
      if 'ECUTRHO' in self.pwFile.fileData[i].upper():
        self.pwFile.fileData[i] = 'ecutrho = '+str(ecutrho)+','

  def changeDegauss(self, degauss):
    for i in range(0,self.pwFile.lineCount):
      if 'DEGAUSS' in self.pwFile.fileData[i].upper():
        self.pwFile.fileData[i] = 'degauss = '+str(degauss)+','

  def changeAtomSpecies(self, atomLabels, atomMasses, atomPPs):
    fileData_new = []
    j = 0
    for i in range(0,self.pwFile.lineCount):
      if(i>=self.pwFile.lineCount):
        break
      if(j>=self.pwFile.lineCount):
        break
      if 'ATOMIC_SPECIES' in self.pwFile.fileData[j].upper():
        fileData_new.append(self.pwFile.fileData[j])
        for k in range(j,self.pwFile.lineCount):
          j = j + 1
          if 'ATOMIC_POSITIONS' in self.pwFile.fileData[j].upper():
            break
          if self.pwFile.fileData[j] == "":
            break
        for k in range(0,len(atomLabels)):
          atomRow = atomLabels[k]+"   "+str(atomMasses[k])+"   "+atomPPs[k]
          fileData_new.append(atomRow)
      else:
        fileData_new.append(self.pwFile.fileData[j])
        j = j + 1
    self.pwFile.fileData = fileData_new
    self.pwFile.lineCount = len(fileData_new)
    # Update number of atom types
    self.changeNtyp(len(atomLabels))
    # Extract data
    self.extractData()

  #def changeAtomSpecies(self, atomLabels, atomMasses, atomPPs):

  def changeKpoints(self, kpoints):
    # Check input
    kpoints = str(kpoints)
    kpoints = oStrings.removeDouble(kpoints," ")
    kpoints = oStrings.trimEnds(kpoints)
    kpoints = kpoints
    kpointsArr = kpoints.split(" ")
    print (len(kpointsArr))
    if(len(kpointsArr)==6):
      kpoints_row = str(kpointsArr[0])+" "+str(kpointsArr[1])+" "+str(kpointsArr[2])+" "+"   1 1 1"
    if(len(kpointsArr)==3):
        kpoints_row = str(kpointsArr[0])+" "+str(kpointsArr[1])+" "+str(kpointsArr[2])+" "+"   1 1 1"
    if(len(kpointsArr)==1):
      kpoints_row = str(kpointsArr[0])+" "+str(kpointsArr[0])+" "+str(kpointsArr[0])+"   1 1 1"

    for i in range(0,self.pwFile.lineCount):
      if 'K_POINTS' in self.pwFile.fileData[i].upper():
        self.pwFile.fileData[i] = 'K_POINTS automatic'
        self.pwFile.fileData[i+1] = kpoints_row

  def changeAlat(self, aLat, updateValue=False):
    for i in range(0,self.pwFile.lineCount):
      if 'celldm(1)' in self.pwFile.fileData[i]:
        self.pwFile.fileData[i] = 'celldm(1) = '+str(aLat)+','
    if updateValue == True:
      self.aLat = float(aLat)

  def changeCellParameters(self, newCell):
    for i in range(0,self.pwFile.lineCount):
      testStr = self.pwFile.fileData[i];
      testStr = testStr.upper()
      # Update data file
      if 'CELL_PARAMETERS' in testStr:
        for j in range(0,3):
          i = i + 1
          self.pwFile.fileData[i] = str(newCell[j][0])+" "+str(newCell[j][1])+" "+str(newCell[j][2])
      # Update cell
      self.cellParameters = newCell

  def changePrefix(self, prefix):
    for i in range(0,self.pwFile.lineCount):
      testStr = self.pwFile.fileData[i];
      testStr = testStr.upper()
      # Update data file
      if ('PREFIX =' in testStr or 'PREFIX=' in testStr):
        self.pwFile.fileData[i] = "prefix = \""+str(prefix)+"\","

  def changeNat(self, atomCount):
    for i in range(0,self.pwFile.lineCount):
      testStr = self.pwFile.fileData[i];
      testStr = testStr.upper()
      # Update data file
      if ('NAT =' in testStr or 'NAT=' in testStr):
        self.pwFile.fileData[i] = "nat = "+str(atomCount)+","

  def changeNtyp(self, atomTypes):
    for i in range(0,self.pwFile.lineCount):
      testStr = self.pwFile.fileData[i];
      testStr = testStr.upper()
      # Update data file
      if ('NTYP =' in testStr or 'NTYP=' in testStr):
        self.pwFile.fileData[i] = "ntyp = "+str(atomTypes)+","

  def changeAtoms(self, newStructure):
    tempFile = oFileData()
    tempFile.fileName = self.pwFile.fileName
    tempFile.onDisk = self.pwFile.onDisk
    n = 0
    skip = 0
    for row in self.pwFile.fileData:
      testStr = row.upper()
      rowWritten = False
      if testStr[0:3] == 'NAT':
        n = n + 1
        rowWritten = True
        tempFile.fileData.append("nat = "+str(len(newStructure))+",")
      if testStr[0:16] == 'ATOMIC_POSITIONS':
        n = n + 1
        rowWritten = True
        tempFile.fileData.append("ATOMIC_POSITIONS crystal")
        for j in range(0,len(newStructure)):
          n = n + 1
          newStr = str(newStructure[j][0])+"     "+str(newStructure[j][1])+"     "
          newStr = newStr + str(newStructure[j][2])+"     "+str(newStructure[j][3])
          tempFile.fileData.append(newStr)
        skip = self.nat
      if rowWritten == False:
        if skip == 0:
          n = n + 1
          tempFile.fileData.append(row)
        else:
          skip = skip - 1
    tempFile.lineCount = n
    self.pwFile = tempFile   # overwrite file
    #for row in self.pwFile.fileData:
    #  print (row)
    self.extractData()

  def setMixingTF(self):
    rowCounter = 0
    for i in range(0,self.pwFile.lineCount):
      testStr = self.pwFile.fileData[i];
      testStr = testStr.upper()
      # Update data file
      if ('MIXING_MODE' in testStr):
        self.pwFile.fileData[i] = "mixing_mode = 'TF',"

  def setMixingLocalTF(self):
    rowCounter = 0
    for i in range(0,self.pwFile.lineCount):
      testStr = self.pwFile.fileData[i];
      testStr = testStr.upper()
      # Update data file
      if ('MIXING_MODE' in testStr):
        self.pwFile.fileData[i] = "mixing_mode = 'local-TF',"

  def setMagnetic0(self):
    rowCounter = 0
    for i in range(0,self.pwFile.lineCount):
      testStr = self.pwFile.fileData[i];
      testStr = testStr.upper()
      # Update data file
      if ('STARTING_MAGNETIZATION' in testStr):
        self.removeLine(rowCounter)
      else:
        rowCounter = rowCounter + 1
      if(rowCounter==self.pwFile.lineCount-1):
        break

  def setMagnetic2(self):
    nspinSet = False
    self.readAtomSpecies()
    for i in range(0,self.pwFile.lineCount):
      testStr = self.pwFile.fileData[i];
      testStr = testStr.upper()
      # Update data file
      if ('NSPIN =' in testStr or 'NSPIN=' in testStr):
        self.pwFile.fileData[i] = "nspin = 2,"
        nspinSet = True
    if(not nspinSet):
      for i in range(0,self.pwFile.lineCount):
        testStr = self.pwFile.fileData[i];
        testStr = testStr.upper()
        if ('NTYP =' in testStr or 'NTYP=' in testStr):
          self.insertLine(i+1,"nspin = 2,")
          for j in range(0,len(self.atomSpecies_symbol)):
            self.insertLine(i+j+2,"starting_magnetization("+str(j+1)+") = 0.3")
          break

  def makeStructure(self,structure,copyX=1,copyY=None,copyZ=None,perturb=None):
    if(copyY is None):
      copyY = copyX
    if(copyZ is None):
      copyZ = copyX
    if(perturb is None):
      perturb = 0.0
    self.readAtomSpecies()
    primitive = structures.primitiveStructure(structure,self.atomSpecies_symbol)
    expanded = structures.buildStructure(primitive, copyX, copyY, copyZ, perturb)
    self.changeAtoms(expanded)

  def insertLine(self,fileRow,line):
    # Append extra line
    self.pwFile.fileData.append("")
    self.pwFile.lineCount = self.pwFile.lineCount  + 1
    # Move all rows down
    row = self.pwFile.lineCount-2
    while(row>=fileRow):
      self.pwFile.fileData[row+1] = self.pwFile.fileData[row]
      row = row - 1
    self.pwFile.fileData[fileRow] = line

  def removeLine(self,fileRow):
    tempList = []
    for i in range(0,self.pwFile.lineCount):
      if(i!=fileRow):
        tempList.append(self.pwFile.fileData[i])
    self.pwFile.fileData = []
    for i in range(0,len(tempList)):
      self.pwFile.fileData.append(tempList[i])
    self.pwFile.lineCount = self.pwFile.lineCount - 1

  def normaliseParameters(self):
# Set the cell_parameter matrix such that (1,1) = 1.0 and adjust aLat accordingly
    aLat = 0.0e0
    cellParameters = [[0 for x in range(3)] for y in range(3)]
    for i in range(0,self.pwFile.lineCount):
      if 'celldm(1)' in self.pwFile.fileData[i]:
        aLat = self.getInput(self.pwFile.fileData[i])
      if 'CELL_PARAMETERS' in self.pwFile.fileData[i]:
        for j in range(0,3):
          i = i + 1
          fileRow = oStrings.removeDouble(self.pwFile.fileData[i]," ")
          fileRow = fileRow.split(" ")
          for k in range(0,3):
            cellParameters[j][k] = float(fileRow[k])
    topVal = cellParameters[0][0]
    # Change aLat
    aLat = aLat * topVal
    # Change cell
    for j in range(0,3):
      for k in range(0,3):
        if(j==0 and k==0):
          cellParameters[j][k] = 1.0
        else:
          cellParameters[j][k] = cellParameters[j][k] / topVal
    self.changeAlat(aLat)
    self.changeCellParameters(cellParameters)


  def fileName(self, fileNameIn):
    self.pwFile.fileName = fileNameIn

  def outputFile(self, fileNameIn = None):
    if(fileNameIn is not None):
      self.pwFile.fileName = fileNameIn
    self.pwFile.writeFile()

  def displayAlat(self):
    print ()
    print("aLat:                  ",self.aLat)
    for i in range(0,3):
      #print("Cell Parameters:       ",end="")
      for j in range(0,3):
        #print(self.cellParameters[i][j],end=" ")
        print(self.cellParameters[i][j])
      print ()
    print ()

  def displayFile(self):
    self.pwFile.printFile()


  def transform(self, transformVec, inputUnitVec=None):
    output = [[0 for x in range(3)] for y in range(3)]
    cpFlag = False
    if(inputUnitVec is None):
      cpFlag = True
      inputUnitVec = self.cellParameters
    for i in range (0,3):
      for j in range (0,3):
        output[i][j] = 0.0e0
    for i in range (0,3):
      for j in range (0,3):
        for k in range(0,3):
          output[i][j] = output[i][j] + transformVec[i][k] * inputUnitVec[k][j]
    if cpFlag is True:
      self.cellParameters = output        # update matrix
      self.changeCellParameters(output)   # update file data
    return output

  def scale(self, sVal, inputUnitVec=None):
    transformVec = [[0 for x in range(3)] for y in range(3)]
    cpFlag = False
    if(inputUnitVec is None):
      cpFlag = True
      inputUnitVec = self.cellParameters
    for i in range (0,3):
      for j in range (0,3):
        if (i == j):
          transformVec[i][j] = float(sVal)
        else:
          transformVec[i][j] = 0.0e0
    output = self.transform(transformVec, inputUnitVec)
    if cpFlag is True:
      self.cellParameters = output        # update matrix
      self.changeCellParameters(output)   # update file data

  def resetCellParameters(self):
    iMatrix = strains.identityMatrix()
    self.cellParameters = iMatrix        # update matrix
    self.changeCellParameters(iMatrix)   # update file data


  ###############################################
  # Specific templates
  ###############################################

  def makeIsolated(self):
    tempFile = oFileData()
    tempFile.fileName = "isolated.in"
    tempFile.onDisk = False

    # Control
    tempFile.fileData.append("&CONTROL")
    tempFile.fileData.append("restart_mode = 'from_scratch',")
    tempFile.fileData.append("calculation = \"scf\",")
    tempFile.fileData.append("etot_conv_thr = 1.0E-4,")
    tempFile.fileData.append("forc_conv_thr = 1.0D-3,")
    tempFile.fileData.append("nstep = 40,")
    tempFile.fileData.append("tprnfor=.true.,")
    tempFile.fileData.append("tstress=.true.,")
    tempFile.fileData.append("disk_io='low',")
    tempFile.fileData.append("prefix=\"isolatedAtom\",")
    tempFile.fileData.append("outdir = \"\",")
    tempFile.fileData.append("pseudo_dir = \"\",")
    tempFile.fileData.append("/")
    tempFile.fileData.append("")

    # System
    tempFile.fileData.append("&SYSTEM")
    tempFile.fileData.append("ibrav = 0,")
    tempFile.fileData.append("celldm(1) = 20.0,")
    tempFile.fileData.append("nat = 1,")
    tempFile.fileData.append("ntyp = 1,")
    tempFile.fileData.append("ecutwfc = 50,")
    tempFile.fileData.append("ecutrho = 200,")
    tempFile.fileData.append("occupations = 'smearing',")
    tempFile.fileData.append("smearing = 'mv',")
    tempFile.fileData.append("degauss = 0.02,")
    tempFile.fileData.append("nosym = .TRUE.,")
    tempFile.fileData.append("starting_magnetization(1) = 0.7,")
    tempFile.fileData.append("nspin = 2,")
    tempFile.fileData.append("/")
    tempFile.fileData.append("")

    # Electrons
    tempFile.fileData.append("&ELECTRONS")
    tempFile.fileData.append("mixing_beta = 3.0000000E-01,")
    tempFile.fileData.append("mixing_ndim = 10,")
    tempFile.fileData.append("diagonalization='david',")
    tempFile.fileData.append("mixing_mode = 'TF',")
    tempFile.fileData.append("conv_thr = 1.0D-6,")
    tempFile.fileData.append("/")
    tempFile.fileData.append("")

    # Ions
    tempFile.fileData.append("&IONS")
    tempFile.fileData.append("ion_dynamics='bfgs',")
    tempFile.fileData.append("/")
    tempFile.fileData.append("")

    # Cell
    tempFile.fileData.append("&CELL")
    tempFile.fileData.append("cell_dynamics='bfgs',")
    tempFile.fileData.append("press = 0.0,")
    tempFile.fileData.append("cell_factor = 2.0,")
    tempFile.fileData.append("/")
    tempFile.fileData.append("")

    # Atomic Species
    tempFile.fileData.append("ATOMIC_SPECIES")
    tempFile.fileData.append("n 50.0 n.pbe-nl-kjpaw_psl.1.0.0.UPF")
    tempFile.fileData.append("")

    # Atomic positions
    tempFile.fileData.append("ATOMIC_POSITIONS crystal")
    tempFile.fileData.append("n     0.5     0.5     0.5")
    tempFile.fileData.append("")

    # K Points
    tempFile.fileData.append("K_POINTS automatic")
    tempFile.fileData.append("10 10 10    1 1 1")
    tempFile.fileData.append("")

    # Cell Parameters
    tempFile.fileData.append("CELL_PARAMETERS alat")
    tempFile.fileData.append("1.0    0.0    0.0")
    tempFile.fileData.append("0.0    1.0    0.0")
    tempFile.fileData.append("0.0    0.0    1.0")
    tempFile.fileData.append("")

    # Line count
    tempFile.lineCount = len(tempFile.fileData)

    # Store file and extract data
    self.pwFile = tempFile              # make new oFileData object
    self.extractData()


  def makeTemplate(self):
    tempFile = oFileData()
    tempFile.fileName = "template.in"
    tempFile.onDisk = False

    # Control
    tempFile.fileData.append("&CONTROL")
    tempFile.fileData.append("restart_mode = 'from_scratch',")
    tempFile.fileData.append("calculation = \"scf\",")
    tempFile.fileData.append("etot_conv_thr = 1.0E-4,")
    tempFile.fileData.append("forc_conv_thr = 1.0D-3,")
    tempFile.fileData.append("nstep = 40,")
    tempFile.fileData.append("tprnfor = .true.,")
    tempFile.fileData.append("tstress = .true.,")
    tempFile.fileData.append("disk_io = 'low',")
    tempFile.fileData.append("prefix = \"template\",")
    tempFile.fileData.append("outdir = \"\",")
    tempFile.fileData.append("pseudo_dir = \"\",")
    tempFile.fileData.append("/")
    tempFile.fileData.append("")

    # System
    tempFile.fileData.append("&SYSTEM")
    tempFile.fileData.append("ibrav = 0,")
    tempFile.fileData.append("celldm(1) = 20.0,")
    tempFile.fileData.append("nat = 0,")
    tempFile.fileData.append("ntyp = 0,")
    tempFile.fileData.append("ecutwfc = 50,")
    tempFile.fileData.append("ecutrho = 200,")
    tempFile.fileData.append("occupations = 'smearing',")
    tempFile.fileData.append("smearing = 'mv',")
    tempFile.fileData.append("degauss = 0.02,")
    #tempFile.fileData.append("nosym = .TRUE.,")
    #tempFile.fileData.append("starting_magnetization(1) = 0.7,")
    #tempFile.fileData.append("nspin = 2,")
    tempFile.fileData.append("/")
    tempFile.fileData.append("")

    # Electrons
    tempFile.fileData.append("&ELECTRONS")
    tempFile.fileData.append("mixing_beta = 3.0000000E-01,")
    tempFile.fileData.append("mixing_ndim = 10,")
    tempFile.fileData.append("diagonalization='david',")
    tempFile.fileData.append("mixing_mode = 'TF',")
    tempFile.fileData.append("conv_thr = 1.0D-6,")
    tempFile.fileData.append("/")
    tempFile.fileData.append("")

    # Ions
    tempFile.fileData.append("&IONS")
    tempFile.fileData.append("ion_dynamics='bfgs',")
    tempFile.fileData.append("/")
    tempFile.fileData.append("")

    # Cell
    tempFile.fileData.append("&CELL")
    tempFile.fileData.append("cell_dynamics='bfgs',")
    tempFile.fileData.append("press = 0.0,")
    tempFile.fileData.append("cell_factor = 2.0,")
    tempFile.fileData.append("/")
    tempFile.fileData.append("")

    # Atomic Species
    tempFile.fileData.append("ATOMIC_SPECIES")
    tempFile.fileData.append("")

    # Atomic positions
    tempFile.fileData.append("ATOMIC_POSITIONS crystal")
    tempFile.fileData.append("")

    # K Points
    tempFile.fileData.append("K_POINTS automatic")
    tempFile.fileData.append("10 10 10    1 1 1")
    tempFile.fileData.append("")

    # Cell Parameters
    tempFile.fileData.append("CELL_PARAMETERS alat")
    tempFile.fileData.append("1.0    0.0    0.0")
    tempFile.fileData.append("0.0    1.0    0.0")
    tempFile.fileData.append("0.0    0.0    1.0")
    tempFile.fileData.append("")

    # Line count
    tempFile.lineCount = len(tempFile.fileData)

    # Store file and extract data
    self.pwFile = tempFile              # make new oFileData object
    self.extractData()


################################################################################
