#!/bin/python3
################################################################################
#  Requires:
#
#

import os


###########################
# Gnuplot
###########################
class gnuplot:

  def __init__(self):
    self.reset()

  def reset(self):
    self.lastPlotID = None
    self.constructed = 1
    self.plotName = "plot"
    self.dir = None
    self.plotTitle = ""
    self.x1Label = None
    self.x2Label = None
    self.y1Label = None
    self.y2Label = None
    self.plots = 1        # Number of individual plots output
    self.dataSets = 0     # Number of data sets
    self.cellParameters = [[0 for y in range(100)] for x in range(1000)]
    self.dataPointCount = [0 for y in range(100)]
    self.dataPointMultiplier = [0 for y in range(100)]
    self.dataAxes = [0 for y in range(100)]
    self.dataName = [0 for y in range(100)]
    self.dataCircles = [0 for y in range(100)]
    self.plotType = [0 for y in range(100)]
    self.maxRows = 0
    for x in range(0,1000):
      for y in range(0,100):
        self.cellParameters[x][y] = None
    for y in range(0,100):
      self.dataPointCount[y] = 0
    for y in range(0,100):
      self.dataPointMultiplier[y] = 1.0
    for y in range(0,100):
      self.dataAxes[y] = "x1y1"
    for y in range(0,100):
      self.dataName[y] = " "
    for y in range(0,100):
      self.dataCircles[y] = None
    for y in range(0,100):
      self.plotType[y] = "linespoints"
    # output file header
    self.plot = "#################################################################################\n"
    self.plot = self.plot + "# Gnuplot\n"
    self.plot = self.plot + "#\n"
    self.plot = self.plot + "#################################################################################\n"
    # batch list
    self.plotBatchList = ""

  def outputPlot(self, fileName):
    self.plotName = fileName

  def title(self, titleIn):
    self.plotTitle = titleIn

  def setDir(self, inDir=None):
    self.dir = inDir

  def axisLabel(self, axis, label):
    if(axis=="x1"):
      self.x1Label = label
    if(axis=="x2"):
      self.x2Label = label
    if(axis=="y1"):
      self.y1Label = label
    if(axis=="y2"):
      self.y2Label = label

  def addCircle(self, circleX, circleY, axis=1):
    if(axis==1):
      line = "set object circle at first "
    if(axis==2):
      line = "set object circle at second "
    line = line + str(circleX) + "," + str(circleY)
    line = line + " radius char 0.7"
    self.dataCircles.append(line)

  def addPlot(self, dataX, dataY, axis=None, dataName="", xMult=1.0, yMult=1.0, plot=None):
    if(axis is None):
      axis = "x1y1"
    if(plot is None):
      plot = 1
    # Set x and y col
    xCol = 2 * self.dataSets
    yCol = 2 * self.dataSets + 1
    # plot multiplier
    self.dataPointMultiplier[xCol] = xMult;
    self.dataPointMultiplier[yCol] = yMult;
    self.dataName[xCol] = dataName
    # Store x data points
    row = 0
    for dataRow in dataX:
      self.cellParameters[xCol][row] = dataRow
      row = row + 1
    self.dataPointCount[xCol] = row
    self.dataAxes[xCol] = axis
    if(row>self.maxRows):
      self.maxRows = row
    row = 0
    # Store y data points
    for dataRow in dataY:
      self.cellParameters[yCol][row] = dataRow
      row = row + 1
    self.dataPointCount[yCol] = row
    # Increment data set counter
    self.dataSets = self.dataSets + 1
    # Last added
    self.lastPlotID = self.dataSets - 1

  def getLastPlotID(self):
    return self.lastPlotID

  def setPlotType(self, plotID, plotType):
    self.plotType[plotID] = plotType

  def makePlot(self):
    ##
    ## Set file names
    ##
    if(self.dir is not None):
      # Make dir
      cmdIn = "mkdir -p "+self.dir
      os.system(cmdIn)

      self.dataFileName = self.dir+"/"+self.plotName+".csv"
      self.gplotFileName = self.dir+"/"+self.plotName+".plot"
      self.plotFileName = self.dir+"/"+self.plotName+".eps"
      self.plotBatchName = self.dir+"/"+self.plotName+".sh"

      self.dataFileName_R = self.plotName+".csv"
      self.gplotFileName_R = self.plotName+".plot"
      self.plotFileName_R = self.plotName+".eps"
      self.plotBatchName_R = self.plotName+".sh"
      self.plotFilePNG_R = self.plotName+".png"
    else:
      self.dataFileName = self.plotName+".csv"
      self.gplotFileName = self.plotName+".plot"
      self.plotFileName = self.plotName+".eps"
      self.plotBatchName = self.plotName+".sh"

      self.dataFileName_R = self.plotName+".csv"
      self.gplotFileName_R = self.plotName+".plot"
      self.plotFileName_R = self.plotName+".eps"
      self.plotBatchName_R = self.plotName+".sh"
      self.plotFilePNG_R = self.plotName+".png"

    ##
    ## Data File
    ##

    wFile = open(self.dataFileName, 'w')
    for row in range (0,self.maxRows):
      fileLine = ""
      i = 0
      for dataSet in range(0,self.dataSets):
        if(i>0):
          fileLine = fileLine + ","
        col = 2 * dataSet
        fileLine = fileLine + str(self.cellParameters[col][row])
        fileLine = fileLine + ","
        fileLine = fileLine + str(self.cellParameters[col+1][row])
        i = i + 1
      wFile.write(str(fileLine)+'\n')
    wFile.close()
    ##
    ## Plot File
    ##

    ## Make plot file
    #self.plot = self.plot + "set terminal postscript eps monochrome enhanced blacktext  size 6.6,3.6\n"
    self.plot = self.plot + "set terminal postscript eps enhanced color size 6.6,3.6\n"
    self.plot = self.plot + "set output \""+self.plotFileName_R+"\"\n"
    self.plot = self.plot + "#\n"
    self.plot = self.plot + "# Set multiple plot layout\n"
    self.plot = self.plot + "#============================================\n"
    self.plot = self.plot + "set multiplot layout 1,1 rowsfirst\n"
    self.plot = self.plot + "#\n"
    self.plot = self.plot + "# Data file\n"
    self.plot = self.plot + "#============================================\n"
    self.plot = self.plot + "set datafile separator \",\"\n"
    ## One for each plot, if multiple plots
    self.plot = self.plot + "#============================================\n"
    self.plot = self.plot + "# Plot 1\n"
    self.plot = self.plot + "#============================================\n"
    self.plot = self.plot + "# Title \n"
    self.plot = self.plot + "set title \""+self.plotTitle+"\"\n"
    self.plot = self.plot + "# Grid settings \n"
    self.plot = self.plot + "#set grid xtics lc rgb \"#CCCCCC\" lw 0.2 lt 1 \n"
    self.plot = self.plot + "#set grid ytics lc rgb \"#CCCCCC\" lw 0.2 lt 1 \n"
    self.plot = self.plot + "# Key settings \n"
    self.plot = self.plot + "set key box opaque \n"
    self.plot = self.plot + "set border back \n"
    self.plot = self.plot + "# Axis \n"
    if(self.x1Label is not None):
      self.plot = self.plot + "set xlabel \""+self.x1Label+"\"\n"
    if(self.x2Label is not None):
      self.plot = self.plot + "set x2label \""+self.x2Label+"\"\n"
    if(self.y1Label is not None):
      self.plot = self.plot + "set ylabel \""+self.y1Label+"\"\n"
    if(self.y2Label is not None):
      self.plot = self.plot + "set y2label \""+self.y2Label+"\"\n"
    self.plot = self.plot + "#set xtics 10\n"
    self.plot = self.plot + "set ytics nomirror tc lt 1\n"
    self.plot = self.plot + "set y2tics nomirror tc lt 1\n"
    self.plot = self.plot + "# Circles \n"
    for circles in self.dataCircles:
      if circles is not None:
        self.plot = self.plot + circles+"\n"
    self.plot = self.plot + "# Plot \n"
    self.plot = self.plot + "plot \\\n"
    y = 0
    for dataSet in range(0,self.dataSets):
      col = 2 * dataSet
      self.plot = self.plot + "'"+self.dataFileName_R+"' "
      self.plot = self.plot + "using "
      self.plot = self.plot + "("+str(self.dataPointMultiplier[col])+" * $"+str(col+1)+")"
      self.plot = self.plot + ":"
      self.plot = self.plot + "("+str(self.dataPointMultiplier[col+1])+" * $"+str(col+2)+") "
      self.plot = self.plot + " title '"+self.dataName[col]+"' "
      self.plot = self.plot + " with "+str(self.plotType[y])+" axes "
      self.plot = self.plot + self.dataAxes[col]
      y = y + 1
      if(dataSet<self.dataSets-1):
        self.plot = self.plot + ", \\\n"

    wFile = open(self.gplotFileName, 'w')
    wFile.write(self.plot)
    wFile.close()
    #time.sleep(0.5)
    ## Run
    wFile = open(self.plotBatchName, 'w')
    wFile.write("#!/bin/bash\n")
    wFile.write("gnuplot "+self.gplotFileName_R+"\n")
    wFile.write("convert -density 300 "+self.plotFileName_R+" "+self.plotFilePNG_R+" \n")
    wFile.close()
    #time.sleep(0.2)
    os.system("chmod +x "+self.plotBatchName)
    #time.sleep(0.2)
    if(self.dir is not None):
      os.system("cd "+self.dir+" && ./"+self.plotBatchName_R)
    else:
      os.system("./"+self.plotBatchName)




################################################################################
