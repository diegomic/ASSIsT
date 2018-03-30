# -*- coding: utf-8 -*-
from PyQt4 import QtGui, QtCore, Qt
from PyQt4.QtCore import Qt
from PyQt4.QtGui import QApplication, QMainWindow, QDialog, QFileDialog, QTextCursor, QMessageBox, QColor
try:
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
except ImportError:
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import os
import sys
import webbrowser
import subprocess
from aboutDialog import Ui_aboutDialog
from exportDialog import Ui_exportDialog
from loadFilesDialog import Ui_Dialog
from mainWin import Ui_MainWindow
from paramDialog import Ui_paramsDialog
import mainProcess as snp_filter_pipeline
import snp_classes

#Import the SNP_filter_pipeline
class inFilesDialog(QDialog, Ui_Dialog):
        def __init__(self):
            super(inFilesDialog, self).__init__(None)
            self.ui = Ui_Dialog()
            self.ui.setupUi(self)
            self.foldName = None
            #Connect buttons.
            self.ui.dataTableSelect.clicked.connect(lambda : self.openFileChooser(self.ui.dataTableText, "TXT files (*.txt)"))
            #self.ui.snpTableSelect.clicked.connect(lambda : self.openFileChooser(self.ui.snpTableText,"TXT files (*.txt)"))
            self.ui.dnaReportSelect.clicked.connect(lambda : self.openFileChooser(self.ui.dnaReportText,"CSV files (*.csv)"))
            self.ui.pedigreeSelect.clicked.connect(lambda : self.openFileChooser(self.ui.pedigreeText,"TXT files (*.txt)"))
            self.ui.mapFileSelect.clicked.connect(lambda : self.openFileChooser(self.ui.mapFileText,"TXT files (*.txt)"))        
            
            
            
        def openFileChooser(self,qtext,filetype):
            if (self.foldName == None):
                self.foldName = QtCore.QDir.homePath()
                
            if (len(qtext.text()) > 0):
                self.foldName = qtext.text()
            str= QFileDialog.getOpenFileName(self, "Choose file", self.foldName, filetype)
            if str:
                qtext.setText(str)
                self.foldName=str
                
        def getDataFiles(self):
            res=dict()
            t=self.ui.dataTableText.text()
            res["finalreport"]=str(t)
            #t=self.ui.snpTableText.text()
            #res["snpTable"]=str(t)
            t=self.ui.dnaReportText.text()
            res["dnareport"]=str(t)
            t=self.ui.pedigreeText.text()
            res["pedigreefile"]=str(t)
            t=self.ui.mapFileText.text()
            res["mapfile"]=str(t)
            return res
        
        def parsePedigreesFile(self):
            pedigrees={}
            pf=open(self.ui.pedigreeText.text(),'r')
            cnt=0 #counter to skip first element (the header)
            errorMsg=""
            try:
                for line in pf:
                    if(cnt ==0):
                        cnt+=1
                    else:
#                         line = line.decode('utf-8')
                        line=line.strip('\n')
                        line=line.strip('\r')
                        els=line.split('\t')
                        if(els[0] not in pedigrees):
                            pedigrees[els[0]]=[els[1],els[2]]
            except:
                raise
                errorMsg="Error: unable to parse pedigree file!"
            return pedigrees,errorMsg    

class paramsDialog(QDialog, Ui_paramsDialog):
        def __init__(self):
            super(paramsDialog, self).__init__(None)
            self.ui = Ui_paramsDialog()
            self.ui.setupUi(self)
            #Connect buttons.
            self.ui.populationCombo.currentIndexChanged.connect(self.toggleView)
            #self.ui.genTrainLine.editingFinished.connect(lambda : self.checkInputData(self.ui.genTrainLine,"float",[0,1],"genTrainScore"))
            #self.ui.gcLine.editingFinished.connect(lambda : self.checkInputData(self.ui.gcLine,"float",[0,1],"10% GC range"))
            self.ui.missingLine.editingFinished.connect(lambda : self.checkInputData(self.ui.missingLine,"float",[0,1],"proportion Missing"))
            self.ui.toleranceLine.editingFinished.connect(lambda : self.checkInputData(self.ui.toleranceLine,"float",[0,1],"call rate tolerance"))
            self.ui.pvalueLine.editingFinished.connect(lambda : self.checkInputData(self.ui.pvalueLine,"float",[0,1],"pValue"))
            self.ui.unexpGenotypeLine.editingFinished.connect(lambda : self.checkInputData(self.ui.unexpGenotypeLine,"float",[0,1],"unexpected Genotype"))
            self.ui.unexpGenotypeLine_snp.editingFinished.connect(lambda : self.checkInputData(self.ui.unexpGenotypeLine_snp,"float",[0,1],"unexpected Genotype"))
            self.ui.rareLine.editingFinished.connect(lambda : self.checkInputData(self.ui.rareLine,"float",[0,1],"frequency rare allele"))
            self.ui.chrLine.editingFinished.connect(lambda : self.checkInputData(self.ui.chrLine,"int",[0,999999],"chromosomes"))
            
    
            
        def checkInputData(self,what,datatype,range,dataLabel):
            val=what.text()
            try :
                if(datatype == "float"):
                    val=float(val)
                else:
                    val=int(val)
                if(range[0] > val or range[1]< val):
                    warn=QMessageBox()
                    warn.setWindowTitle("Input Error: " +dataLabel) 
                    warn.setText("Warning "+dataLabel + " value (" + str(val)+") out of range.\nCorrect range is ["+str(range[0]) + ","+str(range[1])+ "].")
                    warn.exec_()
                    what.setFocus(True)
            except:
                warn=QMessageBox()
                warn.setWindowTitle("Input Error: " +dataLabel)
                warn.setText("Warning "+dataLabel + " value (" + str(val)+") is not numeric.\nPlease provide a value in correct range ["+str(range[0]) + ","+str(range[1])+ "].")
                warn.exec_()
                what.setFocus(True)
                
        def toggleView(self):
            if(self.ui.populationCombo.currentText() == "Germplasm"):
                self.ui.parentsCombo.setEnabled(False)
                self.ui.nullAlleleCombo.setEnabled(False)
                self.ui.unexpGenotypeLine_snp.setEnabled(False)
                self.ui.label_13.setEnabled(False)
                self.ui.rareLine.setEnabled(True)
                self.ui.label_9.setEnabled(True)
            elif (self.ui.populationCombo.currentText() == "F2"):
                self.ui.nullAlleleCombo.setEnabled(False)
                self.ui.unexpGenotypeLine_snp.setEnabled(True)
                self.ui.label_13.setEnabled(True)
                self.ui.rareLine.setEnabled(False)
                self.ui.label_9.setEnabled(False)
            else:
                self.ui.parentsCombo.setEnabled(True)
                self.ui.nullAlleleCombo.setEnabled(True)
                self.ui.unexpGenotypeLine_snp.setEnabled(True)
                self.ui.label_13.setEnabled(True)
                self.ui.rareLine.setEnabled(False)
                self.ui.label_9.setEnabled(False)
                
        def getParams(self):
            res=dict()
            #===================================================================
            # t=self.ui.genTrainLine.text()
            # res["gentrain"]=float(t)
            # t=self.ui.gcLine.text()
            # res["gencall"]=float(t)
            #===================================================================
            t=self.ui.missingLine.text()
            res["missing"]=float(t)
            t=self.ui.toleranceLine.text()
            res["tolerance"]=float(t)
            t=self.ui.pvalueLine.text()
            res["pvalue"]=float(t)
            
                
            t=self.ui.unexpGenotypeLine.text()
            res["thresh_out"]=float(t)
            if self.ui.populationCombo.currentText() == 'Germplasm':
                t=self.ui.rareLine.text()
                res["rare"]=float(t)
            else:
                t=self.ui.unexpGenotypeLine_snp.text()
                res["rare"]=float(t)
            t=self.ui.populationCombo.currentText()
            if t == 'CP (F1)' or t == 'BC':
                t = 'F1'
            res["pop"]=str(t)
            t=self.ui.parentsCombo.currentText()
            try:
                [p1,p2]=str(t).split(' x ')
            except  ValueError:
                [p1,p2] = str('-'), str('-')
            res["p1"]=str(p1)
            res["p2"]=str(p2)
            outcrList=[str(v.text()) for v in self.ui.outcrossList.selectedItems()]
            res["outcross"]=outcrList
            t=self.ui.chrLine.text()
            res["chr"]=int(t)
            t=self.ui.nullAlleleCombo.currentText()
            #TODO: at the end unify recover in snp_filter_pipeline or separate efXeg null in GUI. 
            res["recover-efxeg"]=t
            res["recover-null"]=t
            
            #DEFAULT parameters!!!!
            res["null_R"]=0.3
            return res
        
        
        
class exportDialog(QDialog, Ui_exportDialog):
        def __init__(self,startFolder):
            super(exportDialog, self).__init__(None)
            self.ui = Ui_exportDialog()
            self.ui.setupUi(self)
            #Connect buttons.
            self.ui.outfolderSelect.clicked.connect(lambda : self.openDirChooser(self.ui.outfolderText,""))
            self.ui.outfolderText.editingFinished.connect(lambda : self.checkInputData(self.ui.outfolderText,"outputfolder"))
            self.ui.outfileText.editingFinished.connect(lambda : self.checkInputData(self.ui.outfileText,"outputfile"))
            self.ui.outfolderText.setText(startFolder)
            #Set Focus to enforce specifying outfolder
            self.ui.outfolderText.setFocus(True)


        def checkInputData(self,what,dataLabel):
            val=str(what.text())
            if(len(val) == 0):
                    warn=QMessageBox()
                    warn.setWindowTitle("Input Error: " +dataLabel)
                    warn.setText("Warning "+dataLabel + " cannot be empty. Please specify a value.")
                    warn.exec_()
                    if(dataLabel == "outputfile"):
                        self.ui.outfileText.setFocus(True)
            else:

                if(dataLabel == "outputfolder"):
                    if(not os.path.isdir(val) or (os.path.isdir(val) and not os.access(val,os.W_OK))):
                        try:
                            os.makedirs(val)
                        except:
                            warn=QMessageBox()
                            warn.setWindowTitle("Input Error: " +dataLabel)
                            if(os.path.isdir(val)):
                                warn.setText("Cannot WRITE to output folder "+val + ". Please specify a different one.")
                            else:
                                warn.setText("Cannot CREATE output folder "+val + ". Please specify a different one.")
                            warn.exec_()

                            self.openDirChooser(self.ui.outfolderText,"")

        #def saveResultsToFile(self,loc,plink,structure,flexQTL,hapmap, table):
        def saveResultsToFile(self,tmp_out_type, table):
            '''
            Save data to files'''
            export_states = {} # export_states[snp] = False(0)/True(2)
            num_rows, num_cols = table.rowCount(), table.columnCount()
            for row in range(num_rows):
                try:
                    state = table.item(row, 0).checkState()
                    snp_id = str(table.item(row, 1).text())
                    if state == 2:
                        export_states[snp_id] = True
                    else:
                        export_states[snp_id] = False
                except AttributeError:
                    continue
            #=======================================================================================
            # tmp_out_type = {'joinmap': loc, 'flexqtl': flexQTL,'hapmap': hapmap,
            #             'plink': plink, 'structure': structure}
            #=======================================================================================
            out_type = {}
            for k, v in tmp_out_type.items():
                if v == 2:
                    out_type[k] = True
                else: 
                    out_type[k] = False
            return export_states, out_type    
                      
        def openDirChooser(self,qtext,filetype):
            rootDir=QtCore.QDir.homePath()
            if (len(qtext.text()) > 0):
                rootDir = qtext.text()
            myDial = QFileDialog()
            myDial.setFileMode(2)
            str= myDial.getExistingDirectory(self, "Choose output directory")
            if str:
                qtext.setText(str)
                
                
class aboutDialog(QDialog, Ui_aboutDialog):
        def __init__(self):
            super(aboutDialog, self).__init__(None)
            self.ui = Ui_aboutDialog()
            self.ui.setupUi(self)
            #Connect buttons.


class programWindow(QMainWindow,Ui_MainWindow):
  
        
    def __init__(self):
                    
        super(programWindow, self).__init__(None)

        # Set up the user interface from Designer.
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        
        self.dialog=None
        self.params=None
        self.export=None
        self.about=None
        self.pedigrees={}
        self.paramsArray={}
        self.inFilesArray={}
        self.ui.consoleTextEdit.setText("Welcome to SNP Filter...\nStart Loading Input Files")
        
        
        #Connect buttons.
        self.ui.inFilesButton.clicked.connect(lambda : self.showLoadFiles(self.inFilesArray))
        self.ui.paramButton.clicked.connect(lambda : self.showSetParams(self.paramsArray,self.pedigrees))
        self.ui.runButton.clicked.connect(lambda : self.runProgram(self.inFilesArray,self.paramsArray))
        self.ui.exportButton.clicked.connect(lambda : self.showExport(self.inFilesArray,self.paramsArray))
        self.ui.findButton.clicked.connect(self.findSNP)
        self.ui.snpidText.returnPressed.connect(self.findSNP)
        #Connect actions within menu
        self.ui.actionLoadFiles.triggered.connect(lambda : self.showLoadFiles(self.inFilesArray))
        self.ui.actionSetParameters.triggered.connect(lambda : self.showSetParams(self.paramsArray,self.pedigrees))
        self.ui.actionExportResults.triggered.connect(lambda : self.showExport(self.inFilesArray,self.paramsArray))
        self.ui.actionRunAnalysis.triggered.connect(lambda : self.runProgram(self.inFilesArray,self.paramsArray))
        self.ui.actionAbout.triggered.connect(self.showAbout)
        
        self.ui.actionReference_Manual.triggered.connect(self.showManual)
        
        self.memmap_file_order, self.memmap_file = None, None
        self.offsprings, self.ind_order = None, None
        #Connect table selection change
        self.ui.tableWidget.itemSelectionChanged.connect(lambda : self.plotSNP(self.inFilesArray,
                                                                               self.paramsArray, 
                                                                               self.ui.tableWidget.item(self.ui.tableWidget.currentRow(),1).text(),
                                                                               self.ui.tableWidget.currentRow()))
        self.ui.progressBar.setVisible(False)
        self.ui.progressLabel.setVisible(False)
        mpl_toolbar = NavigationToolbar(self.ui.plotWidget.canvas, self)
        self.ui.gridLayout_3.addWidget(mpl_toolbar,1,0)
        self.ui.gridLayout.setAlignment(mpl_toolbar,QtCore.Qt.AlignCenter)

        
    #Methods.
    #load input files dialog
    def showLoadFiles(self,inputsArray):
        pedigrees={}
        if self.dialog is None:
            self.dialog= inFilesDialog()
        self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+"\nEditing Input files...")
        res=self.dialog.exec_()
        self.inFilesArray=self.dialog.getDataFiles()
        errorMsg=""
        outRes=self.inFilesArray
        for el in outRes:
            if(len(outRes[el])== 0 and (el != 'mapfile')):
                errorMsg+=el +" file cannot be empty\n"
                
            else:
                if(el=='mapfile' and len(outRes[el])==0):
                    continue
                if( not os.path.exists(outRes[el])):
                    errorMsg+=el +" file \'"+outRes[el]+"' does not exist!\n"
        
        if(errorMsg== ""):
            self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+" done.\nFiles looking OK.")
            #LOAD the pedigree file into a suitable data structure
            [pedigrees,msg]=self.dialog.parsePedigreesFile()
            if(len(msg)==0):
                self.ui.paramButton.setEnabled(True)
                self.ui.actionSetParameters.setEnabled(True)
                self.ui.paramLabel.setEnabled(True)
                self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+"\nCan now proceed to setting parameters.")
                self.ui.runButton.setEnabled(False)
                self.ui.runLabel.setEnabled(False)
                self.ui.exportButton.setEnabled(False)
                self.ui.exportLabel.setEnabled(False)
                self.ui.actionRunAnalysis.setEnabled(False)
            else:
                self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+"\n"+msg) 
                self.ui.paramButton.setEnabled(False)
                self.ui.paramLabel.setEnabled(False)
                self.ui.inFilesLabel.setEnabled(False)
                self.ui.runButton.setEnabled(False)
                self.ui.runLabel.setEnabled(False)
                self.ui.exportButton.setEnabled(False)
                self.ui.exportLabel.setEnabled(False)
                self.ui.actionRunAnalysis.setEnabled(False)
                self.ui.actionExport.setEnabled(False)   
        else:
            self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+" done.\n")
            cc = self.ui.consoleTextEdit.textColor()
            self.ui.consoleTextEdit.setTextColor(QColor("red"))
            self.ui.consoleTextEdit.append("ERRORS FOUND.\n"+errorMsg)
            self.ui.consoleTextEdit.setTextColor(cc)
            self.ui.consoleTextEdit.append("")
            self.ui.paramButton.setEnabled(False)
            self.ui.actionSetParameters.setEnabled(False)
            self.ui.runButton.setEnabled(False)
            self.ui.runLabel.setEnabled(False)
            self.ui.actionRunAnalysis.setEnabled(False)
            self.ui.exportButton.setEnabled(False)
            self.ui.exportLabel.setEnabled(False)
            self.ui.actionExportResults.setEnabled(False)
            self.ui.actionRunAnalysis.setEnabled(False)
        self.ui.consoleTextEdit.moveCursor(QTextCursor.End)    
        self.pedigrees=pedigrees


        
    #load set parameters dialog
    def showSetParams(self,paramsArray,pedigrees):
        if self.params is None:
            self.params= paramsDialog()
        labels= set()
        samples=set()
        
        for el in pedigrees:
            samples.add(el)
            l=pedigrees[el][0]+" x " + pedigrees[el][1]
            if(l != " x " and l != "x"  and l != "- x -"and l != "0 x 0" ):
                labels.add(pedigrees[el][0]+" x " + pedigrees[el][1])
        sortLab=list(labels)
        sortLab.sort()
        self.params.ui.parentsCombo.clear()
        for el in sortLab:
            self.params.ui.parentsCombo.addItem(el)
        
        sampSort=list(samples)
        sampSort.sort()
        self.params.ui.outcrossList.clear()
        for el in sampSort:
            self.params.ui.outcrossList.addItem(el)
        self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+"\nEditing Parameters...")    
        res=self.params.exec_()
        self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+"done.\n")
        self.paramsArray=self.params.getParams()
        self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+"Parameters look OK.\nAll set to RUN the analysis.")
        self.ui.runButton.setEnabled(True)
        self.ui.runLabel.setEnabled(True)
        self.ui.actionRunAnalysis.setEnabled(True)
        self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
        
    #run program
    def runProgram(self,inputsArray,paramsArray):
        self.ui.progressBar.setValue(0)
        self.ui.progressBar.setVisible(True)
        self.ui.progressLabel.setVisible(True)
        self.ui.findButton.setEnabled(False)
        #Clear Table and Summary when Run start
        self.ui.inFilesLabel_2.clear()
        self.ui.tableWidget.clearContents()
        self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+"\nRUNNING...")
        self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
        self.ui.tableWidget.setRowCount(0)
        self.ui.tableWidget.setSortingEnabled(False)

        paramsArray.update(inputsArray)
        
        #Let's run the pipeline
        #SNPs number
        snp_number = snp_filter_pipeline.RunAnalysis(paramsArray,self.ui).final_report_header()
        (self.output_summary, self.snptable, self.snptable_index, self.memmap_file_order, self.memmap_file,
          self.recovered_gtypes, self.ui.progressBar, self.offsprings, self.ind_order, self.unexpected_gt) = \
                        snp_filter_pipeline.RunAnalysis(paramsArray,self.ui).display_output(self.ui.progressBar, self.ui.inFilesLabel_2)
        #populate Summary
        
        #self.ui.inFilesLabel_2.setText(self.output_summary.replace('\n','<br />'))
        # Populate Table with SNPs
        self.ui.tableWidget.setColumnCount(15)
        self.ui.tableWidget.setRowCount(snp_number)

        # Cycle over snptable and fill the tableWidget
        if paramsArray['pop'] == 'Germplasm':
            exported_class = ['Robust',                       # 'Robust'
                              'OneHomozygRare_HWE',      # 'OneHomoRare_HWE'
                              'OneHomozygRare_NotHWE',  # 'OneHomoRare_NotHWE'
                              'DistortedAndUnexpSegreg']                    # 'DistortedAndUnexpSegreg'
           
            non_exported_class = ['Monomorphic',                    # 'FalseSNP'
                                  'Failed',                       # 'Failed'
                                  'ShiftedHomo',                  # 'ShiftedHomo'
                                  'NullAllele-Failed',            # 'NullAllele'
                                  'NotInformative']#             # 'NonInformative'
                              
        else:
            exported_class =['Robust',                       # 'Robust'
                             'AB_2_sub-clusters',           # 'Recovered_2Het'
                             'Null_2_Cluster',  # '2ClusterNull'
                             'Null_4_Cluster'] 
            
            non_exported_class = ['Monomorphic',                    # 'FalseSNP'
                                  'Failed',                       # 'Failed'
                                  'DistortedAndUnexpSegreg',                    # 'DistortedAndUnexpSegreg'
                                  'NullAllele-Failed',            # 'NullAllele'
                                  'ShiftedHomo',                  # 'ShiftedHomo'
                                  'NotInformative']              # 'NonInformative'

        for ln, snp in enumerate(self.snptable):
            #chkBoxItem = QtGui.QTableWidgetItem()
            chkBoxItem=CheckBoxItem()
            chkBoxItem.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
            if snp[3] in exported_class:
                chkBoxItem.setCheckState(QtCore.Qt.Checked)
                self.ui.tableWidget.setItem(ln, 0, chkBoxItem)
            else:
                chkBoxItem.setCheckState(QtCore.Qt.Unchecked)
                self.ui.tableWidget.setItem(ln, 0, chkBoxItem)
            for c_num in range(1,15):
                #item = QtGui.QTableWidgetItem(QtCore.QString("%1").arg(snp[c_num-1]))
                item=QtGui.QTableWidgetItem()
                try:
                    v=float(snp[c_num-1])
                    item.setData(Qt.DisplayRole,v)

                except:

                    try:
                        v=int(snp[c_num-1])
                        item.setData(Qt.DisplayRole,v)

                    except:
                        item.setData(Qt.DisplayRole,snp[c_num-1])

                item.setData(Qt.DisplayRole,snp[c_num-1])
                self.ui.tableWidget.setItem(ln, c_num, item)

#            item=QtGui.QTableWidgetItem()
#            item.setData(Qt.DisplayRole,ln)
#            self.ui.tableWidget.setItem(ln, 16,item)

        self.ui.tableWidget.setSortingEnabled(True)
        self.ui.exportButton.setEnabled(True)
        self.ui.exportLabel.setEnabled(True)
        self.ui.actionExportResults.setEnabled(True)    
        self.ui.exportLabel.setEnabled(True)
        
        self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+"done.")
        self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
        self.ui.findButton.setEnabled(True)
        self.ui.progressBar.setVisible(False)
        self.ui.progressLabel.setVisible(False)


            
    #load export dialog
    def showExport(self,inputsArray,paramsArray):
        paramsArray.update(inputsArray)
        if self.export is None:
            curFiles=self.inFilesArray
            startFold=os.path.dirname(curFiles["finalreport"])
            self.export= exportDialog(startFold)
        self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+"\nEditing Export Parameters...")    
        res=self.export.exec_()
        self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+"done.")
        self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
        if(res == 1):
            fileExpCheckState = {'joinmap':self.export.ui.locCheckBox.checkState(),
                                 'plink': self.export.ui.plinkCheckBox.checkState(),
                                 'structure':self.export.ui.structureCheckBox.checkState(),
                                 'flexqtl': self.export.ui.flexQTLCheckBox.checkState(), 
                                 'hapmap': self.export.ui.hapmapCheckBox.checkState(),
                                 'snpinfo': self.export.ui.snpinfoCheckBox.checkState(),
                                 'sum': self.export.ui.sumCheckBox.checkState(),
                                 'mendel': self.export.ui.menderrCheckBox.checkState(),
                                 'gtypes': self.export.ui.gtypesCheckBox.checkState()}
            export_states, out_type = self.export.saveResultsToFile(fileExpCheckState,
                                                                    self.ui.tableWidget)
            # export_states, out_type = self.export.saveResultsToFile(self.export.ui.locCheckBox.checkState(),
                                                      #-- self.export.ui.plinkCheckBox.checkState(),
                                                      # self.export.ui.structureCheckBox.checkState(),
                                                      # self.export.ui.flexQTLCheckBox.checkState(),
                                                      #- self.export.ui.hapmapCheckBox.checkState(),
                                                      #------------------------ self.ui.tableWidget)
#         print str(self.export.ui.outfolderText.text()), str(self.export.ui.outfileText.text())
            file_exp = snp_filter_pipeline.FileExporter(paramsArray, self.memmap_file_order, 
                                                        self.memmap_file,
                                                        self.recovered_gtypes, self.unexpected_gt,
                                                        self.offsprings, self.ind_order,
                                                        str(self.export.ui.outfolderText.text()),
                                                        str(self.export.ui.outfileText.text()),
                                                        self.ui)


            file_exp.export_output(out_type, self.snptable, self.snptable_index, export_states, 
                                   self.output_summary)
            self.ui.consoleTextEdit.setText(self.ui.consoleTextEdit.toPlainText()+
                                            "\nResults SUCCESSFULLY written to " + 
                                            self.export.ui.outfolderText.text())
            self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
        

    #load about dialog
    def showAbout(self):
        if self.about is None:
            self.about= aboutDialog()
        res=self.about.exec_()
        
    def plotSNP(self, paramsArray, inputsArray, SNPname, tableRowN):
        #self.ui.tableWidget.setItemSelected(self.ui.tableWidget.item(tableRowN,0),True)
        paramsArray.update(inputsArray)
        #print "SNP is:"+str(SNPname)
        main = snp_filter_pipeline.GslikePlot(paramsArray, self.memmap_file_order, self.memmap_file, self.offsprings, self.ind_order)
        self.ui.plotWidget.canvas = main.make_plot(SNPname, self.ui.plotWidget.canvas)


        #self.ui.plotWidget.pl

    def findSNP(self):
        snpid=self.ui.snpidText.text()
        found=self.ui.tableWidget.findItems(snpid,QtCore.Qt.MatchContains) #2 == MatchContains
        #print found[0].setSelected(True)
        if(len(found)>0):
            row=found[0].row()
            self.ui.tableWidget.setCurrentItem(found[0])
            self.ui.tableWidget.scrollToItem(found[0])
            self.ui.tableWidget.setFocus()
            
    def showManual(self):
        filepath = os.path.join(os.curdir, 'docs','ASSIsT_Reference_Manual.pdf')
        if sys.platform.startswith('darwin'):
            subprocess.call(('open', filepath))
        elif os.name == 'nt':
            os.system("start "+filepath)
        elif os.name == 'posix':
            subprocess.call(('xdg-open', filepath))
        #webbrowser.open(os.path.join(os.curdir, 'docs','_build','html','index.html'))#'./docs/_build/html/index.html')
        
class CheckBoxItem(QtGui.QTableWidgetItem):
    def __lt__(self, other):
        if ( isinstance(other, QtGui.QTableWidgetItem) ):
            my_value=self.checkState()
            other_value=other.checkState()
            return my_value < other_value

        return super(CheckBoxItem, self).__lt__(other)
        
if __name__ == "__main__":
    #visualize the GUI        
    app = QApplication(sys.argv)
    ui = programWindow()
    app.setWindowIcon(QtGui.QIcon('assist.ico'))
    ui.setWindowIcon(QtGui.QIcon('assist.ico'))
    ui.show()
    app.exec_()