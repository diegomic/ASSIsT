# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mainWindow.ui'
#
# Created: Fri Sep 19 08:16:44 2014
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.setWindowModality(QtCore.Qt.NonModal)
        MainWindow.resize(1324, 813)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        MainWindow.setBaseSize(QtCore.QSize(800, 600))
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.groupBox = QtGui.QGroupBox(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setMinimumSize(QtCore.QSize(320, 421))
        self.groupBox.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.groupBox.setFont(font)
        self.groupBox.setAlignment(QtCore.Qt.AlignCenter)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayout_2 = QtGui.QGridLayout(self.groupBox)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.inFilesLabel = QtGui.QLabel(self.groupBox)
        self.inFilesLabel.setEnabled(True)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.inFilesLabel.setFont(font)
        self.inFilesLabel.setTextFormat(QtCore.Qt.RichText)
        self.inFilesLabel.setObjectName(_fromUtf8("inFilesLabel"))
        self.gridLayout_2.addWidget(self.inFilesLabel, 0, 0, 1, 3)
        self.inFilesButton = QtGui.QPushButton(self.groupBox)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(True)
        font.setWeight(75)
        self.inFilesButton.setFont(font)
        self.inFilesButton.setObjectName(_fromUtf8("inFilesButton"))
        self.gridLayout_2.addWidget(self.inFilesButton, 0, 3, 1, 2)
        self.paramLabel = QtGui.QLabel(self.groupBox)
        self.paramLabel.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.paramLabel.setFont(font)
        self.paramLabel.setTextFormat(QtCore.Qt.PlainText)
        self.paramLabel.setObjectName(_fromUtf8("paramLabel"))
        self.gridLayout_2.addWidget(self.paramLabel, 1, 0, 1, 3)
        self.paramButton = QtGui.QPushButton(self.groupBox)
        self.paramButton.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(True)
        font.setWeight(75)
        self.paramButton.setFont(font)
        self.paramButton.setObjectName(_fromUtf8("paramButton"))
        self.gridLayout_2.addWidget(self.paramButton, 1, 3, 1, 2)
        self.runLabel = QtGui.QLabel(self.groupBox)
        self.runLabel.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.runLabel.setFont(font)
        self.runLabel.setTextFormat(QtCore.Qt.PlainText)
        self.runLabel.setObjectName(_fromUtf8("runLabel"))
        self.gridLayout_2.addWidget(self.runLabel, 2, 0, 1, 3)
        self.runButton = QtGui.QPushButton(self.groupBox)
        self.runButton.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(True)
        font.setWeight(75)
        self.runButton.setFont(font)
        self.runButton.setObjectName(_fromUtf8("runButton"))
        self.gridLayout_2.addWidget(self.runButton, 2, 3, 1, 2)
        self.exportLabel = QtGui.QLabel(self.groupBox)
        self.exportLabel.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.exportLabel.setFont(font)
        self.exportLabel.setTextFormat(QtCore.Qt.PlainText)
        self.exportLabel.setObjectName(_fromUtf8("exportLabel"))
        self.gridLayout_2.addWidget(self.exportLabel, 3, 0, 1, 3)
        self.exportButton = QtGui.QPushButton(self.groupBox)
        self.exportButton.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(True)
        font.setWeight(75)
        self.exportButton.setFont(font)
        self.exportButton.setObjectName(_fromUtf8("exportButton"))
        self.gridLayout_2.addWidget(self.exportButton, 3, 3, 1, 2)
        self.label = QtGui.QLabel(self.groupBox)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Ubuntu"))
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout_2.addWidget(self.label, 4, 0, 1, 4)
        self.consoleTextEdit = QtGui.QTextEdit(self.groupBox)
        font = QtGui.QFont()
        font.setPointSize(7)
        font.setBold(False)
        font.setWeight(50)
        self.consoleTextEdit.setFont(font)
        self.consoleTextEdit.setFocusPolicy(QtCore.Qt.NoFocus)
        self.consoleTextEdit.setAcceptDrops(False)
        self.consoleTextEdit.setAutoFormatting(QtGui.QTextEdit.AutoNone)
        self.consoleTextEdit.setLineWrapMode(QtGui.QTextEdit.NoWrap)
        self.consoleTextEdit.setLineWrapColumnOrWidth(0)
        self.consoleTextEdit.setAcceptRichText(False)
        self.consoleTextEdit.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        self.consoleTextEdit.setObjectName(_fromUtf8("consoleTextEdit"))
        self.gridLayout_2.addWidget(self.consoleTextEdit, 5, 0, 1, 5)
        self.progressLabel = QtGui.QLabel(self.groupBox)
        self.progressLabel.setObjectName(_fromUtf8("progressLabel"))
        self.gridLayout_2.addWidget(self.progressLabel, 6, 0, 1, 2)
        self.progressBar = QtGui.QProgressBar(self.groupBox)
        self.progressBar.setEnabled(True)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setTextVisible(True)
        self.progressBar.setInvertedAppearance(False)
        self.progressBar.setObjectName(_fromUtf8("progressBar"))
        self.gridLayout_2.addWidget(self.progressBar, 6, 2, 1, 3)
        self.progressLabel_2 = QtGui.QLabel(self.groupBox)
        self.progressLabel_2.setObjectName(_fromUtf8("progressLabel_2"))
        self.gridLayout_2.addWidget(self.progressLabel_2, 7, 0, 1, 1)
        self.snpidText = QtGui.QLineEdit(self.groupBox)
        self.snpidText.setObjectName(_fromUtf8("snpidText"))
        self.gridLayout_2.addWidget(self.snpidText, 7, 1, 1, 3)
        self.findButton = QtGui.QPushButton(self.groupBox)
        self.findButton.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(True)
        font.setWeight(75)
        self.findButton.setFont(font)
        self.findButton.setObjectName(_fromUtf8("findButton"))
        self.gridLayout_2.addWidget(self.findButton, 7, 4, 1, 1)
        self.gridLayout.addWidget(self.groupBox, 0, 0, 1, 1)
        self.groupBox_4 = QtGui.QGroupBox(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_4.sizePolicy().hasHeightForWidth())
        self.groupBox_4.setSizePolicy(sizePolicy)
        self.groupBox_4.setMinimumSize(QtCore.QSize(431, 421))
        self.groupBox_4.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.groupBox_4.setFont(font)
        self.groupBox_4.setObjectName(_fromUtf8("groupBox_4"))
        self.inFilesLabel_2 = QtGui.QLabel(self.groupBox_4)
        self.inFilesLabel_2.setEnabled(True)
        self.inFilesLabel_2.setGeometry(QtCore.QRect(10, 30, 401, 421))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.inFilesLabel_2.sizePolicy().hasHeightForWidth())
        self.inFilesLabel_2.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.inFilesLabel_2.setFont(font)
        self.inFilesLabel_2.setText(_fromUtf8(""))
        self.inFilesLabel_2.setTextFormat(QtCore.Qt.RichText)
        self.inFilesLabel_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.inFilesLabel_2.setObjectName(_fromUtf8("inFilesLabel_2"))
        self.gridLayout.addWidget(self.groupBox_4, 0, 2, 1, 1)
        self.groupBox_2 = QtGui.QGroupBox(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_2.sizePolicy().hasHeightForWidth())
        self.groupBox_2.setSizePolicy(sizePolicy)
        self.groupBox_2.setMinimumSize(QtCore.QSize(541, 421))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(10)
        self.groupBox_2.setFont(font)
        self.groupBox_2.setAlignment(QtCore.Qt.AlignCenter)
        self.groupBox_2.setFlat(False)
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.groupBox_2)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.gridLayout_3 = QtGui.QGridLayout()
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.plotWidget = matplotlibWidget(self.groupBox_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plotWidget.sizePolicy().hasHeightForWidth())
        self.plotWidget.setSizePolicy(sizePolicy)
        self.plotWidget.setObjectName(_fromUtf8("plotWidget"))
        self.gridLayout_3.addWidget(self.plotWidget, 0, 0, 1, 1)
        self.horizontalLayout_2.addLayout(self.gridLayout_3)
        self.gridLayout.addWidget(self.groupBox_2, 0, 1, 1, 1)
        self.tableWidget = QtGui.QTableWidget(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tableWidget.sizePolicy().hasHeightForWidth())
        self.tableWidget.setSizePolicy(sizePolicy)
        self.tableWidget.setMinimumSize(QtCore.QSize(1200, 264))
        self.tableWidget.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.tableWidget.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.tableWidget.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.tableWidget.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.tableWidget.setAlternatingRowColors(True)
        self.tableWidget.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.tableWidget.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.tableWidget.setRowCount(0)
        self.tableWidget.setObjectName(_fromUtf8("tableWidget"))
        self.tableWidget.setColumnCount(16)
        item = QtGui.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(2, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(3, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(4, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(5, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(6, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(7, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(8, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(9, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(10, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(11, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(12, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(13, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(14, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(15, item)
        self.tableWidget.horizontalHeader().setVisible(True)
        self.tableWidget.verticalHeader().setVisible(False)
        self.gridLayout.addWidget(self.tableWidget, 1, 0, 1, 3)
        self.horizontalLayout.addLayout(self.gridLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1324, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuHelp = QtGui.QMenu(self.menubar)
        self.menuHelp.setObjectName(_fromUtf8("menuHelp"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionLoadFiles = QtGui.QAction(MainWindow)
        self.actionLoadFiles.setObjectName(_fromUtf8("actionLoadFiles"))
        self.actionSetParameters = QtGui.QAction(MainWindow)
        self.actionSetParameters.setEnabled(False)
        self.actionSetParameters.setObjectName(_fromUtf8("actionSetParameters"))
        self.actionRunAnalysis = QtGui.QAction(MainWindow)
        self.actionRunAnalysis.setEnabled(False)
        self.actionRunAnalysis.setObjectName(_fromUtf8("actionRunAnalysis"))
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionQuit = QtGui.QAction(MainWindow)
        self.actionQuit.setObjectName(_fromUtf8("actionQuit"))
        self.actionReference_Manual = QtGui.QAction(MainWindow)
        self.actionReference_Manual.setObjectName(_fromUtf8("actionReference_Manual"))
        self.actionAbout = QtGui.QAction(MainWindow)
        self.actionAbout.setObjectName(_fromUtf8("actionAbout"))
        self.actionExportResults = QtGui.QAction(MainWindow)
        self.actionExportResults.setEnabled(False)
        self.actionExportResults.setObjectName(_fromUtf8("actionExportResults"))
        self.menuFile.addAction(self.actionLoadFiles)
        self.menuFile.addAction(self.actionSetParameters)
        self.menuFile.addAction(self.actionRunAnalysis)
        self.menuFile.addAction(self.actionExportResults)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionQuit)
        self.menuHelp.addAction(self.actionReference_Manual)
        self.menuHelp.addAction(self.actionAbout)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QObject.connect(self.actionQuit, QtCore.SIGNAL(_fromUtf8("triggered()")), MainWindow.close)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        MainWindow.setTabOrder(self.inFilesButton, self.paramButton)
        MainWindow.setTabOrder(self.paramButton, self.runButton)
        MainWindow.setTabOrder(self.runButton, self.exportButton)
        MainWindow.setTabOrder(self.exportButton, self.tableWidget)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "ASSisT: An automatic SNP Scoring Tool in outcross populations v.1.00", None))
        MainWindow.setToolTip(_translate("MainWindow", "This tool filters SNPs", None))
        self.groupBox.setTitle(_translate("MainWindow", " Console", None))
        self.inFilesLabel.setText(_translate("MainWindow", "<html><head/><body><p>Select Input Files</p></body></html>", None))
        self.inFilesButton.setToolTip(_translate("MainWindow", "Opens Select Input Files dialog", None))
        self.inFilesButton.setText(_translate("MainWindow", "Select", None))
        self.paramLabel.setText(_translate("MainWindow", "Set Parameters", None))
        self.paramButton.setToolTip(_translate("MainWindow", "Opens Set Parameters dialog", None))
        self.paramButton.setText(_translate("MainWindow", "Set", None))
        self.runLabel.setText(_translate("MainWindow", "Run Analysis", None))
        self.runButton.setToolTip(_translate("MainWindow", "Runs the analysis", None))
        self.runButton.setText(_translate("MainWindow", "Run", None))
        self.exportLabel.setText(_translate("MainWindow", "Export Results", None))
        self.exportButton.setToolTip(_translate("MainWindow", "Exports the results in various formats", None))
        self.exportButton.setText(_translate("MainWindow", "Export", None))
        self.label.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:8pt;\">Console Messages:</span></p></body></html>", None))
        self.consoleTextEdit.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Arial\'; font-size:7pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:6pt;\"><br /></p></body></html>", None))
        self.progressLabel.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:9pt;\">Processing progress</span></p><p><span style=\" font-size:9pt;\"><br/></span></p></body></html>", None))
        self.progressLabel_2.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:9pt;\">SNP ID</span></p></body></html>", None))
        self.findButton.setToolTip(_translate("MainWindow", "Opens Select Input Files dialog", None))
        self.findButton.setText(_translate("MainWindow", "Find", None))
        self.groupBox_4.setTitle(_translate("MainWindow", "Results Summary", None))
        self.groupBox_2.setTitle(_translate("MainWindow", "SNP Viewer", None))
        self.plotWidget.setToolTip(_translate("MainWindow", "<html><head/><body><p>SNP clusters will be visualized here</p></body></html>", None))
        self.tableWidget.setSortingEnabled(True)
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Export", None))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "SNP id", None))
        item = self.tableWidget.horizontalHeaderItem(2)
        item.setText(_translate("MainWindow", "Chr", None))
        item = self.tableWidget.horizontalHeaderItem(3)
        item.setText(_translate("MainWindow", "Pos", None))
        item = self.tableWidget.horizontalHeaderItem(4)
        item.setText(_translate("MainWindow", "SNP class", None))
        item = self.tableWidget.horizontalHeaderItem(5)
        item.setText(_translate("MainWindow", "Pop class", None))
        item = self.tableWidget.horizontalHeaderItem(6)
        item.setText(_translate("MainWindow", "Missing", None))
        item = self.tableWidget.horizontalHeaderItem(7)
        item.setText(_translate("MainWindow", "Null", None))
        item = self.tableWidget.horizontalHeaderItem(8)
        item.setText(_translate("MainWindow", "Hom1", None))
        item = self.tableWidget.horizontalHeaderItem(9)
        item.setText(_translate("MainWindow", "Het1", None))
        item = self.tableWidget.horizontalHeaderItem(10)
        item.setText(_translate("MainWindow", "Het2", None))
        item = self.tableWidget.horizontalHeaderItem(11)
        item.setText(_translate("MainWindow", "Hom2", None))
        item = self.tableWidget.horizontalHeaderItem(12)
        item.setText(_translate("MainWindow", "Chi2 pVal", None))
        item = self.tableWidget.horizontalHeaderItem(13)
        item.setText(_translate("MainWindow", "GT Par1", None))
        item = self.tableWidget.horizontalHeaderItem(14)
        item.setText(_translate("MainWindow", "GT Par2", None))
        item = self.tableWidget.horizontalHeaderItem(15)
        item.setText(_translate("MainWindow", "Maf", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.menuHelp.setTitle(_translate("MainWindow", "Help", None))
        self.actionLoadFiles.setText(_translate("MainWindow", "Load Files...", None))
        self.actionSetParameters.setText(_translate("MainWindow", "Set Parameters", None))
        self.actionRunAnalysis.setText(_translate("MainWindow", "Run Analysis", None))
        self.actionExit.setText(_translate("MainWindow", "Exit", None))
        self.actionQuit.setText(_translate("MainWindow", "Exit", None))
        self.actionReference_Manual.setText(_translate("MainWindow", "Reference Manual", None))
        self.actionAbout.setText(_translate("MainWindow", "About", None))
        self.actionExportResults.setText(_translate("MainWindow", "Export Results", None))

from matplotlibWidget import matplotlibWidget
