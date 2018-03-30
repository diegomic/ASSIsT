# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'exportDialog.ui'
#
# Created: Thu Apr  2 14:36:31 2015
#      by: PyQt4 UI code generator 4.9.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_exportDialog(object):
    def setupUi(self, exportDialog):
        exportDialog.setObjectName(_fromUtf8("exportDialog"))
        exportDialog.setWindowModality(QtCore.Qt.ApplicationModal)
        exportDialog.resize(496, 295)
        self.buttonBox = QtGui.QDialogButtonBox(exportDialog)
        self.buttonBox.setGeometry(QtCore.QRect(150, 260, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.outfolderSelect = QtGui.QPushButton(exportDialog)
        self.outfolderSelect.setGeometry(QtCore.QRect(460, 190, 31, 30))
        self.outfolderSelect.setObjectName(_fromUtf8("outfolderSelect"))
        self.outfolderText = QtGui.QLineEdit(exportDialog)
        self.outfolderText.setGeometry(QtCore.QRect(170, 190, 291, 29))
        self.outfolderText.setObjectName(_fromUtf8("outfolderText"))
        self.label_2 = QtGui.QLabel(exportDialog)
        self.label_2.setGeometry(QtCore.QRect(30, 190, 131, 21))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.outfileText = QtGui.QLineEdit(exportDialog)
        self.outfileText.setGeometry(QtCore.QRect(170, 220, 291, 29))
        self.outfileText.setObjectName(_fromUtf8("outfileText"))
        self.label_3 = QtGui.QLabel(exportDialog)
        self.label_3.setGeometry(QtCore.QRect(30, 220, 131, 21))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.structureCheckBox = QtGui.QCheckBox(exportDialog)
        self.structureCheckBox.setGeometry(QtCore.QRect(340, 120, 168, 22))
        self.structureCheckBox.setObjectName(_fromUtf8("structureCheckBox"))
        self.plinkCheckBox = QtGui.QCheckBox(exportDialog)
        self.plinkCheckBox.setGeometry(QtCore.QRect(340, 90, 168, 22))
        self.plinkCheckBox.setChecked(False)
        self.plinkCheckBox.setObjectName(_fromUtf8("plinkCheckBox"))
        self.hapmapCheckBox = QtGui.QCheckBox(exportDialog)
        self.hapmapCheckBox.setGeometry(QtCore.QRect(180, 130, 168, 22))
        self.hapmapCheckBox.setObjectName(_fromUtf8("hapmapCheckBox"))
        self.sumCheckBox = QtGui.QCheckBox(exportDialog)
        self.sumCheckBox.setGeometry(QtCore.QRect(20, 50, 120, 41))
        self.sumCheckBox.setChecked(True)
        self.sumCheckBox.setObjectName(_fromUtf8("sumCheckBox"))
        self.gtypesCheckBox = QtGui.QCheckBox(exportDialog)
        self.gtypesCheckBox.setGeometry(QtCore.QRect(20, 80, 131, 41))
        self.gtypesCheckBox.setChecked(True)
        self.gtypesCheckBox.setObjectName(_fromUtf8("gtypesCheckBox"))
        self.snpinfoCheckBox = QtGui.QCheckBox(exportDialog)
        self.snpinfoCheckBox.setGeometry(QtCore.QRect(20, 110, 151, 61))
        self.snpinfoCheckBox.setChecked(True)
        self.snpinfoCheckBox.setObjectName(_fromUtf8("snpinfoCheckBox"))
        self.menderrCheckBox = QtGui.QCheckBox(exportDialog)
        self.menderrCheckBox.setGeometry(QtCore.QRect(180, 40, 141, 61))
        self.menderrCheckBox.setChecked(True)
        self.menderrCheckBox.setObjectName(_fromUtf8("menderrCheckBox"))
        self.label = QtGui.QLabel(exportDialog)
        self.label.setGeometry(QtCore.QRect(21, 21, 301, 21))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setWordWrap(False)
        self.label.setObjectName(_fromUtf8("label"))
        self.locCheckBox = QtGui.QCheckBox(exportDialog)
        self.locCheckBox.setGeometry(QtCore.QRect(180, 90, 120, 41))
        self.locCheckBox.setChecked(True)
        self.locCheckBox.setObjectName(_fromUtf8("locCheckBox"))
        self.flexQTLCheckBox = QtGui.QCheckBox(exportDialog)
        self.flexQTLCheckBox.setGeometry(QtCore.QRect(340, 60, 151, 22))
        self.flexQTLCheckBox.setObjectName(_fromUtf8("flexQTLCheckBox"))

        self.retranslateUi(exportDialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), exportDialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), exportDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(exportDialog)

    def retranslateUi(self, exportDialog):
        exportDialog.setWindowTitle(QtGui.QApplication.translate("exportDialog", "Export Results", None, QtGui.QApplication.UnicodeUTF8))
        self.outfolderSelect.setToolTip(QtGui.QApplication.translate("exportDialog", "Select the Full Data Table input file", None, QtGui.QApplication.UnicodeUTF8))
        self.outfolderSelect.setText(QtGui.QApplication.translate("exportDialog", "...", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("exportDialog", "Output Folder", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("exportDialog", "Output File Prefix", None, QtGui.QApplication.UnicodeUTF8))
        self.structureCheckBox.setText(QtGui.QApplication.translate("exportDialog", "STRUCTURE", None, QtGui.QApplication.UnicodeUTF8))
        self.plinkCheckBox.setText(QtGui.QApplication.translate("exportDialog", "PLINK (.ped, .map)", None, QtGui.QApplication.UnicodeUTF8))
        self.hapmapCheckBox.setText(QtGui.QApplication.translate("exportDialog", "HapMap", None, QtGui.QApplication.UnicodeUTF8))
        self.sumCheckBox.setText(QtGui.QApplication.translate("exportDialog", "Summary", None, QtGui.QApplication.UnicodeUTF8))
        self.gtypesCheckBox.setText(QtGui.QApplication.translate("exportDialog", "Custom gtypes", None, QtGui.QApplication.UnicodeUTF8))
        self.snpinfoCheckBox.setText(QtGui.QApplication.translate("exportDialog", "Custom SNP \n"
"information table", None, QtGui.QApplication.UnicodeUTF8))
        self.menderrCheckBox.setText(QtGui.QApplication.translate("exportDialog", "Custom Mendel\n"
"error report", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("exportDialog", "Export Results as", None, QtGui.QApplication.UnicodeUTF8))
        self.locCheckBox.setText(QtGui.QApplication.translate("exportDialog", "JoinMap (.loc)", None, QtGui.QApplication.UnicodeUTF8))
        self.flexQTLCheckBox.setText(QtGui.QApplication.translate("exportDialog", "FQ_DataPrepper", None, QtGui.QApplication.UnicodeUTF8))

