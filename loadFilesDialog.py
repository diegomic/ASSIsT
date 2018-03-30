# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'loadFilesDialog.ui'
#
# Created: Wed Jan 15 14:51:09 2014
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.setWindowModality(QtCore.Qt.ApplicationModal)
        Dialog.resize(693, 191)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        Dialog.setMaximumSize(QtCore.QSize(693, 228))
        self.buttonBox = QtGui.QDialogButtonBox(Dialog)
        self.buttonBox.setGeometry(QtCore.QRect(330, 150, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.dataTableText = QtGui.QLineEdit(Dialog)
        self.dataTableText.setGeometry(QtCore.QRect(250, 20, 391, 29))
        self.dataTableText.setObjectName(_fromUtf8("dataTableText"))
        self.dnaReportText = QtGui.QLineEdit(Dialog)
        self.dnaReportText.setGeometry(QtCore.QRect(250, 50, 391, 29))
        self.dnaReportText.setObjectName(_fromUtf8("dnaReportText"))
        self.label = QtGui.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(20, 20, 241, 21))
        self.label.setObjectName(_fromUtf8("label"))
        self.dataTableSelect = QtGui.QPushButton(Dialog)
        self.dataTableSelect.setGeometry(QtCore.QRect(650, 20, 31, 30))
        self.dataTableSelect.setObjectName(_fromUtf8("dataTableSelect"))
        self.label_4 = QtGui.QLabel(Dialog)
        self.label_4.setGeometry(QtCore.QRect(20, 80, 241, 21))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setGeometry(QtCore.QRect(20, 50, 241, 21))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.dnaReportSelect = QtGui.QPushButton(Dialog)
        self.dnaReportSelect.setGeometry(QtCore.QRect(650, 50, 31, 30))
        self.dnaReportSelect.setObjectName(_fromUtf8("dnaReportSelect"))
        self.pedigreeText = QtGui.QLineEdit(Dialog)
        self.pedigreeText.setGeometry(QtCore.QRect(250, 80, 391, 29))
        self.pedigreeText.setObjectName(_fromUtf8("pedigreeText"))
        self.pedigreeSelect = QtGui.QPushButton(Dialog)
        self.pedigreeSelect.setGeometry(QtCore.QRect(650, 80, 31, 30))
        self.pedigreeSelect.setObjectName(_fromUtf8("pedigreeSelect"))
        self.label_5 = QtGui.QLabel(Dialog)
        self.label_5.setGeometry(QtCore.QRect(20, 110, 241, 21))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.mapFileSelect = QtGui.QPushButton(Dialog)
        self.mapFileSelect.setGeometry(QtCore.QRect(650, 110, 31, 30))
        self.mapFileSelect.setObjectName(_fromUtf8("mapFileSelect"))
        self.mapFileText = QtGui.QLineEdit(Dialog)
        self.mapFileText.setGeometry(QtCore.QRect(250, 110, 391, 29))
        self.mapFileText.setObjectName(_fromUtf8("mapFileText"))

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), Dialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)
        Dialog.setTabOrder(self.dataTableText, self.dnaReportText)
        Dialog.setTabOrder(self.dnaReportText, self.pedigreeText)
        Dialog.setTabOrder(self.pedigreeText, self.mapFileText)
        Dialog.setTabOrder(self.mapFileText, self.dataTableSelect)
        Dialog.setTabOrder(self.dataTableSelect, self.dnaReportSelect)
        Dialog.setTabOrder(self.dnaReportSelect, self.pedigreeSelect)
        Dialog.setTabOrder(self.pedigreeSelect, self.mapFileSelect)
        Dialog.setTabOrder(self.mapFileSelect, self.buttonBox)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Select Input Files", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("Dialog", "Select Final Report File", None, QtGui.QApplication.UnicodeUTF8))
        self.dataTableSelect.setToolTip(QtGui.QApplication.translate("Dialog", "Select the Full Data Table input file", None, QtGui.QApplication.UnicodeUTF8))
        self.dataTableSelect.setText(QtGui.QApplication.translate("Dialog", "...", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("Dialog", "Select Pedigree File", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("Dialog", "Select DNA Report File", None, QtGui.QApplication.UnicodeUTF8))
        self.dnaReportSelect.setToolTip(QtGui.QApplication.translate("Dialog", "Select the DNA Report  input file", None, QtGui.QApplication.UnicodeUTF8))
        self.dnaReportSelect.setText(QtGui.QApplication.translate("Dialog", "...", None, QtGui.QApplication.UnicodeUTF8))
        self.pedigreeSelect.setToolTip(QtGui.QApplication.translate("Dialog", "Select the DNA Report  input file", None, QtGui.QApplication.UnicodeUTF8))
        self.pedigreeSelect.setText(QtGui.QApplication.translate("Dialog", "...", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("Dialog", "Select Map File", None, QtGui.QApplication.UnicodeUTF8))
        self.mapFileSelect.setToolTip(QtGui.QApplication.translate("Dialog", "Select the DNA Report  input file", None, QtGui.QApplication.UnicodeUTF8))
        self.mapFileSelect.setText(QtGui.QApplication.translate("Dialog", "...", None, QtGui.QApplication.UnicodeUTF8))

