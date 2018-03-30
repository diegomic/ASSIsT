# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'errorsDialog.ui'
#
# Created: Fri Feb 21 11:32:48 2014
#      by: PyQt4 UI code generator 4.9.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_errorDialog(object):
    def setupUi(self, errorDialog):
        errorDialog.setObjectName(_fromUtf8("errorDialog"))
        errorDialog.resize(432, 234)
        self.buttonBox = QtGui.QDialogButtonBox(errorDialog)
        self.buttonBox.setGeometry(QtCore.QRect(80, 200, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.errorConsole = QtGui.QTextEdit(errorDialog)
        self.errorConsole.setGeometry(QtCore.QRect(10, 10, 411, 171))
        self.errorConsole.setObjectName(_fromUtf8("errorConsole"))

        self.retranslateUi(errorDialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), errorDialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), errorDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(errorDialog)

    def retranslateUi(self, errorDialog):
        errorDialog.setWindowTitle(QtGui.QApplication.translate("errorDialog", "Error Message", None, QtGui.QApplication.UnicodeUTF8))

