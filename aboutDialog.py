# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'aboutDialog.ui'
#
# Created: Tue Sep  1 10:47:57 2015
#      by: PyQt4 UI code generator 4.11.3
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

class Ui_aboutDialog(object):
    def setupUi(self, aboutDialog):
        aboutDialog.setObjectName(_fromUtf8("aboutDialog"))
        aboutDialog.setWindowModality(QtCore.Qt.ApplicationModal)
        aboutDialog.resize(653, 596)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(aboutDialog.sizePolicy().hasHeightForWidth())
        aboutDialog.setSizePolicy(sizePolicy)
        aboutDialog.setMaximumSize(QtCore.QSize(726, 700))
        self.okButton = QtGui.QPushButton(aboutDialog)
        self.okButton.setGeometry(QtCore.QRect(510, 550, 122, 30))
        self.okButton.setObjectName(_fromUtf8("okButton"))
        self.label = QtGui.QLabel(aboutDialog)
        self.label.setGeometry(QtCore.QRect(130, 10, 501, 21))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName(_fromUtf8("label"))
        self.image = QtGui.QLabel(aboutDialog)
        self.image.setGeometry(QtCore.QRect(10, 10, 131, 151))
        self.image.setText(_fromUtf8(""))
        self.image.setPixmap(QtGui.QPixmap(_fromUtf8("assist.png")))
        self.image.setScaledContents(True)
        self.image.setObjectName(_fromUtf8("image"))
        self.label_7 = QtGui.QLabel(aboutDialog)
        self.label_7.setGeometry(QtCore.QRect(250, 40, 101, 21))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.label_8 = QtGui.QLabel(aboutDialog)
        self.label_8.setGeometry(QtCore.QRect(150, 110, 471, 51))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.label_9 = QtGui.QLabel(aboutDialog)
        self.label_9.setGeometry(QtCore.QRect(51, 171, 542, 21))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.label_10 = QtGui.QLabel(aboutDialog)
        self.label_10.setGeometry(QtCore.QRect(70, 190, 565, 21))
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.label_11 = QtGui.QLabel(aboutDialog)
        self.label_11.setGeometry(QtCore.QRect(70, 210, 461, 21))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.label_12 = QtGui.QLabel(aboutDialog)
        self.label_12.setGeometry(QtCore.QRect(70, 230, 561, 21))
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.label_13 = QtGui.QLabel(aboutDialog)
        self.label_13.setGeometry(QtCore.QRect(50, 300, 571, 221))
        self.label_13.setTextFormat(QtCore.Qt.RichText)
        self.label_13.setOpenExternalLinks(True)
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.label_14 = QtGui.QLabel(aboutDialog)
        self.label_14.setGeometry(QtCore.QRect(70, 250, 561, 21))
        self.label_14.setObjectName(_fromUtf8("label_14"))

        self.retranslateUi(aboutDialog)
        QtCore.QObject.connect(self.okButton, QtCore.SIGNAL(_fromUtf8("clicked()")), aboutDialog.accept)
        QtCore.QMetaObject.connectSlotsByName(aboutDialog)

    def retranslateUi(self, aboutDialog):
        aboutDialog.setWindowTitle(_translate("aboutDialog", "Credits", None))
        self.okButton.setText(_translate("aboutDialog", "OK", None))
        self.label.setText(_translate("aboutDialog", "ASSIsT:  An Automatic SNP ScorIng Tool in outcross populations", None))
        self.label_7.setText(_translate("aboutDialog", "Version 1.02", None))
        self.label_8.setText(_translate("aboutDialog", "ASSIsT is a tool for efficient filtering and scoring of SNP markers\n"
" from Illumina Infinium/BeadExpress Assays.", None))
        self.label_9.setText(_translate("aboutDialog", "ASSIsT should be cited as Di Guardo and Micheletti et al. 2015 and referenced as:", None))
        self.label_10.setText(_translate("aboutDialog", "Di Guardo M, Micheletti D, Bianco L, Koehorst-van Putten HJJ, Longhi S, Costa F,", None))
        self.label_11.setText(_translate("aboutDialog", "Aranzana MJ,  Velasco R, Arús P, Troggio M, van de Weg EW (2015) ", None))
        self.label_12.setText(_translate("aboutDialog", "ASSIsT: An Automatic SNP ScorIng Tool for in- and out-breeding species.", None))
        self.label_13.setText(_translate("aboutDialog", "<html><head/><body><p>ASSIsT was developed with financial support from the Commission of European <br/>Communities, Seventh Framework Program, Project ‘‘FruitBreedomics: Integrated <br/>approach for increasing breeding efficiency in fruit tree crops’’ <br/>(Grant #FP7-265582; <a href=\"http://www.fruitbreedomics.com\"><span style=\" text-decoration: underline; color:#0000ff;\">www.fruitbreedomics.com</span></a>). <br/>The views and approaches underlying this work are the sole responsibility of the <br/>authors and do not necessary reflect the views of the European Commission.</p><p>For questions contacts Eric Van De Weg at <a href=\"mailto:eric.vandeweg@wur.nl\"><span style=\" text-decoration: underline; color:#0000ff;\">eric.vandeweg@wur.nl</span></a></p><p><span style=\" font-family:\'arial,sans-serif\'; color:#222222;\"><br/>©Fondazione Edmund Mach (</span><a href=\"http://www.fmach.it\"><span style=\" text-decoration: underline; color:#0000ff;\">www.fmach.it</span></a>) - San Michele all\'Adige (TN) - Italy</p></body></html>", None))
        self.label_14.setText(_translate("aboutDialog", "Bioinformatics, DOI: 10.1093/bioinformatics/btv446", None))

