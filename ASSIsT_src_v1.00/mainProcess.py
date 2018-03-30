# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 15:24:42 2014
"""
__author__ = 'Diego Micheletti'
__version__ = '0.1'

    
from PyQt4.QtGui import QTextCursor, QColor
#from matplotlib import pyplot as plt
from operator import itemgetter
#from scipy.cluster.vq import vq, kmeans2, whiten
from scipy.stats import chisquare as chisq
from scipy.stats import shapiro
#from scipy.stats import ttest_ind
from tempfile import mkdtemp
from textwrap import dedent
import argparse
import csv
import os
import random
import re
import snp_classes
import string
import sys
import warnings

import numpy as np
#from mainWin import Ui_MainWindow
global jmgtypes
jmgtypes = {}
#Plot in QT canvas
#from PyQt4 import QtGui
#from PyQt4.QtGui import QApplication, QMainWindow, QDialog, QFileDialog, QTextCursor, QMessageBox
# from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
# from errorsDialog import Ui_errorDialog
# Hide warnings
warnings.filterwarnings("ignore")


class OutcrossAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))


def get_args():
    '''
    Parse the arguments from the command line.
    
    Returns:
        ::Args: NameSpace containing all the arguments and options to perform the analysis
    '''
    parser = argparse.ArgumentParser(formatter_class =
                        argparse.RawTextHelpFormatter,
                    description=dedent('''\
                      **WARNING:**\n The use of this script is deprecated. To run ASSIsT use the
                      ing provided GUI that can be accessed by running the ASSIsT_1_00.py. 
                      Script to filter useful markers from GenomeGtudio outputs.
                      The script take Full Data Table, SNP Table, and DNA Report
                      from genome studio. Additionally a Pedigree and a map file
                      have to be provided. Several population type are accepted:
                      F1, CP, BC, F2, and unrelated germplasm.
                      ----------------------------------------------------------
                      '''))
    parser.add_argument('--finalreport', type=str, metavar='<Final_Report.txt>',
                    help=dedent('''\
                        The Final Report (.txt format) need to contain the
                        following columns:
                            - SNP Name
                            - Sample ID
                            - Allele1 - Top
                            - Allele2 - Top
                            - GC Score
                            - Theta
                            - R
                            - GT Score'''))
    parser.add_argument('--dnareport', type=str, metavar='<DNA Report>',
                    help='DNA Report in csv format.')
    parser.add_argument('--pedigreefile', type=str, metavar='<Pedigree>',
                    help=dedent('''\
                        File containing the pedigree of the analyzed population.
                        The file have to respect the following structure:
                            IndidualID    Mother    Father
                        The field separator is the tabulation'''))
    parser.add_argument('--mapfile', type=str, metavar='<Map File>',
                    help=dedent('''\
                        File with the marker position. The file have to
                        contain the following 3 columns: SNP_ID, Chromosome and
                        Position,  separated by a tabulation.'''))
    parser.add_argument('-o', '--out', type=str, metavar='<output file>',
                    default='snp_filter', help='Output file prefix')
    parser.add_argument('--outpath', type=str, metavar='<output path>',
                    default=os.getcwd(), help='Output file prefix')
    parser.add_argument('-t', '--gentrain', default = 0.0, metavar='[0:1]', type=float,
                    help='GenTrain Treshold [default: 0.4]')
    parser.add_argument('-c', '--gencall', default = 0.0, metavar='[0:1]',type=float,
                    help='Threshold GeneCall 10%% (GC10) [Default: 0.2]')
    parser.add_argument('-n', '--missing', default = 0.05, metavar='[0:1]',type=float,
                    help='proportion of allowed missing data. [Default: 0.05 ]')
    parser.add_argument('-d', '--tolerance', default = 0.05, metavar='[0:1]',type=float,
                    #action=MaxMinAction,
                    help=dedent('''\
                        Tolerance for the call rate variation from the
                        average call rate of the considered population.
                        Default: 0.1'''))
    parser.add_argument('--rare', dest='rare',default = 0.05, metavar='[0:1]',type=float,
                    help='Frequency to define an allele as rare')
    parser.add_argument('-e', '--pvalue', default = 0.001, metavar='[0:1]',type=float,
                    help= 'p-value for the Chi-square test [Default = 0.001].')
    parser.add_argument('-u', '--thresh_out', default = 0.001, metavar='[0:1]',type=float,
                    help= 'No. of accepted unexpected genotypes. [Default: 0.001]')
    parser.add_argument('--null_R', dest='null_R', default = 0.3,type=float,
                    metavar='[0:1]',help= dedent('''\
                        Maximum intensity of R to call the genotype as
                        Null-allele. [Default: 0.3]'''))
    parser.add_argument('--pop', choices= ['F1', 'F2', 'Germplasm'],
                    help = dedent('''\
                        Population type. Choices: F1, F2 or Germplasm.
                        For BCx or CP population type use F1.'''))
    parser.add_argument('--p1', type=str, help='Parent 1')
    parser.add_argument('--p2', type=str, help='Parent 2')
    parser.add_argument('--uninformative', action='store_true',
                    help='Group identical marker beforr to create the input file')
    parser.add_argument('--dist_regions', action='store_true',
                    help=dedent('''\
                        Identify regions with more than two consecutive
                        markers distorted'''))
    parser.add_argument('--recover-efxeg', action='store_true',#StoreTrueIf, nargs=0,
                    dest='efxeg', help=dedent('''\
                        Find and score ef x eg markers. Valid only associated
                        with F1/CP.'''))
    parser.add_argument('--recover-null',  action='store_true',#StoreTrueIf,, nargs=0,
                    dest='null',help=dedent('''\
                        Find and score markers with the null-allele. Valid only
                        associated with F1/CP and F2 population.'''))
    parser.add_argument('--outcross', nargs='*', #action=OutcrossAction,
                    default = [], help=dedent('''\
                        List of outcrossing provided as: --outcross ind1 .. indN
                        without space. If not provided the script looks for
                        outcrossing based on No. of accepted unexpected
                        genotypes'''))
    parser.add_argument('--chr', type=int, default=17,
                    help='Number of chromosome of the species. Default: 17')
    parser.add_argument('--logfile', type=str, default = os.path.join(os.getcwd(), 'snp_filter_pipeline.log'))
    return parser


class InputReader(object):#Ui_errorDialog):
    '''
    Class that parse of the input files. '''
    def __init__(self, par, ui):
        '''
        Definition of the initial value of the variable valid for all the def in this class
        '''
#         self.ui = Ui_errorDialog()
#         self.ui.setupUi(self)
            
        self.args = par
        self.snps_number = 0
        self.samples = 0
        self.ind_order = {}
        self.offsprings = {} #possible values: 0 #offspring, 1 #outcross, 2 #other, 3 #failed, 4 #parent
        self.genome_pos = dict()
        self.subpop_s_e = [0,0]
        self.othe_s_e = [0,0]
        #Text console from the Assist Main windows 
        self.ui = ui
        
        #self.consoleTextEdit = consoleTextEdit


    def read_pedigree(self):
        '''
        Read the pedigree file and return a dictionary with the classification of the 
        individuals. The possible class are:
            ::0: offspring, germplasm to analyze
            ::1: outcross
            ::2: other
            ::3: poor quality DNA (defined from "dna_report")
            ::4: parent (defined from arguments)
        If the population type is "Germplasm" all the individuals not in outcrossing are included 
        in the analysis. '''
        try:
            for line in open(self.args.pedigreefile):
    #             line = line.decode('utf-8')
                line = line.rstrip('\r\n').split('\t')
                #skip empty line and header if start with Name
                if len(line) == 0 or line[0] == 'Name':
                    continue
                if self.args.pop in ('F1', 'F2'):
                    if line[1] == self.args.p1 and line[2] == self.args.p2:
                        if line[0] in self.args.outcross:
                            self.offsprings[line[0]]= 1 #outcross
                        else:
                            self.offsprings[line[0]]= 0 #offspring
                    elif line[0] in (self.args.p1, self.args.p2):
                        self.offsprings[line[0]] = 4 #parent
                    else:
                        if line[0] in self.args.outcross:
                            self.offsprings[line[0]]= 1 #outcross
                        else:
                            self.offsprings[line[0]] = 2 #other
                elif self.args.pop == 'Germplasm':
                    if line[0] in self.args.outcross:
                        self.offsprings[line[0]]= 1 #outcross
                    else:
                        self.offsprings[line[0]] = 0 #germplasm to analyze
                        #print line[0], type(line[0]), line[0].encode('utf-8')
        except Exception as inst:
            msg = ('\n--> ERROR:\nError while reading pedigree file (' + inst.message + '). \n'
                   'Check the pedigree file format, re-import the file and rerun.')
            self.ui.consoleTextEdit.setTextColor(QColor("red"))
            self.ui.consoleTextEdit.append(msg)
            self.ui.consoleTextEdit.setTextColor(QColor("black"))
            self.ui.consoleTextEdit.append("")

            self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
            raise
        

    def read_map_file(self, delimiter='\t', comment='#'):
        '''
        Create a dictionary with the positions of all the SNPs present in
        the SNP array.
        The map file have to contain at least three fields "SNPid", "chromosome" and "position". 
        The first line of the file is always skipped. Additionally lines starting with "#" or "\\" 
        are skipped. Field separator: tabulation.
        
        Returns: 
            ::genome_pos: Dictionary where keys are the SNPid and the values are a list with
                          chromosome and position'''
        try:
            for i, line in enumerate(open(self.args.mapfile)):
                line = line.rstrip('\r\n').split(delimiter)
                if 'snpid' in line[0].lower() or i == 0 or comment in line[0]:
                    continue
                self.genome_pos[line[0]] = [float(line[1]), float(line[2])]
            return self.genome_pos
        except Exception as inst:
            msg = ('\n--> ERROR:\nError while reading Map file (' + inst.message + '). \n'
                   'Check the Map file format, re-import the file and rerun.')
            self.ui.consoleTextEdit.setTextColor(QColor("red"))
            self.ui.consoleTextEdit.append(msg)
            self.ui.consoleTextEdit.setTextColor(QColor("black"))
            self.ui.consoleTextEdit.append("")
            self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
            raise


    def dna_report(self):
        '''
        Read the DNA Report from Genome Studio and classifies the DNA that were
        genotyped based on call rate and gencall 10%.
        
        Returns:
            ::offsprings: dictionary where keys are the individuals ID and the values are a list 
                          with individual class (see "read_pedigree"), #No_Calls, #Calls, Call_Rate,
                          10%_GC_Score.'''
        try:
            c_rates = list()
            dnas = dict()
            dnar = csv.reader(open(self.args.dnareport), delimiter=',',
                              quoting=csv.QUOTE_NONE)
            for i, line in enumerate(dnar):
                if i == 2:
                    if line[1] == "DNA_ID" and line[2] == "DNA_Name":
                        delta = 0
                    elif not "DNA_Name" in line:
                        delta = -1
                        
                elif re.compile('^\d').search(line[0]):#line start with a number
                    if (line[-2] in ('NaN', 'nan', 'Non un numero reale')): # Completely failed DNAs
                        self.offsprings[line[1]] = [3, self.snps_number, '0', '0'] #failed
                        c_rates.append(0.)
                    else:
                        try:
                            callrate = float(line[4+delta]) / self.snps_number
                            c_rates.append(callrate)
                            if float(line[-2]) > 1.:
                                gc10 = str(line[-3]) + '.' + str(line[-2])
                            else:
                                gc10 =line[-2]
                            try:
                                tmp_v = self.offsprings[line[1]]
                            except KeyError:
                                tmp_v = 2 #other
                            self.offsprings[line[1]] = [tmp_v,int(line[3+delta]),int(line[4+delta]),
                                                        float(callrate),float(gc10)]
                        except ValueError:
                            self.offsprings[line[1]] = [3, self.snps_number, '0', '0'] #failed
            #print 'prima',self.offsprings
            avgcallrate = sum(c_rates) / len(c_rates)
            del c_rates
            min_c_rate = avgcallrate - self.args.tolerance
            for dna in dnas:
                if (self.offsprings[dna][3] < min_c_rate
                     or self.offsprings[dna][4] < self.args.gencall):
                    self.offsprings[dna][0] = 3 #failed
            for k, v in self.offsprings.items():
                try:
                    iter(v)
                except TypeError:
                    del self.offsprings[k]
            #print 'dopo',self.offsprings
            return self.offsprings
        except:
            msg = '\n--> ERROR:\nCheck the DNA Report Format, re-import the file and rerun.'
            self.ui.consoleTextEdit.setTextColor(QColor("red"))
            self.ui.consoleTextEdit.append(msg)
            self.ui.consoleTextEdit.setTextColor(QColor("black"))
            self.ui.consoleTextEdit.append("")
            self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
            raise

    def final_report_header(self):
        '''
        Read the final report header. 
        Returns:
            ::snps_number: Number of assayed SNPs.'''
        try:
            with open(self.args.finalreport) as fin:
                for line in fin:
                    line = line.rstrip('\r\n')
                    if 'Num SNPs' in line:
                        self.snps_number = int(line.split()[-1])
                    elif 'Num Samples' in line:
                        self.samples = int(line.split()[-1])
                    elif '[Data]' in line:
                        break
            return self.snps_number
        except:
            msg = '\n--> ERROR:\nCheck the Final Report format, re-import the file and rerun.'
            self.ui.consoleTextEdit.setTextColor(QColor("red"))
            self.ui.consoleTextEdit.append(msg)
            self.ui.consoleTextEdit.setTextColor(QColor("black"))
            self.ui.consoleTextEdit.append("")
            self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
            raise
        
        
    def _get_gtypes(self, line, tmp_gtypes,tmp_theta,tmp_r):
        '''
        Retrieve the individuals genotypes from the Final report.
        
        Returns:
            ::tmp_gtypes: SNP genotype
            ::tmp_theta: Theta value
            ::tmp_r: R value '''
        try:
            gt = ''.join(line[2:4])
            #print(self.ind_order,line, gt, tmp_gtypes)
            tmp_gtypes[self.ind_order[line[1]]] = gt
            tmp_theta[self.ind_order[line[1]]] = float(line[6])
            tmp_r[self.ind_order[line[1]]] = float(line[7])
        except KeyError:
            pass
        except:
            msg =  ('\nError while getting genotypes! \nError type: ' + str(sys.exc_info()[0])  )
            self.ui.consoleTextEdit.setTextColor(QColor("red"))
            self.ui.consoleTextEdit.append(msg)
            self.ui.consoleTextEdit.setTextColor(QColor("black"))
            self.ui.consoleTextEdit.append("")
            self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
            raise
        return tmp_gtypes,tmp_theta,tmp_r


    def _get_ind_order(self):
        '''
        Define the individuals order according to the population type.
        
        Returns:
            ::ind_order: dictionary with individual_ID as keys and individual position as value
            ::subpop_s_e: index of first and last offspring/accession to analyze
            ::othe_s_e: index of first and last individuals to exclude'''
        
        if self.args.pop == 'Germplasm':
            ind_order, subpop_s_e, othe_s_e = get_individuals_order('-', '--', self.offsprings)
        else:
            ind_order, subpop_s_e, othe_s_e = get_individuals_order(self.args.p1, self.args.p2, self.offsprings)
        return ind_order, subpop_s_e, othe_s_e


    def final_report_data(self, offsprings):
        '''
        Read the FinalReport and make the classification of the SNP, 
        The data are sorted in an order suitable to make the output files.
        The Final Report have to be grouped by SNP.
        The Final Report (.txt format) have to contain the following columns:
            SNP Name, Sample ID, Allele1 - Top, Allele2 - Top, GC Score, GT Score
             Theta, R, [Chr, Position, Cluster Separation]
             
        At each cycle yield the following values:
            snp_id, gtypes, theta, R, gt and gc '''
        try:
            self.offsprings = offsprings
            self.ind_order, self.subpop_s_e, self.othe_s_e = self._get_ind_order()
            gentrain, gencall = list(), list()
            #self.dna_report()
            data_part = False
            #dictionary with gtypes and coord of posively genotyped SNPs:
            #    gtypes[snp] = [p1,p2,off1,..,offn,oth1,..,othn,out1,..outn]
            tmp_gtypes = ['' for _ in self.ind_order]
            tmp_theta = ['nan' for _ in self.ind_order]
            tmp_r = ['nan' for _ in self.ind_order]
            last_snp = None
            with open(self.args.finalreport) as fin:
                first_line = True
                for line in fin:
                    line = line.rstrip('\r\n')
                    if not data_part:
                        if '[Data]' in line:
                            data_part = True
                        continue
                    #Do this only if after [Data]
                    line = line.split('\t')
                    #Check header order and interrupt the execution if somethings wrong
                    if first_line:#line[0] == 'SNP Name':
                        first_line = False
                        expected = ['SNP Name', 'Sample ID', 'Allele1 - Top',
                                    'Allele2 - Top', 'GC Score', 'GT Score',
                                    'Theta', 'R'] #'Chr', 'Position',' Cluster Sep'
                        for i, l in enumerate(expected):
                            if l != line[i]:
                                msg =  ('Unexpected order in Final Report Columns.\n' +
                                             'The expected columns order is:\n' +
                                             ', '.join(expected) + '\n\n' +
                                             'Please correct your Final Report and rerun.')
                                self.ui.consoleTextEdit.setTextColor(QColor("red"))
                                self.ui.consoleTextEdit.append(msg)
                                self.ui.consoleTextEdit.setTextColor(QColor("black"))
                                self.ui.consoleTextEdit.append("")
                                self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
                            

                        continue

                    # Read and store data
                    if not last_snp:
                        last_snp = line[0]
                        gentrain = [line[4]]
                        gencall = [line[5]]
                    if last_snp == line[0]:
                        tmp_gtypes,tmp_theta,tmp_r = self._get_gtypes(line,tmp_gtypes,tmp_theta,tmp_r)
                        try:
                            if float(line[4]) > 0.:
                                gentrain.append(line[4])
                        except:
                            pass
                        try:
                            if float(line[5]) > 0.:              
                                gencall.append(line[5])
                        except:
                            pass
                    elif last_snp != line[0]:
                        try:
                            gt = np.nanmean(np.array(gentrain).astype(np.float64)),
                            gc = np.percentile(np.array(gencall).astype(np.float64), [10])
                        except AttributeError:
                            gt = np.mean(np.array(gentrain).astype(np.float64)),
                            gc = np.percentile(np.array(gencall).astype(np.float64), [10])
    #                     if gt < 0.4:
    #                         print last_snp, gentrain
    #                     if gc < 0.2:
    #                         print last_snp, gencall
                        yield last_snp, tmp_gtypes, tmp_theta, tmp_r, gt, gc                                                                                                
                        gentrain = [line[4]]
                        gencall = [line[5]]
                        tmp_gtypes = ['' for _ in self.ind_order]
                        tmp_theta = ['nan' for _ in self.ind_order]
                        tmp_r = ['nan' for _ in self.ind_order]
                        tmp_gtypes,tmp_theta,tmp_r = self._get_gtypes(line,tmp_gtypes,tmp_theta,tmp_r)
                        last_snp = line[0]
                try:
                    gt = np.nanmean(np.array(gentrain).astype(np.float64)),
                    gc = np.percentile(np.array(gencall).astype(np.float64), [10])
                except AttributeError:
                    gt = np.mean(np.array(gentrain).astype(np.float64)),
                    gc = np.percentile(np.array(gencall).astype(np.float64), [10])
                    #gt, gc = np.nanmean(np.array(gentrain).astype(np.float64)), np.percentile(np.array(gencall).astype(np.float64), [10])
                yield last_snp, tmp_gtypes, tmp_theta, tmp_r, gt, gc 
        except:
            msg = '\n--> ERROR:\nCheck the Final Report format, re-import the file and rerun. '
            self.ui.consoleTextEdit.setTextColor(QColor("red"))
            self.ui.consoleTextEdit.append(msg)
            self.ui.consoleTextEdit.setTextColor(QColor("black"))
            self.ui.consoleTextEdit.append("")
            self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
            raise
        
        
class SnpFilter(InputReader):
    def __init__(self, par,offsprings, ui):
        super(SnpFilter, self).__init__(par, ui)
        self.offsprings = offsprings
        self.ind_order, self.subpop_s_e, self.othe_s_e = self._get_ind_order()
        self.all_miss_pop = float(self.subpop_s_e[1] - self.subpop_s_e[0]) * self.args.missing
        self.all_miss_ful = float(len(self.offsprings)) * self.args.missing
        self.rare_pop = float(self.subpop_s_e[1] - self.subpop_s_e[0]) * self.args.rare
        self.rare_ful = float(len(self.offsprings)) * self.args.rare
        #self.snps_number = snps_number

    
    def NullAllele_value(self, gtypes_w_p, theta_w_p, r_w_p, snp, a1, a2, max_null_value=0.4):
        '''
        Define the best threshold for the null allele maximum intensity for each snp.
        
        Returns:
            ::NullAllele-Failed_value: Float'''
        
        try:
            _, _, r = gtypes_w_p, theta_w_p, r_w_p
            #t_p1, t_p2 = theta_w_p[0], theta_w_p[1]
            #r_p1, r_p2 = r_w_p[0], r_w_p[1]
            # sort individuals by r
            tmp_r = dict((i,r[i]) for i in range(len(r)))
            sorted_tmp_r = sorted(tmp_r.iteritems(), key=itemgetter(1))
            sort_r =  list(i[1] for i in sorted_tmp_r)
            sort_ind = list(i[0] for i in sorted_tmp_r)
            trim_d = 0 #int(len(sort_r)*0.1)
            trim_up = int(len(sort_r)*0.8)
            trimmed_r = np.array(sort_r[trim_d:trim_up])
            trimmed_rind = np.array(sort_ind[trim_d:trim_up])#[trim:-trim])
            #derivate
            derivate = trimmed_r[1:] - trimmed_r[0:-1]
            ind_pairs = list((trim_d+i, trim_d+i-1) for i in range(1,len(trimmed_rind[:])))
            out_threshold = np.percentile(derivate, 95) * 3. # Gap bigger than 3 times the 95th percentile of the distances
            d_otl, lower_ind,  higher_ind = list(), list(), list()
            start, stop = 0, len(derivate)#* 0.7# * 0.3, len(derivate) * 0.7
            for i, d in enumerate(derivate):
                if start < i < stop and d > out_threshold and d >0.1:
                    d_otl.append(d)
                    lower_ind.append(ind_pairs[i][1])
                    higher_ind.append(ind_pairs[i][0])
            if len(d_otl) >0:
                for u,l in zip(higher_ind, lower_ind):
                    upper_r = trimmed_r[u]
                    lower_r = trimmed_r[l]
                    if lower_r <= max_null_value :#and (lower_r + upper_r)/2 <= max_null_value: 
                        return (lower_r + upper_r)/2
                    else:
                        return np.nanmin(r)-0.2
            elif len(d_otl) == 0:
                return np.nanmin(r)-0.2
        except:
            msg = '\n--> ERROR:\nWhile calculating the null allele threshold'
            self.ui.consoleTextEdit.setTextColor(QColor("red"))
            self.ui.consoleTextEdit.append(msg)
            self.ui.consoleTextEdit.setTextColor(QColor("black"))
            self.ui.consoleTextEdit.append("")
            self.ui.consoleTextEdit.moveCursor(QTextCursor.End)   
            raise
   
    
    def count_genotypes(self, gtypes):
        '''
        count the gtypes for each class in the subpopulation and in the full
        dataset
        
        Returns:
            ::gtypes_count_sub: List with the number of Hom1, Het1, Het2, Hom2, null, miss for the 
                                subpopulation
            ::gtypes_count_ful: List with the number of Hom1, Het1, Het2, Hom2, null, miss for the
                                in all the sample
            ::a1: allele 1
            ::a2: allele 2
            '''
        gtypes_count_sub = [0, 0, 0, 0, 0, 0] #Hom1, Het1, Het2, Hom2, null, miss
        gtypes_count_ful = [0, 0, 0, 0, 0, 0] #Hom1, Het1, Het2, Hom2, null, miss
        gt_class = sorted(set(gtypes))
        #define alleles
        a1, a2 = None, None
        tmp = []
        for c in gt_class:
            try:
                if c[0] != c[1]:
                    a1, a2 = c[0], c[1]
            except IndexError:
                continue
            if not c in ('--', 'OO'):
                tmp.append(c)
        if len(tmp) == 1 and tmp[0][0] == tmp[0][1]:
            a1 = a2 = tmp[0][0]
        elif (len(tmp) == 0) or (a1 and a2):
            pass
        elif len(tmp) == 2 and tmp[0] != tmp[1] and tmp[0][0] == tmp[0][1] and tmp[1][0] == tmp[1][1]:
            a1, a2 = tmp[0][0], tmp[1][0]
        else:
            msg =  ('\nFinal Report Error:\nUnexpected Genotypes in count_genotypes.\n'
                         + str(', '.join(map(str,[a1,a2,set(gtypes), tmp]))))
            self.ui.consoleTextEdit.setTextColor(QColor("red"))
            self.ui.consoleTextEdit.append(msg)
            self.ui.consoleTextEdit.setTextColor(QColor("black"))
            self.ui.consoleTextEdit.append("")
            self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
        
        for i, gt in enumerate(gtypes):
            if gt == '':
                continue
            if gt == '--':
                if self.subpop_s_e[0] < i < self.subpop_s_e[1]:
                    gtypes_count_sub[-1] += 1
                gtypes_count_ful[-1] += 1
            elif gt == 'OO':
                if self.subpop_s_e[0] < i < self.subpop_s_e[1]:
                    gtypes_count_sub[-2] += 1
                gtypes_count_ful[-2] += 1
            elif a1 and a2:
                if gt == a1 + a1:
                    if self.subpop_s_e[0] < i < self.subpop_s_e[1]:
                        gtypes_count_sub[0] += 1
                    gtypes_count_ful[0] += 1
                elif gt == a1 + a2:
                    if self.subpop_s_e[0] < i < self.subpop_s_e[1]:
                        gtypes_count_sub[1] += 1
                    gtypes_count_ful[1] += 1
                elif gt == a2 + a2:
                    if self.subpop_s_e[0] < i < self.subpop_s_e[1]:
                        gtypes_count_sub[3] += 1
                    gtypes_count_ful[3] += 1
                else:
                    msg =  ('Unexpected genotypes at count_genotypes\n' +
                                 'all1 = ' + str(a1) + ', all2 = ' + str(a2) + 
                                 'Current gtypes = ' + gt + 
                                 '\nset(gtypes) = ' + str(set(gtypes))) 
                    self.ui.consoleTextEdit.setTextColor(QColor("red"))
                    self.ui.consoleTextEdit.append(msg)
                    self.ui.consoleTextEdit.setTextColor(QColor("black"))
                    self.ui.consoleTextEdit.append("")
                    self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
                    
        return gtypes_count_sub, gtypes_count_ful, a1, a2


    def mendel_errors(self, gtypes_count_sub, snp, gtypes, off_ind, unexpected_gt):
        '''
        Trio consistency analysis. Individuals passing the accepted error rate
        are classified as outcrossing and excluded from the following analysis.
        Only for F1/BC and F2 population. With the Germplasm the pedigree is not
        taken in account
        
        Returns:
            ::unexpected_gt: dictinary with key the individual id and values a list a tuple with 
                             SNP id and genotype. '''
        try:
            AA = gtypes_count_sub[0]
            AB = gtypes_count_sub[1]
            BB = gtypes_count_sub[3]
            if ((AA > 0 and AA <= self.all_miss_pop)
                  or (AB > 0 and AB <= self.all_miss_pop)
                  or (BB > 0 and BB <= self.all_miss_pop)):
                gtp1 = gtypes[0]
                gtp2 = gtypes[1]
                if not gtp1 in ('--', 'OO') and not gtp2 in ('--', 'OO'):
                    al_par= set([gtp1[0],gtp1[1], gtp2[0], gtp2[1]])
                    for i, gt in enumerate(gtypes[2:self.subpop_s_e[1]+1], 2):
                        if (not gt in ('--', 'OO') and
                            ((self.args.pop == 'F1' and ((not gt[0] in al_par and not gt[1] in al_par)
                              or (gtp1[0] != gtp1[1] and gtp2[0] == gtp2[1] and gt[0] == gt[1] != gtp2[0])
                              or (gtp1[0] == gtp1[1] and gtp2[0] != gtp2[1] and gt[0] == gt[1] != gtp1[0])))
                              or (self.args.pop in ('F1', 'F2') and (gtp1 == gtp2 and gtp1[0] == gtp1[1] and gtp1!=gt)))):
                            #or (self.args.pop == 'F2' and ((AA < self.rare_pop or BB < self.rare_pop) and AB < self.rare_pop))):
                            #print(AA,AB, BB)
                            if off_ind[i] in unexpected_gt:
                                unexpected_gt[off_ind[i]].append((snp,gt))
                            else:
                                unexpected_gt[off_ind[i]] = [(snp,gt)]
            return unexpected_gt
        except Exception as inst:
            msg = ('\n--> ERROR:\nError while checking for Mendel error. \n' +
                   inst.message )
            self.ui.consoleTextEdit.setTextColor(QColor("red"))
            self.ui.consoleTextEdit.append(msg)
            self.ui.consoleTextEdit.setTextColor(QColor("black"))
            self.ui.consoleTextEdit.append("")
            self.ui.consoleTextEdit.moveCursor(QTextCursor.End)


    def exclude_outcrossing(self, unexpected_gt, non_failed):
        '''
        Exclude the outcrossing from the exported individuals. 
        
        Returns:
            ::offsprings: dictionary with the individual classification
            ::ind_order: list with the ordered individuals
            ::subpop_s_e: start and end index of the subpopulation individuals in ind_order
            ::othe_s_e: start and end index of the individuals not in subpopulation in ind_order
            '''
        allowed_out = float(non_failed) * self.args.thresh_out
        #print(allowed_out)
        out_cr = []
        new_order = dict()
        sorted_ind = sorted(self.ind_order.items(), key=lambda x:x[1])
        for ind, _ in sorted_ind:
            if self.args.p1 in unexpected_gt or self.args.p2 in unexpected_gt:
                errorMsg =  ('Unexpected gt fall in parents\n' +
                               ', '.join(map(str,[ind, unexpected_gt[ind], self.args.p1, self.args.p2]))) 
                self.ui.consoleTextEdit.setTextColor(QColor("red"))
                self.ui.consoleTextEdit.append(errorMsg)
                self.ui.consoleTextEdit.setTextColor(QColor("black"))
                self.ui.consoleTextEdit.append("")
                self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
               
            if ind in unexpected_gt and len(unexpected_gt[ind])> allowed_out:
                self.offsprings[ind][0] = 1 #outcross
                out_cr.append(ind)
                del self.ind_order[ind]
        #print(out_cr, unexpected_gt)
        if len(out_cr) == 0:
            return self.offsprings, self.ind_order,self.subpop_s_e, self.othe_s_e
        sorted_ind = sorted(self.ind_order.items(), key=lambda x:x[1])
        for i, (ind, _) in enumerate(sorted_ind):
            new_order[ind] = i
        for i, ind in enumerate(out_cr,1):
            new_order[ind] = max(new_order.values()) + i
        self.subpop_s_e[1] = self.subpop_s_e[1]-len(out_cr)
        self.othe_s_e[0] = self.othe_s_e[0]-len(out_cr)
        self.ind_order = new_order
        return self.offsprings, self.ind_order, self.subpop_s_e,self.othe_s_e


    def classification(self, snp_class, gtypes_count_ful, snp, gtypes, gtrain, gcall):
        '''
        classify the SNPs based on the full germplasm.
        The SNPs are classified in:
            - Failed: if gentrain or gencall to low, missing > allowed,
                      freq het > 0.75
            - Null Allele: OO not rare but gentrain and gencall good
            - False SNP: only 1 genotype
            - ShiftedHomo: only 2 cluster useful in F1, F2
            - multiABCluster: efxeg.
            - Robust: 3 cluster, few missing and null, good genetrain and
                      genecall
        Returns:
            ::snp_class: dictionary with snp as key and [aa, ab, bb, nu, nc ,'CLASS'] as value'''
        #snp_cl_list = snp_classes.snp_classes()

        aa, ab, bb = float(gtypes_count_ful[0]), float(gtypes_count_ful[1]), float(gtypes_count_ful[3])
        nc, nu = float(gtypes_count_ful[-1]), float(gtypes_count_ful[-2])
        if (nc > self.all_miss_ful or nu > float(len(self.offsprings)) * 0.75
                or gcall < self.args.gencall):
            #or gtrain < self.args.gentrain or gcall < self.args.gencall):
            snp_class[snp] = [aa, ab, bb, nu, nc ,'Failed']
        elif (nc <= self.all_miss_ful and nu >= self.all_miss_ful):
            snp_class[snp] = [aa, ab, bb, nu, nc ,'NullAllele-Failed']
        elif (aa < self.rare_ful and bb < self.rare_ful and ab > self.rare_ful):
            snp_class[snp] = [aa, ab, bb, nu, nc ,'Paralogs']
        elif ((ab < self.rare_ful and aa > self.rare_ful and bb < self.rare_ful)
               or (ab < self.rare_ful and aa < self.rare_ful and bb > self.rare_ful)):
            snp_class[snp] = [aa, ab, bb, nu, nc ,'False SNP']
        elif ((ab > self.rare_ful and aa > self.rare_ful and bb < self.rare_ful)
               or (ab > self.rare_ful and aa < self.rare_ful and bb > self.rare_ful)):
            snp_class[snp] = [aa, ab, bb, nu, nc , 'ShiftedHomo']
        elif (ab < self.rare_ful and aa > self.rare_ful and bb > self.rare_ful):
            snp_class[snp] = [aa, ab, bb, nu, nc , 'MissingHetero']
        elif (ab > self.rare_ful and aa > self.rare_ful and bb > self.rare_ful):
            snp_class[snp] = [aa, ab, bb, nu, nc , 'Robust']
        else:
            snp_class[snp] = [aa, ab, bb, nu, nc ,'Failed']
#             errorMsg = ' '.join(map(str, ['You lost the following SNP: ', snp, '\nwhere AA = ',  aa,
#                                           ', AB = ', ab, ', BB = ', bb, ', miss= ', nc , ', null = ', nu,
#                                           '\nNo. of Individuals = ', len(self.offsprings)]))
#             self.ui.consoleTextEdit.setTextColor(QColor("red"))
#             self.ui.consoleTextEdit.append(errorMsg)
#             self.ui.consoleTextEdit.setTextColor(QColor("black"))
#             self.ui.consoleTextEdit.append("")
#             self.ui.consoleTextEdit.moveCursor(QTextCursor.End)                              
        
        return snp_class


    def chi_square2class(self, gtypes_count):
        '''
        Make the chi-square with one degree of freedom for the SNP with 2 genotypic class.
        
        Returns:
            ::(chisq, p): chi-squared test statistic and p-value of the test. 
        '''
        if gtypes_count[0] > self.rare_pop:
            homo = gtypes_count[0]
        elif gtypes_count[3] > self.rare_pop:
            homo = gtypes_count[3]
        f_obs = np.array([homo, gtypes_count[1]])
        return chisq(f_obs)


    def chi_square3class(self, gtypes_count):
        '''
        Make the chi-square with one degree of freedom for the SNP with 3 genotypic class.
        Returns:
            ::(chisq, p): chi-squared test statistic and p-value of the test. '''
        f_obs = np.array(gtypes_count[:2] + [gtypes_count[3]])
        f_exp = np.array([sum(f_obs)/4, sum(f_obs)/2,sum(f_obs)/4])
        return chisq(f_obs, f_exp)


    def try_efeg_recovery(self, gtypes_w_p, theta_w_p, r_w_p, snp, a1, a2, gap_multpl = 3.):
        '''
        Check if discontinuity present in theta AB.
        Discontinuity check using the derivates of the distance between
        contigouos data point. If less than 3 point over 95th percentile and they are consecutive
        the snp is classified as efXeg otherwise the abXab classification is mantained.
        
        Returns:
            ::classification: efXeg or abXab,
            ::recoded_genotypes: list with the genotypes
            ::gt_count: list with the nummber of ind in each genotypic class
            ::chi: chi-squared test statistic
            ::p_val: p-value of the test'''
        
        gtypes, theta, r = gtypes_w_p[2:self.subpop_s_e[1]+1], theta_w_p[2:self.subpop_s_e[1]+1], r_w_p[2:self.subpop_s_e[1]+1]
        gt_p1, gt_p2 = gtypes_w_p[0], gtypes_w_p[1]
        t_p1, t_p2 = theta_w_p[0], theta_w_p[1]
        #r_p1, r_p2 = r_w_p[0], r_w_p[1]

        #get theta of heterozygous ind and sort them by theta
        tmp_ab = dict((i,theta[i]) for i in range(len(theta)) if gtypes[i][0] != gtypes[i][1])
        sorted_tmp_ab = sorted(tmp_ab.iteritems(), key=itemgetter(1))
        t_ab =  list(i[1] for i in sorted_tmp_ab)
        ind_ab = list(i[0] for i in sorted_tmp_ab)
        #test normal distribution witin the ab using the Shapiro–Wilk test
        if shapiro(t_ab)[1] > 0.05:
            return 'abxab', gtypes_w_p, None, None, None
        #trim the 10% of the data at each side of the derivate
        trim = int(len(t_ab)*0.1)
        trimmed_t_ab = np.array(t_ab[trim:-trim])
        trimmed_tind_ab = np.array(ind_ab[trim:-trim])
        #derivate
        derivate = trimmed_t_ab[1:] - trimmed_t_ab[0:-1]
        ind_pairs = list((trim+i, trim+i-1) for i in range(1, len(trimmed_tind_ab[:])))
        out_threshold = np.percentile(derivate, 95) * gap_multpl # Gap bigger than 2 times the 95th percentile of the distances
        dist_p1p2 = abs(t_p1-t_p2)
        d_otl, left_ind,  right_ind = [], [], []
        start, stop = len(derivate) * 0.3, len(derivate) * 0.7
        if dist_p1p2 > out_threshold:
            for i, d in enumerate(derivate):
                if start < i < stop and d > out_threshold:
                    d_otl.append(d)
                    left_ind.append(ind_pairs[i][1])
                    right_ind.append(ind_pairs[i][0])
        # number of allowed point between cluster (equivalent to half ofmissing allowed in population)
        miss = self.args.missing * len(theta) / 2.
        theta_left, theta_right= None, None
        if 0 < len(d_otl) <= miss and 0 < abs(min(right_ind)-max(left_ind)) <= miss: 
            #less than miss betwen first and last outlier
            theta_left = t_ab[max(left_ind)]
            theta_right = t_ab[min(right_ind)]
        if theta_left and theta_right:
            tmp_gt = []
            new_p1, new_p2 = gt_p1, gt_p2
            #mark as "ef" the cluster with mother and as "eg" clust with father
            if theta_left < t_p1 < theta_right or theta_left < t_p2 < theta_right:
                return 'abxab', gtypes_w_p, None, None, None
            elif (t_p1 <= theta_left and t_p2 <= theta_left) or (t_p1 >= theta_right and t_p2 >= theta_right): 
                return 'abxab', gtypes_w_p, None, None, None   
            elif  t_p1 <= theta_left:
                tnew_p1, tnew_p2= gt_p1.lower(),gt_p2.upper()
                new_p1, new_p2 = 'ef', 'eg'
                left, right = 'ef', 'eg'
            elif  t_p1 >= theta_right:
                tnew_p1, tnew_p2= gt_p1.upper(),gt_p2.lower()
                new_p1, new_p2 = 'ef', 'eg'  
                left, right = 'eg', 'ef'
            else:
                return 'abxab', gtypes_w_p, None, None, None
            tmp_r_aa = list(r[i] for i in range(len(r)) if gtypes[i] == a1+a1)
            r_mean_aa = np.mean(tmp_r_aa)
            tmp_r_bb = [r[_] for _ in range(len(r)) if gtypes[_] == a2+a2]
            r_mean_bb = np.mean(tmp_r_bb)
            if r_mean_aa > r_mean_bb:
                ee = a1 + a1
                fg = a2 + a2
            elif r_mean_aa < r_mean_bb:
                ee = a2 + a2
                fg = a1 + a1
            else:
                return 'abxab', gtypes_w_p, None, None, None
            jmgtypes[snp] = [new_p1, new_p2]
            for i, gt in enumerate(gtypes):
                if gt in ('--', 'OO'):
                    jmgtypes[snp].append('--')
                    tmp_gt.append(gt)
                elif gt == ee:
                    jmgtypes[snp].append('ee')
                    tmp_gt.append(ee)
                elif gt == fg:
                    jmgtypes[snp].append('fg')
                    tmp_gt.append(fg)
                elif theta_left < theta[i] < theta_right:
                    tmp_gt.append('--')
                elif theta[i] <= theta_left:
                    jmgtypes[snp].append(left)
                    tmp_gt.append(gt.lower())
                elif theta[i] >= theta_right:
                    jmgtypes[snp].append(right)
                    tmp_gt.append(gt.upper())
                else:
                    del jmgtypes[snp]
                    return 'abxab', gtypes_w_p, None, None, None
            #count gt for each class
            gt_count = {}
            for g in jmgtypes[snp][2:]:#tmp_gt:
                if not g in ('OO', '--'):
                    if g in gt_count:
                        gt_count[g] += 1
                    else:
                        gt_count[g] = 1
            # try to make a chi-square with the 4 class (3 df) to keep only the SNPs that segregate in a proper way
            f_obs = np.array(gt_count.values())
            chi, p_val = chisq(f_obs)
            if p_val < self.args.pvalue:
                del jmgtypes[snp]
                return 'abxab',gtypes_w_p, None, None, None
            else:
                return 'efXeg', [tnew_p1, tnew_p2] + tmp_gt, gt_count, chi, p_val
        if snp in jmgtypes:
            del jmgtypes[snp]
        return 'abxab', gtypes_w_p, None, None, None

    #@print_timing
    

    def recover_null(self, gtypes_count_sub, gtypes_w_p, theta_w_p, r_w_p, snp, a1, a2):
        '''
        If Null Allele is present, this function is trying to classify the SNP as Recovered_Null
        deriving the AA, AB, BB, OO from the A*, OO
        Discontinuity check using the derivates of the distance between
        contigouos data point. If less than 3 point over 95th percentile and they are consecutive
        the snp is classified as Null_4_Cluster (efXeg) or Null_2_Cluster(dominant marker)
        otherwise the NullAllele-Failed or Failed classification is used.
        Returns:
            ::classification: Null_4_Cluster,Null_2_Cluster, Failed or NullAllele-Failed
            ::recoded_genotypes: list with the genotypes
            ::gt_count: list with the nummber of ind in each genotypic class
            ::chi: chi-squared test statistic
            ::p_val: p-value of the test'''
         
        gt_p1, gt_p2, gtypes = gtypes_w_p[0], gtypes_w_p[1], gtypes_w_p[2:self.subpop_s_e[1]+1]
        theta, _ = theta_w_p[2:self.subpop_s_e[1]+1], r_w_p[2:self.subpop_s_e[1]+1]
        #r_p1, r_p2 = r_w_p[0], r_w_p[1]
        aa, bb, _, _= self._get_genomestudio_alleles(gtypes, theta, a1, a2)
        if (gtypes_count_sub[0] >= self.rare_pop
              and gtypes_count_sub[1] >= self.rare_pop
              and gtypes_count_sub[3] >= self.rare_pop
              and gtypes_count_sub[-2] >= self.rare_pop
              and not 'OO' in (gt_p1, gt_p2) and gt_p1 != gt_p2
              and gt_p1[0] == gt_p1[1] and gt_p2[0] == gt_p2[1]) :#4 cluster parents not OO and parents AA and BB
            tmp_gtypes = []
            jmgtypes[snp] = ['ef', 'eg']
            for gt in gtypes:
                if  gt =='OO': #r[i] < self.args.null_R or
                    jmgtypes[snp].append('ee')
                    tmp_gtypes.append(gt)
                elif gt[0] != gt[1]:
                    jmgtypes[snp].append('fg')
                    tmp_gtypes.append(gt)
                elif gt == '--':
                    jmgtypes[snp].append(gt)
                    tmp_gtypes.append(gt)
                elif gt_p1 == aa:
                    tmp_gtypes.append(gt[0] +'O')
                    if gt == aa:
                        jmgtypes[snp].append('ef')
                    elif gt == bb:
                        jmgtypes[snp].append('eg')
                    else:
                        print 'gt mother AA', gt_p1, gt_p2, gt
                        sys.exit()
                elif gt_p1 == bb:
                    tmp_gtypes.append(gt[0] +'O')
                    if gt == aa:
                        jmgtypes[snp].append('eg')
                    elif gt == bb:
                        jmgtypes[snp].append('ef')
                    else:
                        print 'gt mother BB', gt_p1, gt_p2, gt
                        sys.exit()
                else:
                    del jmgtypes[snp]
                    return 'NullAllele-Failed', gtypes_w_p, None, None, None
                    #print 'gt_p1, gt_p2, gt:', gt_p1, gt_p2, gt
            #count gt for each class
            gt_count = {}
            for g in jmgtypes[snp][2:]:#tmp_gtypes:
                if g in ('--'):
                    continue
                if g in gt_count:
                    gt_count[g] += 1
                else:
                    gt_count[g] = 1
            # try to make a chi-square with the 4 class (3 df) to keep only the SNPs that segregate in a proper way
            f_obs = np.array(gt_count.values())
            chi, p_val = chisq(f_obs)
            if p_val < self.args.pvalue:
                del jmgtypes[snp]
                return 'NullAllele-Failed', gtypes_w_p, None, None, None
            else:
                return ('Null_4_Cluster', [gt_p1[0]+'O', gt_p2[0]+'O'] + tmp_gtypes, 
                        gt_count, chi, p_val)
        elif ((gt_p1[0] != gt_p1[1] and gt_p2 in ('--', 'OO')) 
                or (gt_p1 in ('--', 'OO') and gt_p2[0] != gt_p2[1])):
            return 'NullAllele-Failed', gtypes_w_p, None, None, None
        elif (min([gtypes_count_sub[0], gtypes_count_sub[3]]) < self.rare_pop
              and gtypes_count_sub[1]  < self.rare_pop
              and gtypes_count_sub[-2] >= self.rare_pop):
            #2 cluster
            gt_count = {}
            for g in gtypes:
                if g in ('OO', '--'):
                    continue
                if g in gt_count:
                    gt_count[g] += 1
                else:
                    gt_count[g] = 1
            new_p1, new_p2 = gt_p1, gt_p2
            if gt_p1 in ('--', 'OO') and gt_p2 in ('--', 'OO'):
                return 'NullAllele-Failed', gtypes_w_p, None, None, None
            elif gt_p1 == 'OO' and not gt_p2 in ('--', 'OO'):
                #if the mother is lower than 0.35 WRITE 9b (lmxll in JM)
                tnew_p1, tnew_p2 = gt_p1, gt_p2[0]+'O'
                new_p1, new_p2 = 'nn', 'np'
                homo, het = 'nn', 'np'
            elif gt_p2 == 'OO' and not gt_p1 in ('--', 'OO'):
                #if the father is lower than 0.35 WRITE 9a (lmxll in JM)
                tnew_p1, tnew_p2 = gt_p1[0]+'O', gt_p2
                new_p1, new_p2 = 'lm', 'll'
                homo, het = 'll', 'lm'
            elif not gt_p1 in ('--', 'OO') and not gt_p2 in ('--', 'OO'):
                #if both are higher than 0.35 WRITE 9c (h-xkk in JM)
                tnew_p1, tnew_p2 = gt_p1[0]+'O', gt_p2[0]+'O'
                new_p1, new_p2 = 'h-', 'h-'
                homo, het = 'kk', 'h-'
                
            else:
                return 'NullAllele-Failed', gtypes_w_p, None, None, None
            if homo and het:
                tmp_gtypes = []
                jmgtypes[snp] = [new_p1, new_p2]
                for gt in gtypes:
                    if  gt =='OO':#r[i] < self.args.null_R or
                        jmgtypes[snp].append(homo)
                        tmp_gtypes.append(gt)
                    elif gt == '--' or gt_count[gt] <  self.rare_pop:
                        jmgtypes[snp].append('--')
                        tmp_gtypes.append('--')
                    elif gt in (aa, bb):
                        jmgtypes[snp].append(het)
                        if het == 'h-':
                            tmp_gtypes.append(gt[0]+'-')
                        else:
                            tmp_gtypes.append(gt[0]+'O')
                    else:
                        #Not recoverable
                        del jmgtypes[snp]
                        return 'NullAllele-Failed', gtypes_w_p, None, None, None
                #count gt for each class
                gt_count = {}
                for g in jmgtypes[snp][2:]:#tmp_gtypes:
                    if g in ('--'):
                        continue
                    if g in gt_count:
                        gt_count[g] += 1
                    else:
                        gt_count[g] = 1
                # try to make a chi-square to keep only the SNPs that segregate
                # in a proper way
                if len(gt_count)> 2:
                    del jmgtypes[snp]
                    return 'NullAllele-Failed', gtypes_w_p, None, None, None
#                     print '2cluster III:', gt_count
#                     sys.exit()
                if 'kk' in gt_count:
                    try:
                        f_obs = [gt_count['kk'], gt_count['h-']]
                    except:
                        print gt_count
                    chi, p_val = chisq(f_obs, f_exp=np.array([sum(f_obs)/4, (sum(f_obs)/4)*3]))
                else:
                    f_obs = np.array(gt_count.values())
                    chi, p_val = chisq(f_obs)
                if p_val < self.args.pvalue:
                    del jmgtypes[snp]
                    return 'NullAllele-Failed', gtypes_w_p, None, None, None
                else:
                    
                    return 'Null_2_Cluster', [tnew_p1, tnew_p2]+tmp_gtypes, gt_count, chi, p_val
        elif (((gt_p1[0] == gt_p1[1] and gt_p2[0] != gt_p2[1]) or (gt_p1[0] != gt_p1[1] and gt_p2[0] == gt_p2[1]))
              and (gtypes_count_sub[0] > self.rare_pop and gtypes_count_sub[1] > self.rare_pop 
              and gtypes_count_sub[3] > self.rare_pop)):
            return 'Failed', gtypes_w_p, None, None, None
        elif(min([gtypes_count_sub[0], gtypes_count_sub[3]]) < self.rare_pop
              and gtypes_count_sub[1]  > self.rare_pop
              and gtypes_count_sub[-2] >= self.rare_pop):
            #tipo snp SNP_FB_0849456 che segrega come [0, 46, 0, 15, 20, 2]e parentali sono AAxAB
            
            return 'Failed', gtypes_w_p, None, None, None
        elif ((gt_p1 == 'OO' and gt_p2[0] != gt_p2[1]) 
              or (gt_p1[0] != gt_p1[1] and gt_p2 == 'OO')):
            #count gt for each class
            gt_count = {}
            for g in gtypes:
                if g in ('--'):
                    continue
                if g in gt_count:
                    gt_count[g] += 1
                else:
                    gt_count[g] = 1
            f_obs = np.array(gt_count.values())
            if len(f_obs) == 2:
                chi, p_val = chisq(f_obs)
                return 'NullAllele-Failed', [new_p1, new_p2]+gtypes, gt_count, chi, p_val
            else:
                return 'NullAllele-Failed', gtypes_w_p, None, None, None
        else:
            return 'NullAllele-Failed', gtypes_w_p, None, None, None


    def _get_genomestudio_alleles(self, gtypes, theta, a1, a2):
        '''define genomeStudio AA and BB using the mean Theta'''
        tmp_th = {}
        for i, gt in enumerate(gtypes):
            if not gt in ('OO', '--') and gt[0] == gt[1]:
                if gt in tmp_th:
                    tmp_th[gt].append(theta[i])
                else:
                    tmp_th[gt] = [theta[i]]
        aa, bb = a1 + a1, a2 + a2
        try:
            tmean_aa = np.mean(tmp_th[aa])
        except KeyError:
            tmean_aa = 0.5
        try:
            tmean_bb = np.mean(tmp_th[bb])
        except KeyError:
            tmean_bb = 0.5
        if tmean_aa > tmean_bb:
            aa, bb = bb, aa
        del tmp_th
        return aa, bb, tmean_aa, tmean_bb


    def infer_parents(self, gtypes_count_sub, gtypes):
        '''
        Infer parents genotypes based on offspring segregation'''
        gt_p1 = gtypes[0]
        gt_p2 = gtypes[1]
        inferred = False
        AA, AB, BB = gtypes_count_sub[0], gtypes_count_sub[1], gtypes_count_sub[3]
        if (AA > self.rare_pop and AB > self.rare_pop and BB > self.rare_pop 
            and self.chi_square3class(gtypes_count_sub)[1] >= self.args.pvalue): # -- x AB  or AB x AB (AB x AB)
            if gt_p1 in ('--', 'OO') and gt_p2[0] != gt_p2[1]:
                gtypes[0] = gt_p2
                inferred = 'par1'
            elif gt_p1[0] != gt_p1[1] and gt_p2 in ('--', 'OO') :
                gtypes[1] = gt_p1
                inferred = 'par2'
            elif gt_p1 in ('--', 'OO') and gt_p2 in ('--', 'OO'):
                for gt in gtypes[2:self.subpop_s_e[1]]:
                    if gt[0] != gt[1]:
                        gtypes[0] = gt
                        gtypes[1] = gt
                        inferred = 'both'
            else:
                inferred = False
        elif (AA < self.rare_pop or BB < self.rare_pop
              and self.chi_square2class(gtypes_count_sub)[1] >= self.args.pvalue):
            if gt_p1 in ('--', 'OO') and gt_p2 in ('--', 'OO'): # -- X --
                inferred = False
            elif gt_p1 in ('--', 'OO') and gt_p2[0] != gt_p2[1]: # -- x AB
                for gt in gtypes[2:self.subpop_s_e[1]]:
                    if gt[0] == gt[1]:
                        gtypes[0] = gt
                        inferred = 'par1'
            elif gt_p1 in ('--', 'OO') and gt_p2[0] == gt_p2[1]: # -- x AA
                for gt in gtypes[2:self.subpop_s_e[1]]:
                    if gt[0] != gt[1]:
                        gtypes[0] = gt
                        inferred = 'par1'
            elif gt_p2 in ('--', 'OO') and gt_p1[0] != gt_p1[1]: # AB x --
                for gt in gtypes[2:self.subpop_s_e[1]]:
                    if gt[0] == gt[1]:
                        gtypes[1] = gt
                        inferred = 'par2'
            elif gt_p2 in ('--', 'OO') and gt_p1[0] == gt_p1[1]: # AA x --
                for gt in gtypes[2:self.subpop_s_e[1]]:
                    if gt[0] != gt[1]:
                        gtypes[1] = gt
                        inferred = 'par2'
            else:
                inferred = False
        else:
            inferred = False
        return gtypes, inferred


    def define_pop_class_f1(self,pop_class,gtypes_count_sub,snp,gtypes,theta,r,a1,a2):
        '''
        Classify each SNP based on the segregation in the subpopulation.
            - Robust: if F1, F2: 2 or 3 cluster (AAxAB and ABxAB)
                      if Germplasm: 3 cluster
            - AB_2_sub-clusters: efXeg from multiple AB cluster
            - Null_4_Cluster: Null 4 cluster
            - Null_2_Cluster: Null 3 cluster
            - Monomorphics: 1 cluster
            - DistortedAndUnexpSegreg: if F1, F2: 1 class missing or Chi-square > thresh
                         if Germplasm: HWE Chi-square > thresh
            - Failed: all the other SNPs'''
        
        store_mod_geno = False   
        if ((gtypes[0] in ('--', 'OO') or gtypes[1]  in ('--', 'OO'))
             and gtypes_count_sub[-2] < self.rare_pop 
             and sorted([gtypes_count_sub[0], gtypes_count_sub[1], gtypes_count_sub[3]])[1] > self.rare_pop):
            gtypes, inferred = self.infer_parents(gtypes_count_sub, gtypes)
        else:
            inferred = False  
        #print snp
        if (gtypes_count_sub[-1] > self.all_miss_pop
              or sum(gtypes_count_sub[-2:]) > sum(gtypes_count_sub)*0.75):
            #Failed for missing data
            pop_class[snp] = ['Failed'] + gtypes_count_sub + ['--', '--'] + gtypes[:2]
        elif gtypes_count_sub[-2] > self.rare_pop:
            #NullAllele-Failed
            
            if self.args.null:
                
                clas, gtypes, gtypes_count_null, chi_null, p_val_null = self.recover_null(gtypes_count_sub, gtypes, theta, r, snp, a1, a2)
                if clas in ('Null_2_Cluster', 'Null_4_Cluster'):
                    store_mod_geno = True
                if clas == 'Null_2_Cluster':
                    if 'nn' in gtypes_count_null:
                        gtypes_count_sub[0] = gtypes_count_null['nn'] #nn
                        gtypes_count_sub[1] = gtypes_count_null['np'] #np
                    elif 'll' in gtypes_count_null:
                        gtypes_count_sub[0] = gtypes_count_null['ll'] #ll
                        gtypes_count_sub[1] = gtypes_count_null['lm'] #np
                    elif 'kk' in gtypes_count_null:
                        try:
                            gtypes_count_sub[0] = gtypes_count_null['kk'] #ll
                            gtypes_count_sub[1] = gtypes_count_null['h-'] #np
                            gtypes_count_sub[3] = 0
                        except KeyError:
                            gtypes_count_sub[0] = gtypes_count_null['hh']
                            gtypes_count_sub[1] = gtypes_count_null['hk']
                            gtypes_count_sub[3] = gtypes_count_null['kk']
                    gtypes_count_sub[2], gtypes_count_sub[4]= 0, 0
                    pop_class[snp] = ['Null_2_Cluster'] + gtypes_count_sub + [chi_null, p_val_null] + gtypes[:2]
                elif clas == 'Null_4_Cluster':
                    gtypes_count_sub[-2] = gtypes_count_null['ee'] #null-null
                    gtypes_count_sub[0] = gtypes_count_null['ef']
                    gtypes_count_sub[1] = gtypes_count_null['eg']
                    gtypes_count_sub[2] = gtypes_count_null['fg']
                    gtypes_count_sub[3] = 0
                    pop_class[snp] = ['Null_4_Cluster'] + gtypes_count_sub + [chi_null, p_val_null] + gtypes[:2]
                elif clas == 'Failed':
                    pop_class[snp] = ['Failed'] + gtypes_count_sub + ['--', '--'] + gtypes[:2]
                else:
                    pop_class[snp] = ['NullAllele-Failed'] + gtypes_count_sub + ['--', '--'] + gtypes[:2]
            else:
                pop_class[snp] = ['NullAllele-Failed'] + gtypes_count_sub + ['--', '--'] + gtypes[:2]
        elif (a1 == a2
              or (sorted([gtypes_count_sub[0], gtypes_count_sub[1], gtypes_count_sub[3]])[0] < self.rare_pop
                  and sorted([gtypes_count_sub[0], gtypes_count_sub[1], gtypes_count_sub[3]])[1] < self.rare_pop)):
            #Monomorphic
            pop_class[snp] = ['Monomorphic'] + gtypes_count_sub + ['--', '--'] + gtypes[:2]
        elif ((gtypes_count_sub[0] < self.rare_pop
                or gtypes_count_sub[1] < self.rare_pop
                or gtypes_count_sub[3] < self.rare_pop)
              and (gtypes_count_sub[-2] < self.rare_pop) 
               and ((gtypes[0][0] != gtypes[0][1] and gtypes[1][0] == gtypes[1][1])
                    or (gtypes[0][0] == gtypes[0][1] and gtypes[1][0] != gtypes[1][1]))):
            #AA x AB or AB x AA
            (chi, pval) = list(self.chi_square2class(gtypes_count_sub))
            if pval >= self.args.pvalue:
                pop_class[snp] = ['Robust'] + gtypes_count_sub + [chi, pval] + gtypes[:2]
            else:
                pop_class[snp] = ['DistortedAndUnexpSegreg'] + gtypes_count_sub + [chi, pval] + gtypes[:2]
        elif ((gtypes_count_sub[0] > self.rare_pop
                and gtypes_count_sub[1] > self.rare_pop
                and gtypes_count_sub[3] > self.rare_pop)
              and (gtypes_count_sub[-2] < self.rare_pop) 
               and ((gtypes[0][0] != gtypes[0][1] and gtypes[1][0] == gtypes[1][1])
                    or (gtypes[0][0] == gtypes[0][1] and gtypes[1][0] != gtypes[1][1]))):
            #AA x AB or AB x AA with extra class in genotyping
            # convert to lm x ll or nn x np 
            aa_count, ab_count, bb_count = gtypes_count_sub[0],gtypes_count_sub[1],gtypes_count_sub[3]
            if aa_count <  bb_count:
                (chi, pval) = chisq(np.array([aa_count + ab_count, bb_count]))
            else:
                (chi, pval) = chisq(np.array([aa_count, ab_count + bb_count])) 
            pop_class[snp] = ['DistortedAndUnexpSegreg'] + gtypes_count_sub + [chi, pval] + gtypes[:2]
        elif (gtypes_count_sub[0] >= self.rare_pop and gtypes_count_sub[1] >= self.rare_pop
              and gtypes_count_sub[3] >= self.rare_pop and gtypes_count_sub[-2] < self.rare_pop
              and (gtypes[0][0] != gtypes[0][1] and gtypes[1][0] != gtypes[1][1])):
            avg_theta_aa = np.mean([theta[i] for i in range(2,len(theta[2:])) if gtypes[i] == (a1 + a1)])
            avg_theta_bb = np.mean([theta[i] for i in range(2,len(theta[2:])) if gtypes[i] == (a2 + a2)])
            #3 cluster to separate in abXab and efxeg
            if avg_theta_aa > avg_theta_bb:
                avg_theta_bb, avg_theta_aa = avg_theta_aa, avg_theta_bb
                a1, a2 = a2, a1
            chi, pval = self.chi_square3class(gtypes_count_sub)
            if pval >= self.args.pvalue:
                if self.args.efxeg :
                    clas, gtypes, gtypes_count_ef_eg, chi_efeg, p_val_efeg = self.try_efeg_recovery(gtypes, theta, r, snp, a1, a2)
                    if clas == 'efXeg':
                        store_mod_geno = True
                        try:
                            gtypes_count_sub[1] = gtypes_count_ef_eg['ef'] #ef
                            gtypes_count_sub[2] = gtypes_count_ef_eg['eg'] #eg
                            pop_class[snp] = (['AB_2_sub-clusters'] + gtypes_count_sub
                                               + [chi_efeg, p_val_efeg] + gtypes[:2])
                        except KeyError:
                            pop_class[snp] = ['Robust'] + gtypes_count_sub + [chi, pval]  + gtypes[:2]
                    else:
                        pop_class[snp] = ['Robust'] + gtypes_count_sub + [chi, pval]  + gtypes[:2]
                else:
                    pop_class[snp] = ['Robust'] + gtypes_count_sub + [chi, pval] + gtypes[:2]
            else:
                pop_class[snp] = ['DistortedAndUnexpSegreg'] + gtypes_count_sub + [chi, pval] + gtypes[:2]

        else:
            pop_class[snp] = ['Failed'] + gtypes_count_sub + ['--', '--'] + gtypes[:2]

        return pop_class, gtypes, store_mod_geno, inferred


    def define_pop_class_f2(self,pop_class,gtypes_count_sub,snp,gtypes,theta,r,a1,a2):
        '''
        Classify each SNP based on the segregation in the subpopulation.
            - Robust: if F1, F2: 2 or 3 cluster (AAxAB and ABxAB)
                      if Germplasm: 3 cluster
            - AB_2_sub-clusters: efXeg from multiple AB cluster
            - Null_4_Cluster: Null 4 cluster
            - Null_2_Cluster: Null 3 cluster
            - Monomorphics: 1 cluster
            - DistortedAndUnexpSegreg: if F1, F2: 1 class missing or Chi-square > thresh
                         if Germplasm: HWE Chi-square > thresh
            - Failed: all the other SNPs'''
        (aa, ab, _, bb, nu, nc) = gtypes_count_sub
        f_obs = np.array([aa, ab, bb])
        n = sum(f_obs)
        f_exp = np.array([n/4., n/2., n/4.])
        chi, pval = chisq(f_obs, f_exp = f_exp)
        if (nu + nc < self.all_miss_pop 
              and min([aa, ab, bb]) >= self.rare_pop): # Few mmissing and null and no rare alleles
            if pval >= self.args.pvalue:
                pop_class[snp] = ['Robust'] + gtypes_count_sub + [chi, pval, gtypes[0], gtypes[1]]
            else:
                pop_class[snp] = ['DistortedAndUnexpSegreg'] + gtypes_count_sub + [chi, pval, gtypes[0], gtypes[1]]
        elif (nc < self.all_miss_pop and nu > self.rare_pop): 
            #and min([aa, ab, bb]) >= self.rare_pop): #null allele and 3 non rare cluster
            pop_class[snp] = ['NullAllele-Failed'] + gtypes_count_sub + [chi, pval, gtypes[0], gtypes[1]]
        elif nc > self.all_miss_pop:
            pop_class[snp] = ['Failed'] + gtypes_count_sub + ['--','--',gtypes[0], gtypes[1]]
        elif (nu + nc < self.all_miss_pop 
              and ((aa == 0 and ab == 0) or  (bb == 0 and ab == 0))):
            pop_class[snp] = ['NotInformative'] + gtypes_count_sub + [chi, pval, gtypes[0], gtypes[1]]
        elif (nu + nc < self.all_miss_pop 
              and min([aa, ab, bb]) < self.rare_pop): # ShiftedHomo
            pop_class[snp] = ['OneClassMissing'] + gtypes_count_sub + [chi, pval, gtypes[0], gtypes[1]]
        else:
            pop_class[snp] = ['Failed'] + gtypes_count_sub + [chi, pval, gtypes[0], gtypes[1]]
        return pop_class, None,None,None#gtypes_mod, store_mod_geno, inferred
        
        
    def define_pop_class_germplasm(self,pop_class,gtypes_count_sub,snp,gtypes,theta,r,a1,a2):
        '''
        Classify each SNP based on the segregation in the subpopulation.
            - Robust: Germplasm: 3 cluster
            - NotHWE: 3 cluster but snp freq far from HWE.
            - Null: Null + 3 cluster (be careful to heterozygous null)
            - Monomorphics: 1 cluster
            - DistortedAndUnexpSegreg: if F1, F2: 1 class missing or Chi-square > thresh
                         if Germplasm: HWE Chi-square > thresh
            - Failed: all the other SNPs'''
        (aa, ab, _, bb, nu, nc) = gtypes_count_sub
        p = (2. * aa + ab) / (2. * (aa + ab + bb))
        q = 1. - p
        if p > q:
            p,q=q,p
        f_obs = np.array([aa, ab, bb])
        n = sum(f_obs)
        f_exp = np.array([p**2 * n, 2. * p * q *n, q**2 * n])
        maf = p
        chi, pval = chisq(f_obs, f_exp = f_exp, ddof = 1)
        if (nu + nc < self.all_miss_pop 
              and min([aa, ab, bb]) >= self.rare_pop): # Few mmissing and null and no rare alleles
            if pval >= self.args.pvalue:
                pop_class[snp] = ['Robust'] + gtypes_count_sub + [chi, pval, maf, '--']
            else:
                pop_class[snp] = ['DistortedAndUnexpSegreg'] + gtypes_count_sub + [chi, pval, maf, '--']
        elif (nc < self.all_miss_pop and nu > self.rare_pop): 
            #and min([aa, ab, bb]) >= self.rare_pop): #null allele and 3 non rare cluster
            pop_class[snp] = ['NullAllele-Failed'] + gtypes_count_sub + [chi, pval, maf, '--']
        elif nc > self.all_miss_pop:
            pop_class[snp] = ['Failed'] + gtypes_count_sub + ['--','--','--','--']
        elif (nu + nc < self.all_miss_pop 
              and min([aa, ab, bb]) < self.rare_pop): # OneHomoRare
            if min([aa, ab, bb]) == 0:
                pop_class[snp] = ['ShiftedHomo'] + gtypes_count_sub + [chi, pval, maf, '--']
            else:
                if pval >= self.args.pvalue:
                    pop_class[snp] = ['OneHomozygRare_HWE'] + gtypes_count_sub + [chi, pval, maf, '--']
                else:
                    pop_class[snp] = ['OneHomozygRare_NotHWE'] + gtypes_count_sub + [chi, pval, maf, '--']
        else:
            pop_class[snp] = ['Failed'] + gtypes_count_sub + [chi, pval, maf, '--']
        return pop_class, None,None,None#gtypes_mod, store_mod_geno, inferred

def fill_order(some_list, dict2fill):
    start = max(dict2fill.values()) + 1
    for i, el in enumerate(some_list, start):
        dict2fill[el] = i
    return dict2fill, start, max(dict2fill.values())


def get_individuals_order(p1, p2, offsprings):
    '''
    get the order of the individuals in the gtypes and coordinates
    dictionary'''

    order = {p1: 0, p2: 1}
    off, out, oth = [], [], []
    for ind, val in offsprings.items():
        if val[0] == 0: #Offspring
            off.append(ind)
        elif val[0] == 1: #outcross
            out.append(ind)
        elif val[0] == 2 and not ind in (p1,p2): #other
            oth.append(ind)
    order, start, end = fill_order(sorted(off), order)
    subpop_s_e = [start, end]
    order, start, end = fill_order(sorted(out), order)
    othe_s_e = [start]
    order, start, end = fill_order(sorted(oth), order)
    othe_s_e.append(end)
    return order, subpop_s_e, othe_s_e


def avg_r_class(val, gtypes):
    ''' calculate the mean r for AA, AB, BB in gtypes '''
    val_class = {}
    
    for gt, r_va in zip(gtypes, val):
        if gt in ('OO', '--', ''):
            continue
        try:
            val_class[gt].append(r_va)
        except KeyError:
            val_class[gt] = [r_va]
    means = []
    gt = []
    for _, c in val_class.items():
        means.append(np.mean(c))
        gt.append(_)
    return means, gt


class RunAnalysis(InputReader):#, SnpFilter):
    def __init__(self, args, ui):
        if isinstance(args,dict):
            tmpA, tmpb = args, []
            for el, val in tmpA.items():
                if val in ('On', 'Off'):
                    if val == 'On':
                        tmpb.extend(["--" +str(el)])
                elif el == 'outcross':
                    if len(val)>0:
                        tmpb.extend(['--'+str(el)] + val)
                elif el in ('p1', 'p2') and el in ('-', '--'):
                    continue 
                else:
                    tmpb.extend(["--" +str(el),str(val)])
            self.args = get_args().parse_args(tmpb)
        else:
            self.args = args
        
        #Text console from the Assist Main windows 
        self.ui = ui
        
        super(RunAnalysis, self).__init__(self.args, self.ui)
        
        #  Individuals with order and position of first last offspring
        #  and aoutcross individuals
        self.snp_number = self.final_report_header()
        self.read_pedigree()
        self.offsprings = self.dna_report()
        self.ind_order, self.subpop_s_e, self.othe_s_e = self._get_ind_order()
        #Read the map file if present or create a dict with all missing data.
        if self.args.mapfile:
            self.genome_pos = self.read_map_file()
        else:
            self.genome_pos = None
        #  snp classification in full dataset and subpopulation
        self.pop_class = {}
        self.snp_class = {}
        #  Memory mapped array with snp, gt, theta and r
        tmp_rnd_ident = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
        self.fname_a =  os.path.join(mkdtemp(), tmp_rnd_ident + '.dat')
        #  create file with the right shape
        self.staff2plot = np.memmap(self.fname_a,mode='w+', dtype=object, shape=(self.snp_number))
        #  order of snp in file to plot
        self.memmap_file_order = dict()
        #  genotypes of the recovered snps
        self.recovered_gtypes = dict()
        #  inferred parents
        self.inferred_snp = dict()
        # Personalized null allele threshold.
        # null_allele_threshold[snp_id] = rvalue
        self.null_allele_threshold = dict() 
        
        
    def recode_null_gtypes(self, null_R, gtypes, r, snp= None):
        '''
        recode the gtypes if the maximum R value of the null allele 
        is different from 0.3'''
        for i, _ in  enumerate(gtypes):
            if float(r[i]) <= null_R:
                gtypes[i] = 'OO'
        return gtypes
             
    
    def recode_nullasmissing_gtypes(self, gtypes):
        '''
        recode the gtypes classified as null allele (OO) as 
        missing data(--)'''
        for i, gt in  enumerate(gtypes):
            if gt == 'OO':
                gtypes[i] = '--'
        return gtypes
        
        
    def main_process(self, progressBar, min_avg_r_call = 0.3):#0.45):    
        '''
        run the filtering and prepare the data for the export'''
        
        ra = SnpFilter(self.args, self.offsprings, self.ui)
        prog_bar_step = 0
        #Read Final Report a first time to find the Outcrossing and to classify the SNPs
        unexpected_gt = {} # {ind_id: [snps1..snpN], ..}
        off_ind = dict((v,k) for k, v in  ra.ind_order.items())
        for snp_count, (snp, gtypes, theta, r, gtrain, gcall) in enumerate(self.final_report_data(self.offsprings)):
            theta = map(np.float64, theta)
            r = map(np.float64, r)

            #populate progress bar
            if snp_count%(self.snp_number/50) == 0 and progressBar:
                prog_bar_step +=1
                progressBar.setValue(prog_bar_step)
            gtypes_count_sub, gtypes_count_ful, a1, a2 = ra.count_genotypes(gtypes)
            
            # Define a dynamic threshold for the null allele when the r 70th percentile is higher than 
            # the default value of the null allele (0.3)
            if np.percentile(r, [70])[0] > self.args.null_R:
                tmp_null = ra.NullAllele_value(gtypes, theta, r, snp, None, None)
                if tmp_null:
                    gtypes = self.recode_null_gtypes(tmp_null, gtypes, r)
                    gtypes_count_sub, gtypes_count_ful, a1, a2 = ra.count_genotypes(gtypes)
                    if tmp_null < 0:
                        self.null_allele_threshold[snp] = 0.
                    else:
                        self.null_allele_threshold[snp] = tmp_null
                else:
                    self.null_allele_threshold[snp] = self.args.null_R
                if self.args.pop != 'Germplasm':
                    unexpected_gt=ra.mendel_errors(gtypes_count_sub, snp, gtypes, off_ind, unexpected_gt)
                self.snp_class = ra.classification(self.snp_class, gtypes_count_ful, snp, gtypes, gtrain, gcall)
            else:
                self.snp_class[snp] = gtypes_count_ful + ['Failed']
                self.null_allele_threshold[snp] = self.args.null_R
        
        non_failed = 0
        for snp, val in self.snp_class.items():
            if val [-1] != 'Failed':
                non_failed += 1
        self.offsprings, self.ind_order, self.subpop_s_e,self.othe_s_e = ra.exclude_outcrossing(unexpected_gt, non_failed)
        
        # Read Final Report a second time on on the light of the snp classification
        # and after the outcrossing exclusion 
        global jmgtypes
        jmgtypes={}
        for snp_counter, (snp, gtypes, theta, r, gtrain, gcall) in enumerate(self.final_report_data(self.offsprings)):
            try:
                theta = map(np.float64, theta)
                r = map(np.float64, r)
                
                #populate progress bar
                if snp_counter%(self.snp_number/50) == 0 and progressBar:
                    prog_bar_step +=1
                    progressBar.setValue(prog_bar_step)
                
                #Apply when necessary a dynamic threshold for the null allele R
                if self.null_allele_threshold[snp] != self.args.null_R:
                    gtypes = self.recode_null_gtypes(self.null_allele_threshold[snp], gtypes, r, snp)
                gtypes_count_sub, gtypes_count_ful, a1, a2 = ra.count_genotypes(gtypes)
                # change spurious null allele to missing data
                if gtypes_count_ful[-2] < ra.rare_ful and gtypes_count_sub[-2] < ra.rare_pop:
                    gtypes_count_ful[-1] += gtypes_count_ful[-2]
                    gtypes_count_ful[-2] = 0
                    gtypes_count_sub[-1] += gtypes_count_sub[-2]
                    gtypes_count_sub[-2] = 0
                    gtypes = self.recode_nullasmissing_gtypes(gtypes)
                # Set to Failed the SNPs with a low signal intensity
                if (self.snp_class[snp][-1] != 'Failed' and min(avg_r_class(r, gtypes)[0]) < min_avg_r_call 
                      and (gtypes_count_sub[-2] > ra.rare_pop
                           or   gtypes_count_sub[-1] > ra.rare_pop)):
                    self.snp_class[snp][-1] = 'Failed'
                
                if self.snp_class[snp][-1] in ('Robust', 'NullAllele-Failed', 'ShiftedHomo'):
                    if self.args.pop in ('F1'):
                        
                        self.pop_class, gtypes_mod, store_mod_geno, inferred = ra.define_pop_class_f1(
                                                self.pop_class,gtypes_count_sub, snp, gtypes, theta, r, a1, a2)
                        if store_mod_geno:
                            self.recovered_gtypes[snp] = gtypes_mod
                        if inferred:
                            gtypes = gtypes_mod
                            self.inferred_snp[snp] = inferred
                    elif self.args.pop in ('F2'):
                        #FIXME fix F2
                        self.pop_class, gtypes_mod, store_mod_geno, inferred = ra.define_pop_class_f2(
                                                self.pop_class,gtypes_count_sub,snp,gtypes,theta,r,a1,a2)
                    elif self.args.pop == 'Germplasm':
                        self.pop_class, gtypes_mod, store_mod_geno, inferred = ra.define_pop_class_germplasm(
                                                self.pop_class, gtypes_count_sub, snp,gtypes, theta, r, a1, a2)
                else:
                    
                    #print snp, self.snp_class[snp][-1]
                    if self.snp_class[snp][-1] == 'False SNP':
                        self.pop_class[snp] = ['Monomorphic'] + gtypes_count_sub + ['--', '--'] + gtypes[:2]
                    else:
                        self.pop_class[snp] = ['Failed'] + gtypes_count_sub + ['--', '--'] + gtypes[:2]
    
                #write gtypes, r, theta to binary memory mappped file
                if self.pop_class[snp][0] == 'AB_2_sub-clusters':
                    #print len(gtypes_mod), len(gtypes)
                    tmp_arr =  [np.array(gtypes_mod + gtypes[len(gtypes_mod):]).astype(np.str_),
                                np.array(theta).astype(np.float64),
                                np.array(r).astype(np.float64)]
                else:
                    tmp_arr =  [np.array(gtypes).astype(np.str_),
                                np.array(theta).astype(np.float64),
                                np.array(r).astype(np.float64)]
                self.memmap_file_order[snp] = snp_counter
                self.staff2plot[snp_counter] = tmp_arr
            except IndexError:
                continue
        return progressBar, unexpected_gt, self.offsprings,self.ind_order
    
    
    def populate_result_summary(self, inFilesLabel_2):#snp_class, pop_class, genome_pos, offsprings):
        '''
        return the summary of the performed analysis'''

        num_outbred, num_offspr, num_fail = [], 0, []
        for ind, val in self.offsprings.items():
            #possible values: 0 #offspring, 1 #outcross, 2 #other, 3 #failed, 4 #parent
            if val[0] == 0:
                num_offspr += 1
            elif val[0] == 1:
                num_outbred.append(ind)
            elif val[0] == 3:
                num_fail.append(ind)
        tmp = ['Number of SNPs: ' + str(len(self.snp_class)),
               'Number of offsprings (Outcross): ' + str(num_offspr) + ' (' + str(len(num_outbred)) + ')',
               'Outcrossing(s): ' + ', '.join(num_outbred),
               'Total number of individual in FinalReport: ' + str(len(self.offsprings)),
               'Number of excluded samples (Poor quality DNA): ' + str(len(num_fail)),
               'Excluded samples: '  + ', '.join(num_fail)
               ]
        html_text = ('<p><b>Number of SNPs:</b> '  + str(len(self.snp_class)) + '</p>' 
                     '<p><b>Number of offsprings (Outcross):</b> '  + str(num_offspr) 
                            + ' (' + str(len(num_outbred)) + ')' + '<p>'
                     '<p><b>Total number of individual in FinalReport:</b>  ' 
                            + str(len(self.offsprings)) + '</p>'
                     '<p><b>Number of excluded samples (Poor quality DNA): </b>' 
                            + str(len(num_fail)) +'</p>'
                     )
        if self.args.pop == 'Germplasm':
            approved_cls = snp_classes.classes_approved_Germplasm_snps()
            discarded_cls = snp_classes.classes_discarded_Germplasm_snps()
        else:
            approved_cls = snp_classes.classes_approved_F1F2_snps()
            discarded_cls = snp_classes.classes_discarded_F1F2_snps()
        # Summary SNP classification
        title = (#'SNP by class in full dataset (' + str(len(self.offsprings)) + ' individuals):',
                 'SNP Classification (' + str(num_offspr) + ' individuals):',)
        ind = (0,)#-1,0)
        for i, dic in enumerate([self.pop_class]):#self.snp_class, self.pop_class]):
            class_summ = {}
            tmp.append(title[i])
            html_text += '<p><b>' + title[i] + '</b></p><table style="cellpadding:2">'
            for k, v in dic.items():
                try:
                    class_summ[v[ind[i]]] += 1
                except KeyError:
                    class_summ[v[ind[i]]] = 1
            tmp.append('%31s%7s%4s'%(' ','#','%'))
            html_text += ('<tr><td style="padding:0px 15px 0px 0px;"> </td>'
                          '<td style="padding:0px 15px 0px 0px;"><b>#</b></td>'
                          '<td><b>%</b></td></tr>')
            #Approved SNPs:
            appr_snp=sum([v for k, v in class_summ.items() if k in approved_cls])
            appr_snp_p = float(appr_snp) / sum(class_summ.values()) * 100.
            tmp.append('{0:<31s}{1:5d}  {2:5.1f}'.format('Approved', appr_snp, appr_snp_p))
            html_text += ('<tr><td style="padding:0px 15px 0px 0px;"><b>Approved</b></td>'
                          '<td style="padding:0px 15px 0px 0px;"><b>'+str(appr_snp)+'</b></td>'
                          '<td><b>'+'{0:.1f}'.format(appr_snp_p)+'</b></td></tr>')
            for cl in approved_cls:
                if cl in class_summ and class_summ[cl] > 0:
                    appr_snp= class_summ[cl]
                    appr_snp_p = float(appr_snp) / sum(class_summ.values()) * 100.
                    tmp.append('{0:<31s}{1:5d}  {2:5.1f}'.format('  ' + cl, appr_snp, appr_snp_p))
                    html_text += ('<tr><td style="padding:0px 15px 0px 0px;">'+'- '+cl+'</td>'
                                  '<td style="padding:0px 15px 0px 0px;">'+str(appr_snp)+'</td>'
                                  '<td>'+'{0:.1f}'.format(appr_snp_p)+'</td></tr>')
            #Discarded SNPs:
            appr_snp=sum([v for k, v in class_summ.items() if k in discarded_cls])
            appr_snp_p =float(appr_snp) / sum(class_summ.values()) * 100.
            tmp.append('{0:<31s}{1:5d}  {2:5.1f}'.format('Discarded', appr_snp, appr_snp_p))
            html_text += ('<tr><td style="padding:10px 15px 0px 0px;"><b>Discarded</b></td>'
                          '<td style="padding:10px 15px 0px 0px;"><b>'+str(appr_snp)+'</b></td>'
                          '<td style="padding:10px 15px 0px 0px;"><b>'+'{0:.1f}'.format(appr_snp_p)
                          +'<b></td></tr>')
            for cl in discarded_cls:
                if cl in class_summ and class_summ[cl] > 0:
                    appr_snp= class_summ[cl]
                    appr_snp_p = float(appr_snp) / sum(class_summ.values()) * 100.
                    tmp.append('{0:<31s}{1:5d}  {2:5.1f}'.format('  ' + cl, appr_snp, appr_snp_p))
                    html_text += ('<tr><td style="padding:0px 15px 0px 0px;">'+'- '+cl+'</td>'
                                  '<td style="padding:0px 15px 0px 0px;">'+str(appr_snp)+'</td>'
                                  '<td>'+'{0:.1f}'.format(appr_snp_p)+'</td></tr>')
            html_text += '</table>'
        inFilesLabel_2.setText(html_text) 
        return os.linesep.join(tmp)
      
    def populate_snptable(self):#args, snp_class, pop_class, genome_pos):
        ''' 
        Populate the table with the SNPs information.
        The coloumns are: SNP id, Chr, Pos, Classification, Missing, Null, 
            Hom1, Het1, Het2, Hom2, Chi-Square p-val, GT Par1, GT Par2, Maf.
        GT Par1 and GT Par2 anre filled when population type is F1, F2 and Maf when 
        population tye is Germplasm'''
        snptable = list()
        # Sort the SNP by chr and position
        if not self.genome_pos:
            self.genome_pos = dict((snp, ('0', '0')) for snp in self.pop_class)
        snp_sorted_by_chr = sorted(self.genome_pos.items(), key = lambda (k,v): int(v[0]))
        sorted_snp = sorted (snp_sorted_by_chr, key = itemgetter(1))
        row = list()
        #snp_classes
#         classes = ['Robust', 'ShiftedHomo', 'NullAllele-Failed', 'False SNP',
#                    'Paralogs', 'MissingHetero', 'Failed']
        # population Classes
        classes = ['AB_2_sub-clusters','Null_2_Cluster', 
                   'Null_4_Cluster', 'OneClassMissing', 'OneHomozygRare_HWE',
                   'OneHomozygRare_NotHWE','Robust', 'ShiftedHomo',
                   'NullAllele-Failed', 'DistortedAndUnexpSegreg', 'Monomorphic','Failed', 'NotInformative']
        snptable_index = {}
        warn_printed = False
        for c in  classes:
            for tmp in sorted_snp:
                snp, pos = tmp[0], tmp[1]
                try:
                    pop_cl = self.pop_class[snp][0]
                except KeyError:
                    if not warn_printed:
                        msg = ('\nWARNING:\n Map file (' + str(len(sorted_snp)) + ' SNPs) '
                               'longer than FinalReport ('+ str(self.snp_number) + ' SNPs).\n')
                        self.ui.consoleTextEdit.setTextColor(QColor("red"))
                        self.ui.consoleTextEdit.append(msg)
                        self.ui.consoleTextEdit.setTextColor(QColor("black"))
                        self.ui.consoleTextEdit.append("")
                        self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
                        warn_printed = True
                    continue
                    #raise
                if pop_cl != c:
                    continue
                try:
                    miss, null = self.pop_class[snp][6], self.pop_class[snp][5]
                    gt_count = self.pop_class[snp][1:5]
                    try:
                        p_val = "{:5.3f}".format(float(self.pop_class[snp][8]))
                    except ValueError:
                        p_val = self.pop_class[snp][8]
                    if self.args.pop in ('F1', 'F2'):
                        maf = '--'
                        if snp in self.inferred_snp:
                            if self.inferred_snp[snp] == 'par1':
                                gtp = [self.pop_class[snp][9] + '^', self.pop_class[snp][10]]
                            elif self.inferred_snp[snp] == 'par2':
                                gtp = [self.pop_class[snp][9], self.pop_class[snp][10] + '^']
                            elif self.inferred_snp[snp] == 'both':
                                gtp = [self.pop_class[snp][9] + '^', self.pop_class[snp][10] + '^']
                        else:
                            gtp = self.pop_class[snp][9:]
                    else:
                        gtp = ['--', '--']
                        maf = self.pop_class[snp][9]#(2.*gt_count[0] + gt_count[1])/(2.*sum(gt_count))
                        
                    row = [snp, pos[0], pos[1], pop_cl, miss, null] #snp_cl, pop_cl, miss, null]
                    row.extend(gt_count)
                    row.extend([p_val] + gtp + [maf])
                except KeyError:
                    #print snp, snp_class[snp]
                    row = [snp, pos[0], pos[1]]#, snp_cl]
                    for _ in range(11):
                        row.append('--')
                snptable.append(map(str,row))
                snptable_index[snp] = len(snptable) - 1
                #print '\t'.join(map(str,row))
        return snptable, snptable_index
    
    def display_output(self, progressBar, inFilesLabel_2):#, snp_class, pop_class, genome_pos, memmap_file_order, staff2plot, offsprings,ind_order, subpop_s_e, othe_s_e):
        '''
        show the analysis output on the GUI like table and make the plot when asked'''
        #print self.args
        progressBar, unexpected_gt, self.offsprings, self.ind_order = self.main_process(progressBar)
        output_summary = self.populate_result_summary(inFilesLabel_2)#snp_class, pop_class, genome_pos,offsprings)
        
        #populate the snptable
        snptable, snptable_index = self.populate_snptable()#args, snp_class, pop_class, genome_pos)
        return (output_summary, snptable, snptable_index, self.memmap_file_order, self.staff2plot,self.recovered_gtypes, progressBar, 
                self.offsprings, self.ind_order, unexpected_gt)
   
   
class GslikePlot(RunAnalysis):
    '''
    make the plot when a snp is selected in the snptable'''
    def __init__(self, args, memmap_file_order, staff2plot, offsprings,ind_order):
        super(GslikePlot, self).__init__(args, None)
        self.snp_order = memmap_file_order
        self.data2plot = staff2plot
        self.offsprings = offsprings
        self.ind_order = ind_order
        self.sorted_ind = sorted (self.ind_order.items(), key=lambda x:x[1])
 
        
    def plot_setting(self, snp, ax):
        '''
        set the plot feature'''
        max_y = np.nanmax(self.data2plot[self.snp_order[snp]][2]) 
        max_y += 0.2 * max_y #add sam extra space to y to hold legend without covert points
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(0.0, max_y)
        ax.set_title(snp)
        ax.set_xlabel('Theta')
        ax.set_ylabel('R')
        ax.grid(True)
        return ax
        
           
    def make_plot(self, snp, canvas):
        '''
        Create the plot that reproduce the one produced by Genome Studio'''
        data = self.data2plot[self.snp_order[str(snp)]]
        if str(snp) in jmgtypes:
            x, y, gt = data[1],  data[2], np.array(list(jmgtypes[str(snp)])+list(data[0][len(jmgtypes[str(snp)]):]))
            if not 'ee' in gt:
                gt = data[0]
            
        else:
            x, y, gt = data[1],  data[2], data[0]
        try:
            means, gt_l = avg_r_class(x, gt)
            aa, bb = gt_l[means.index(min(means))], gt_l[means.index(max(means))]
            if aa == bb and min(means) > 0.5:
                aa = 'plut'
        except ValueError:
            aa, bb = 'AA', 'BB'
        if gt[0] in ('ef', 'eg') and gt[1] in ('ef', 'eg'):
            aa, bb = 'ee', 'fg'
            if means[gt_l.index('ee')] > means[gt_l.index('fg')]:
                aa, bb = bb, aa
        else:
            if aa[0] != aa[1]:
                aa = 'plut'
            if bb[0] != bb[1]:
                bb = 'plut'
        if self.args.pop in ('F1', 'F2'):
            xy_p1, xy_p2  =  [0., 0.], [0., 0.] 
        else:
            xy_p1, xy_p2  =  None, None
        xy_aa,  xy_ab, xy_ab1, xy_bb, xy_nc, xy_ot = [[], []], [[], []], [[], []], [[], []], [[], []], [[], []]
        for i, g in enumerate(gt):
            try:
                if i == 0 and self.args.pop in ('F1', 'F2'):
                    xy_p1 = [x[0], y[0]]
                elif i == 1 and self.args.pop in ('F1', 'F2'):
                    xy_p2 = [x[1], y[1]]
                elif self.offsprings[self.sorted_ind[i][0]][0] == 0:#self.subpop_s_e[0] < i < self.subpop_s_e[1]:
                    if gt[0] in ('ef', 'eg') and gt[1] in ('ef', 'eg'):
                        if g == aa:
                            xy_aa[0].append(x[i])
                            xy_aa[1].append(y[i])
                        elif g == bb:
                            xy_bb[0].append(x[i])
                            xy_bb[1].append(y[i])
                        elif g == 'ef':
                            xy_ab[0].append(x[i])
                            xy_ab[1].append(y[i])
                        elif g == 'eg':
                            xy_ab1[0].append(x[i])
                            xy_ab1[1].append(y[i])
                        else:
                            xy_nc[0].append(x[i])
                            xy_nc[1].append(y[i])
                    elif 'O' in gt[0] or 'O' in gt[1]:
                        if g in ('OO','--'):
                            xy_nc[0].append(x[i])
                            xy_nc[1].append(y[i])
                        else:#if g[0] != 'O' and g[1] == 'O':
                            if min(means) > 0.5:
                                xy_bb[0].append(x[i])
                                xy_bb[1].append(y[i])
                            else:
                                xy_aa[0].append(x[i])
                                xy_aa[1].append(y[i])
                    elif '-' in gt[0] and gt[0][0] != gt[0][1]:
                        if y[i] < self.args.null_R:
                            xy_nc[0].append(x[i])
                            xy_nc[1].append(y[i])
                        else:
                            if min(means) > 0.5:
                                xy_bb[0].append(x[i])
                                xy_bb[1].append(y[i])
                            else:
                                xy_aa[0].append(y[i])
                                xy_aa[1].append(x[i])      
                            
                            
                    else:
                        if g == aa:
                            xy_aa[0].append(x[i])
                            xy_aa[1].append(y[i])
                        elif g == bb:
                            xy_bb[0].append(x[i])
                            xy_bb[1].append(y[i])
                        elif g[0] != g[1]:
                            xy_ab[0].append(x[i])
                            xy_ab[1].append(y[i])
                        else:
                            xy_nc[0].append(x[i])
                            xy_nc[1].append(y[i])
                else:
                    xy_ot[0].append(x[i])
                    xy_ot[1].append(y[i])
            except KeyError:
                pass
            except IndexError:
                print(g)
            
                
        canvas.ax.clear()
        col = ['k', 'g', 'r', 'm', 'b', 'c', '0.75', '0.35'] #p1, p2, aa, ab, bb, nc, other
        _1, = canvas.ax.plot(np.array(xy_ot[0]), np.array(xy_ot[1]), ms = 2, marker = 'o', 
                             mfc = col[-1], mec = col[-1], ls = '', label='Not subpop')
        _2, = canvas.ax.plot(xy_nc[0], xy_nc[1], ms = 4, marker = 'o', mfc = col[-2], 
                             mec = col[-2], ls = '', label='NC and OO')
        if gt[0] in ('ef', 'eg') and gt[1] in ('ef', 'eg'):
            _3, = canvas.ax.plot(xy_bb[0], xy_bb[1], ms = 4, marker = 'o', mfc = col[4],
                                 mec = col[4], ls = '', label=bb)
            _4, = canvas.ax.plot(xy_ab[0], xy_ab[1], ms = 4, marker = 'o', mfc = col[3], 
                                 mec = col[3], ls = '', label='ef')
            _4a, = canvas.ax.plot(xy_ab1[0], xy_ab1[1], ms = 4, marker = 'o', mfc = col[5], 
                                  mec = col[5], ls = '', label='eg')
            _5, = canvas.ax.plot(xy_aa[0], xy_aa[1], ms = 4, marker = 'o', mfc = col[2], 
                                 mec = col[2], ls = '', label=aa)
        else:
            _3, = canvas.ax.plot(xy_bb[0], xy_bb[1], ms = 4, marker = 'o', mfc = col[4], 
                                 mec = col[4], ls = '', label='BB')
            _4, = canvas.ax.plot(xy_ab[0], xy_ab[1], ms = 4, marker = 'o', mfc = col[3], 
                                 mec = col[3], ls = '', label='AB')
            _5, = canvas.ax.plot(xy_aa[0], xy_aa[1], ms = 4, marker = 'o', mfc = col[2], 
                                 mec = col[2], ls = '', label='AA')
        if self.args.pop != 'Germplasm':
            _6, = canvas.ax.plot(xy_p2[0], xy_p2[1], ms = 8, marker = 'v', mfc = col[1], 
                                 mec = col[1], ls = '', label='Parent 2')
            _7, = canvas.ax.plot(xy_p1[0], xy_p1[1], ms = 8, marker = '^', mfc = col[0], 
                                 mec = col[0], ls = '', label='Parent 1')
        
        canvas.ax = self.plot_setting(str(snp), canvas.ax)
        if self.args.pop != 'Germplasm':
            ncol=4
        else:
            ncol=3
        handles, labels = canvas.ax.get_legend_handles_labels()
        canvas.ax.legend(handles[::-1], labels[::-1], loc = 'upper right',
                         bbox_to_anchor=[0.95, 0.99], ncol=ncol, 
                         borderaxespad=0., prop={'size':9})
        canvas.draw()


        return canvas
                        
        
class FileExporter(RunAnalysis):
    def __init__(self, args, memmap_file_order, memmap_file, recovered_gtypes, unexpected_gt, 
                 offsprings, ind_order, out_path=None, out_prefix=None, ui=None):
        
        super(FileExporter, self).__init__(args, ui)
#         self.args = args
#         print args
        self.offsprings = offsprings
        self.ind_order = ind_order
        if not out_path:
            self.out_path = self.args.outpath
        else:
            self.out_path = str(out_path)
        if not out_prefix:
            self.out_prefix = self.args.out
        else:
            self.out_prefix = str(out_prefix)
        self.unexpected_gt = unexpected_gt
        self.memmap_file_order = memmap_file_order # order of the snps in the mmap file
        self.memmap_file = memmap_file
        self.recovered_gtypes = recovered_gtypes
        self.sorted_ind = sorted (self.ind_order.items(), key=lambda x:x[1])
        # Create optional file object
        self.locfile = None
        self.flexfile = None
        self.pedfile = None
        #  Memory mapped array with ped file 
        tmp_rnd_ident = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
        self.fname =  os.path.join(mkdtemp(), tmp_rnd_ident +'.ped' )
        #  create file with the right shape
        self.ped_memmap_order = dict()
        if self.args.pop == 'Germplasm':
            self.memmap_shape = len([1 for ind in self.offsprings if self.offsprings[ind][0] == 0]) + 2
        else:
            self.memmap_shape = len([1 for ind in self.offsprings if self.offsprings[ind][0] in (0, 4)])
        self.ped_memmap = np.memmap(self.fname,mode='w+', dtype=object, shape=self.memmap_shape)
        self.mapfile = None
        self.structfile = None
        self.hapmapfile = None
     
    def fill_pedmemmapfile(self, gtypes, snp):
        '''
        fill the memmap array to print the plink ped file to avoid putting all ped the file in memory.
        since the informormation have to be transposed'''
        for i, tmp in enumerate(self.sorted_ind):
            ind, pos = tmp[0], tmp[1]
            try:
                if self.offsprings[ind][0] in (0, 4): #offsprings or parents
                    self.ped_memmap_order[ind] = i
                    gtypes[pos] = gtypes[pos].replace('OO', '--')
                    try:
                        self.ped_memmap[i].append(str(gtypes[pos][0]) + ' ' + str(gtypes[pos][1]))
                    except AttributeError:
                        self.ped_memmap[i] = [str(gtypes[pos][0]) + ' ' + str(gtypes[pos][1])]
                    except:
                        errMsg = 'exception raised at fill_pedmemmapfile' + str(sys.exc_info())
                        self.ui.consoleTextEdit.setTextColor(QColor("red"))
                        self.ui.consoleTextEdit.append(errMsg)
                        self.ui.consoleTextEdit.setTextColor(QColor("black"))
                        self.ui.consoleTextEdit.append("")
                        self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
            except KeyError as err:
                if ind in ('-', '--'):
                    pass
                else:
                    self.ui.consoleTextEdit.setTextColor(QColor("red"))
                    self.ui.consoleTextEdit.append(err)
                    self.ui.consoleTextEdit.setTextColor(QColor("black"))
                    self.ui.consoleTextEdit.append("")
                    self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
        
    
    def make_plinkped(self,  fid = '1', p1 = '0', p2 = '0', sex = '0',pheno='-9'):
        for iid, ind_pos in self.ped_memmap_order.items():
            row = [fid, iid, p1, p2, sex, pheno] + self.ped_memmap[ind_pos]+ ['\n']
            self.pedfile.write('\t'.join(row))
        self.pedfile.close()
 
 
    def make_plinkmap(self, snp):
        self.mapfile.write('\t'.join([str(int(float(snp[1]))), snp[0], str(0), str(int(float(snp[2]))), '\n']))
    
    
    def make_hapmap(self,snp, gtypes):
        '''
        create the hapmap file with the genotypes of the successfully genotyped snps
        'rs','alleles', 'chrom', 'pos', 'strand', 'assembly', 
                    'center', 'protLSID', 'assayLSID', 'panel', 'QCcode'''
       
        row = [snp[0], '0', str(int(float(snp[1]))), str(int(float(snp[2])))] + ['0' for _ in range(7)]
#         print 'hapmap, zucc', row, type(row)
#         print self.sorted_ind
        
        for tmp in self.sorted_ind:
            ind, pos = tmp[0], tmp[1]
            if ind in ('-','--'):
                continue
            #print ind, pos, self.offsprings[ind][0], tmp, self.sorted_ind
            if self.offsprings[ind][0] in (0, 4): #offsprings or parents
                try:
                    gtypes[pos] = gtypes[pos].replace('OO', '--')
                    row.append(str(gtypes[pos]))
                except:
                    errMsg = 'exception raised at make_hapmapfile' + str(sys.exc_info())
                    self.ui.consoleTextEdit.setTextColor(QColor("red"))
                    self.ui.consoleTextEdit.append(errMsg)
                    self.ui.consoleTextEdit.setTextColor(QColor("black"))
                    self.ui.consoleTextEdit.append("")
                    self.ui.consoleTextEdit.moveCursor(QTextCursor.End)
        self.hapmapfile.write('\t'.join(row) + '\n')
        
        
    def make_structure_file(self):
        '''
        create the input file for structure for the posively genotyped SNPs'''
        for iid, ind_pos in self.ped_memmap_order.items():
            row = [iid ]
            for gt in  self.ped_memmap[ind_pos]:
                gt = gt.replace('A', '1')
                gt = gt.replace('C', '2')
                gt = gt.replace('G', '3')
                gt = gt.replace('T', '4')
                row.append(gt)
            self.structfile.write('\t'.join(row) + '\n')
        self.structfile.close()

        
    def make_joinmap_f1_file(self, snp_id, gtypes, recov, tmp):
        '''
        create the input file for joinmap with F1/CP coding for the posively genotyped SNPs'''
        row = [snp_id]
        # parents and offspring gtypes
        #self.subpop_s_e[1]]
        if snp_id in jmgtypes:
            gt_p1, gt_p2 =  jmgtypes[snp_id][0], jmgtypes[snp_id][1]
            gt_offspr = jmgtypes[snp_id][2:]#self.memmap_shape]
            #print len(jmgtypes[snp_id]), len(gtypes[:self.memmap_shape])
        else:
            gt_p1, gt_p2 =  gtypes[0], gtypes[1]
            gt_offspr = gtypes[2:self.memmap_shape]
        for i,gt in enumerate(gt_offspr):
            if gt == 'OO':
                gt_offspr[i] = '--'
        # Convert gtypes to joinmap notation and define  Genotype codes 
        if recov:
            row = [snp_id + ']']
            if gt_p1 == 'h-' and gt_p2 == 'h-':
                row.extend(['<hkxhk>', '(h-,kk)'] + gt_offspr)
            else:
                row.extend(['<'+gt_p1+'x'+gt_p2+'>', ' '] + gt_offspr)              
        else:
            if snp_id in self.inferred_snp:
                row = [snp_id + '^']
            else:
                row = [snp_id]
            if '--' in (gt_p1, gt_p2):
                #TODO: Deal with SNP with inferred parents.
                pass
                #print snp_id, row, gtypes,  self.inferred_snp
            elif gt_p1[0] != gt_p1[1] and gt_p2[0] != gt_p2[1]:
                row.extend(['<hkxhk>', ' '])
                (h, k) = gt_p1
                for gt in gt_offspr:
                    if gt == h+h:
                        row.append('hh') 
                    elif gt == h+k:
                        row.append('hk') 
                    elif gt == k+k:
                        row.append('kk')
                    elif gt in ('--', 'OO'):
                        row.append('--')
                    else:
                        print '<hkxhk>',snp_id,gt_p1, gt_p2, gt
            elif gt_p1[0] != gt_p1[1] and gt_p2[0] == gt_p2[1]:
                row.extend(['<lmxll>', ' '])
                (l, m) = gt_p1
                if gt_p2 == m+m:
                    l, m = m, l 
                for gt in gt_offspr:
                    if gt == l+l :
                        row.append('ll') 
                    elif gt in (l+m, m+l):
                        row.append('lm') 
                    elif gt in ('--', 'OO') or not gt in (l+l, l+m, m+l):#account for spurius mendel error putting a missing data
                        row.append('--')
                    else:
                        print '<lmxll>', l, m, snp_id,gt_p1, gt_p2, gt
            elif gt_p1[0] == gt_p1[1] and gt_p2[0] != gt_p2[1]:
                row.extend(['<nnxnp>', ' '])
                (n, p) = gt_p2
                if gt_p1 == p+p:
                    n, p = p, n 
                for gt in gt_offspr:
                    if gt == n+n:
                        row.append('nn') 
                    elif gt in(n+p, p+n):
                        row.append('np') 
                    elif gt in ('--', 'OO') or not gt in (n+n, n+p, p+n):
                        row.append('--')
                    else:
                        print '<nnxnp>',snp_id,gt_p1, gt_p2, gt, n, p
            else:
                print snp_id, 'elsed', gtypes
        if len(row)>1:
            self.locfile.write('\t'.join(row)+ '\n')    
            
    
    def make_flexqtl_file(self, snp, gtypes):
        '''
        create the input file for flexqtl coding for the posively genotyped SNPs'''
        #=======================================================================
        # File structure:
        # SNPid,IndId1,..,IndIdN
        # SNP01,A B,..,A A
        # ..
        # SNPxx,A A,..,B B
        # SNPs to include: all the SNPs with trustable genotyping (classes:
        #      'Null_2_Cluster','Null_4_Cluster','AB_2_sub-clusters','Robust', 'Monomorphic'(SNP class 
        #      qnot population class)                    
        #=======================================================================
        row = [snp[0]]
        for tmp in self.sorted_ind:
            ind, pos = tmp[0], tmp[1]
            if ind in ('-','--'):
                continue
            if self.offsprings[ind][0] in (0, 4): #offsprings or parents
                try:
                    gtypes[pos] = gtypes[pos].replace('OO', '--')
                    row.append(str(gtypes[pos][0]) + ' ' + str(gtypes[pos][1]))
                except:
                    errMsg = 'exception raised at make_flexqtlfile'
                    print errMsg
        self.flexfile.write(','.join(row) + '\n')   
        
        
    def export_output(self, out_type, snptable, snptable_index, export_states, output_summary):
        '''
        when export button activated export the snptable, the genotypes file and
        the additional file selected'''
        # OUT_TYPE KEYS:
        #'joinmap','plink','structure','flexqtl','hapmap','snpinfo','sum','mendel','gtypes'
        #Print Mendel error to file
        if len(self.unexpected_gt) >0 and out_type['mendel']:
            with open(os.path.join(self.out_path, self.out_prefix + '_mendel_error.tsv'), 'wb') as merr_file:
                #Mendel Error Header
                merr_file.write('\t'.join(['Individual_ID', 'SNP_ID', 'Unexpected_gtypes', 
                                           'gtypes_Parent1', 'gtypes_Parent2']) + '\n')
                for ind, merr in self.unexpected_gt.items():
                    for me in merr:
                        snp_id = me[0]
                        if snp_id in self.recovered_gtypes:
                            gtypes_par = self.recovered_gtypes[snp_id][:2]
                        else:
                            gtypes_par = self.memmap_file[self.memmap_file_order[snp_id]][0][:2]
                        
                        merr_file.write('\t'.join([ind] + list(me) + list(gtypes_par)) + '\n')
        #print summary to file
        if out_type['sum']:
            with open(os.path.join(self.out_path, self.out_prefix + '_summary.txt'), 'wb') as fsum:
                #Print used options
                fsum.write('Input Files:\n')
                fsum.write('Final Report File:\t'+str(self.args.finalreport) + '\n')
                fsum.write('DNA Report File:\t'+str(self.args.dnareport) + '\n')
                fsum.write('Pedigree File:\t'+str(self.args.pedigreefile) + '\n')
                if self.args.mapfile:
                    fsum.write('Map File:\t'+str(self.args.mapfile) + '\n\n')
                else: 
                    fsum.write('Map File:\tNot provided\n\n')
                    
                fsum.write('Parameter Set\n')
                fsum.write('Population type:\t'+str(self.args.pop).replace('F1', 'CP (F1) / BC') + '\n')
                fsum.write('Allowed missing data:\t'+str(self.args.missing)+ '\n')
                fsum.write('Call Rate tolerance:\t'+str(self.args.tolerance) + '\n')
                fsum.write('p-Value (Chi-sq) segregation distortion:\t'+str(self.args.pvalue) + '\n')
                fsum.write('Unexpected genotype threshold:\t'+str(self.args.thresh_out) + '\n')
                fsum.write('Frequency rare allele:\t'+str(self.args.rare) + '\n')
                if self.args.pop != 'Germplasm':
                    fsum.write('Parent 1:\t'+str(self.args.p1) + '\n')
                    fsum.write('Parent 2:\t'+str(self.args.p2) + '\n')
                if len(self.args.outcross) > 0: 
                    fsum.write('Individuals to exclude:\t'+str(self.args.outcross) + '\n')
                fsum.write('Number of chromosomes:\t'+str(self.args.chr) + '\n')
                if self.args.efxeg:
                    fsum.write('AB sub-clusters & Null alleles:\tOn\n\n')
                else:
                    fsum.write('AB sub-clusters & Null alleles:\tOff\n\n')
                fsum.write('Outputs:\n')
                fsum.write('Outputs paths:\t'+str(self.out_path) + '\n')
                fsum.write('Outputs prefix:\t'+str(self.out_prefix) + '\n')
                fsum.write('Exported Outputs:\n\t\t')
                if out_type['sum']: fsum.write(str('Summary, '))
                if out_type['gtypes']: fsum.write(str('Custom gtypes, '))
                if out_type['snpinfo']: fsum.write(str('Custom SNP information table, '))
                if out_type['mendel']: fsum.write(str('Custom Mendel error report, '))
                if out_type['joinmap']: fsum.write(str('JoinMap (.loc), '))
                if out_type['hapmap']: fsum.write(str('HapMap, '))
                if out_type['flexqtl']: fsum.write(str('FQ_DataPrepper, '))
                if out_type['plink']: fsum.write(str('PLINK (.ped, .map), '))
                if out_type['structure']: fsum.write(str('STRUCTURE '))  
                fsum.write('\n\n')        
                fsum.write(output_summary)
        if out_type['snpinfo']:
            snptab_outf = open(os.path.join(self.out_path, self.out_prefix + '_snp_info_table.csv'), 'wb')
            
            snptab_outf.write(','.join(["SNP id","Chr", "Pos", "Exported", "Classification",
                                     "Missing", "Null", "Hom1", "Het1", "Het2","Hom2", 
                                     "Chi-Square p-val", "GT Par1", "GT Par2", "MAF"]) + '\n')
        if out_type['gtypes']:
            gtype_outf = open(os.path.join(self.out_path, self.out_prefix + '_gtypes.csv'), 'wb')
            #Print file header
            inds = []
            first_out = True
            first_out_idx = None
            #print('self.sorted_ind',self.sorted_ind)
            for ind in self.sorted_ind: 
                if not ind[0] in ('-', '--'):
                    if self.offsprings[ind[0]][0] == 1 and first_out:
                        #print('caccamo', ind[0],self.offsprings[ind[0]][0])
                        first_out_idx = ind[1]
                        inds.append('|Outcross|')
                        first_out = False
                    if self.offsprings[ind[0]][0] in (0,1,4):
                        inds.append(ind[0]) 
            if not first_out_idx:
                first_out_idx = 0
                for ind in self.sorted_ind: 
                    if not ind[0] in ('-', '--'):
                        if self.offsprings[ind[0]][0] == 0 and ind[1]+1>first_out_idx:
                            first_out_idx = ind[1]+1
                #ind[0] for ind in self.sorted_ind if not ind[0] in ('-', '--')] 
            gtype_outf.write(','.join(["SNP id", "Chr", "Pos", "Classification","Missing", "Null", 
                                       "Hom1", "Het1", "Het2","Hom2", "Chi-Square p-val"] + inds
                                      
                                      + ['\n']))
        #print('inds',inds, first_out_idx)
        # Create optional file object
        if out_type['joinmap']:
            self.locfile = open(os.path.join(self.out_path, self.out_prefix + '_joinmap.loc'), 'wb')
            if self.args.pop == 'F1':
                # Define parent names and number of individuals
                p1, p2 = 'p1', 'p2'
                for (pid, pos) in self.sorted_ind:
                    if int(pos) == 0:
                        p1 = pid
                    elif int(pos) == 1:
                        p2 = pid
                    if  p1 != 'p1' and p2 != 'p2':
                        break   
                # Count number of loci     
                nloc = 0
                for snp_id in export_states:
                    if export_states[snp_id]:
                        nloc += 1
                # Write Joimap CP header 
                self.locfile.write('name = ' + str(str(p1) + 'X' + str(p2))[:20] + '\n')
                self.locfile.write('popt = CP \n')
                self.locfile.write('nloc = ' + str(nloc) + '\n')
                self.locfile.write('nind = ' + str(self.memmap_shape - 2) + '\n')
                # Write individual order as comment line
                self.locfile.write('\t'.join([';SNPid'] 
                                   + [ sid for (sid, pos) in self.sorted_ind 
                                       if int(pos) < self.memmap_shape]) + '\n')
        if out_type['flexqtl']:
            self.flexfile = open(os.path.join(self.out_path, self.out_prefix + '_FQ_DataPrepper.txt'), 'wb')
            # write header
            self.flexfile.write(','.join(['SNPid'] 
                                + [ sid for (sid, pos) in self.sorted_ind 
                                   if int(pos) < self.memmap_shape]) + '\n')
            
            
        if out_type['hapmap']:
            self.hapmapfile= open(os.path.join(self.out_path, self.out_prefix + '_hpm.txt'), 'wb')
            # write header
            head = ['rs','alleles', 'chrom', 'pos', 'strand', 'assembly', 
                    'center', 'protLSID', 'assayLSID', 'panel', 'QCcode']
            for (ind, pos) in self.sorted_ind:
                try:
                    if self.offsprings[ind][0] in (0, 4): #offsprings or parents
                        head.append(str(ind))
                except KeyError:
                    pass
             
            self.hapmapfile.write('\t'.join(head) + '\n')
        if out_type['plink']:
            self.pedfile = open(os.path.join(self.out_path, self.out_prefix + '_plink_in.ped'), 'wb')
            self.mapfile = open(os.path.join(self.out_path, self.out_prefix + '_plink_in.map'), 'wb')
            #Make plink ped memmap file structure
            for i, tmp in enumerate(self.sorted_ind):
                ind, pos = tmp[0], tmp[1]
                try:
                    if self.offsprings[ind][0] in (0, 4): #offsprings or parents
                        self.ped_memmap_order[ind] = i
                        self.ped_memmap[i] = []
                except KeyError:
                    print 'a', i, ind
                    pass

        if out_type['structure']:
            self.structfile = open(os.path.join(self.out_path, self.out_prefix + '_structure.txt'), 'wb')
        # Sort the SNP by chr and position
        if not self.genome_pos:
            self.genome_pos = dict((snp[0], ('0', '0')) for snp in snptable)
        snp_sorted_by_chr = sorted(self.genome_pos.items(), key = lambda (k,v): int(v[0]))
        sorted_snp = sorted (snp_sorted_by_chr, key = itemgetter(1))
        tmp = []
        for tmp_snp in sorted_snp:
            snp_id = tmp_snp[0]
            try:
                snp = snptable[snptable_index[snp_id]]                
                if out_type['flexqtl'] and snp[3] == 'Monomorphic':
                    try:
                        gtypes = self.memmap_file[self.memmap_file_order[snp_id]][0]
                        self.make_flexqtl_file(snp, gtypes)
                    except IndexError:
                        continue
                if export_states[snp[0]]:
                    if snp_id in self.recovered_gtypes:
                        recov = True
                        
                        gtypes = (list(self.recovered_gtypes[snp_id]) + 
                                  list(self.memmap_file[self.memmap_file_order[snp_id]]
                                                        [0][len(self.recovered_gtypes[snp_id]):]))
                    else:
                        recov = False
                        gtypes = self.memmap_file[self.memmap_file_order[snp_id]][0]
                    if out_type['snpinfo']:
                        snptab_outf.write(','.join(snp[:3] + ['True'] + snp[3:]) + '\n')
                    #Export the SNP to genotypes file 
                    if out_type['gtypes']:
                        if isinstance(first_out_idx, int):
                            if first_out_idx == len(inds):
                                genos = list(gtypes[:first_out_idx])
                            else:
                                genos = list(gtypes[:first_out_idx]) + ['|'] + list(gtypes[first_out_idx:len(inds)])
                        else:
                            genos = gtypes[:len(inds)]
                        gtype_outf.write(','.join(snp[:3] + snp[3:11] + list(genos)) + '\n')
                    #Export to optional file format.
                    if not recov:
                        try:
                            if ((out_type['plink'] or out_type['structure']) 
                                    and 0 < int(float(snp[1])) <= self.args.chr): #only snp in valid chromosomes are printed
                                self.fill_pedmemmapfile(gtypes, snp)
                                if out_type['plink']:
                                    self.make_plinkmap(snp)
                            if out_type['hapmap']:
                                self.make_hapmap(snp, gtypes)
                        except KeyError:
                            pass
                    if out_type['joinmap']:
                        if self.args.pop == 'F1':
                            self.make_joinmap_f1_file(snp_id, gtypes,recov, tmp)
                    if out_type['flexqtl']:
                        self.make_flexqtl_file(snp, gtypes)
                else:
                    if out_type['snpinfo']:
                        snptab_outf.write(','.join(snp[:3] + ['False'] + snp[3:] + [os.linesep]))
            except KeyError:
                pass
        if out_type['joinmap']:
            if self.args.pop == 'F1':
                self.locfile.write(os.linesep + 'individual names:' + os.linesep)
                self.locfile.write('\n'.join([ sid for (sid, pos) in self.sorted_ind[2:]
                                       if int(pos) < self.memmap_shape]) + os.linesep)
        if out_type['plink']:
            self.make_plinkped()
        if out_type['structure']:
            self.make_structure_file()
        # TODO: delete temporary files.
        
        #os.remove(self.fname)
    
if __name__ == "__main__":
    if len(sys.argv)==1:
        print (dedent('''\n
                    **WARNING:**
                     The use of this script is deprecated. To run ASSIsT use the
                     provided GUI that can be accessed by running the ASSIsT_1_00.py.\n
                     '''))
        sys.exit(1)
    else:
        args = get_args().parse_args()
    #output_summary, snptable, snptable_index, self.memmap_file_order, self.staff2plot,self.recovered_gtypes, progressBar, self.offsprings, self.ind_order
    output_summary, snptable, snptable_index, memmap_file_order, memmap_file, recovered_gtypes, progressBar, offsprings, ind_order, unexpected_gt = RunAnalysis(args).display_output(None)
    # Cycle over snptable and fill the tableWidget
    exported_class = ['Null_2_Cluster','Null_4_Cluster','AB_2_sub-clusters','Robust']
    non_exported_class = ['DistortedAndUnexpSegreg','Failed','NullAllele-Failed']
    snp_to_export = dict()
    for ln, snp in enumerate(snptable):
        if snp[4] in exported_class:
            snp_to_export[snp[0]] = True
        else:
            snp_to_export[snp[0]] = False
#     main = GslikePlot(args, self.memmap_file_order, self.memmap_file, self.offsprings, self.ind_order)
#     main.make_plot(SNPname, self.ui.plotWidget.canvas)
    
    out_type = {'joinmap': True, 'flexqtl':False,'hapmap':True,
                'plink':True, 'structure':True}
    print 'Select the file to export in addition to snp_table and gtypes. To export the file enter 1 otherwise 0.'
    print 'Default is 1'
    for ot in out_type:
        msg = 'Do you want to export the ' + str(ot) + ' file: (0/[1])'
        s = raw_input(msg)
        if str(s) in ('1', ''):
            out_type[ot] = True
        else:
            out_type[ot] = False             
    fexpor = FileExporter(args, memmap_file_order, memmap_file, recovered_gtypes, unexpected_gt, offsprings, ind_order)
    fexpor.export_output(out_type, snptable, snptable_index, snp_to_export, output_summary)
    #args, snptable, memmap_file_order, staff2plot, offsprings, ind_order, subpop_s_e)
    #fexpor.export_output()

    
    
