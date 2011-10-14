""" Now that you have the Salis program for calculating relative initiation
rates you don't need your old Translation.py any more *sniff sniff*. This
program reads sequence information into a dictionary, obtains the Salis
translation initiation rates, and then correlates this rate with a measure. """

#import sys
#sys.path.append('/home/jorgsk/phdproject/5primeUTR/')
import os
os.chdir('/home/jorgsk/phdproject/5primeUTR/')
import numpy as np
from numpy import array as ar

import Filereader
import Workhouse
import Energycalc
import MyRBS
import cPickle
import matplotlib.pyplot as plt

# Storage location for cpu-heavy dict
storeAdr = '/home/jorgsk/phdproject/5primeUTR/sequence_data/dictstore'

# Rahmi's dataset begins with +1 and has translation intiation at +32. Feeding
# first 80 elements to RBS calculator.
# Veronika's sets begins with -30 and has translation initiation at +31/32.
# Feeding 20:100 to the RBS calculator.

def ReadAndSaveData():
    """ Reads data and adds RSB info through RSBincorp."""
#    seqdict = Filereader.Rahmi104()
    seqdict = Filereader.NikaCombos()
    seqdict = RBSandRNADNAincorp(seqdict)
    cPickle.dump(seqdict, open(storeAdr, 'w'))

def ReadData():
    seqDict = cPickle.load(open(storeAdr, 'r'))
    return seqDict

def RBSandRNADNAincorp(seqdict):
    """ Add RBS and RNA/DNA energy values to the sequence dictionary. """
    rnadna_limit = 18
    for seq in seqdict:
#        sequence = seqdict[seq]['Sequence'][:80] # first 80 nt Rahmi
        sequence = seqdict[seq]['Sequence'][20:100] # 20 to 100 Veronika
        # RBS
        RBS = MyRBS.MyRBS(sequence)
        seqdict[seq]['TLinitiation'] = RBS
        # RNA/DNA
        RNADNA = Energycalc.RNA_DNAenergy(seq[:rnadna_limit])
        seqdict[seq]['TSrnaDNA'] = RNADNA
    return seqdict

def RateExtracter(seqdict):
    """ Extracting the interesting RBS rates (those starting from nt 32 and
    those with max rate) fro the RBS info in MyRBS. """
    for utr in seqdict:
        scores = seqdict[utr]['TLinitiation'][1]
        scores = ar(scores).transpose()
        pos = scores[0].tolist()
        # Getting the indices of the 29, 30, 31, and 32 initiations
        indxs = [pos.index(el) for el in pos if el in [29,30,31,32]]
        # Saving the initiation point and translation rate
        indxscore = [[pos[ind], scores[1][ind]] for ind in indxs]
        indxscore = indxscore[0] # there was only 1 found! NOTE dangerous 
        seqdict[utr]['32starts'] = indxscore
        # Getting the max score as a separate value. If the max score is for a
        # start site higher than N, choose the first less-than-maxstart start site
        seqdict[utr]['maxscore'] = max(scores[1])
        maxstart = 40
        maxind = scores[1].argmax()
        while scores[0][maxind] > maxstart:
            maxind = maxind - 1
        seqdict[utr]['maxscore_adjusted'] = scores[1][maxind]
    return seqdict

def Correlator(seqdict):
    """ Correlate RBS calculated values with PY. """
    x = [] # Induced
    z = [] # Uninduced
    y1 = [] # Predicted rate
    y2 = [] # Predicted maxrate
    for utr, utrDict in seqdict.iteritems():
        y1.append(utrDict['32starts'][1])
        y2.append(utrDict['maxscore_adjusted'])
        x.append(utrDict['Induced'])
        z.append(utrDict['Uninduced'])
    xycorr = Workhouse.Spearman(x,y2)
    zycorr = Workhouse.Spearman(z,y2)
    print "Induced and RBS rate", xycorr
    print "Uninduced and RBS rate", zycorr
    plt.scatter(y1,x)
    plt.xlabel('RBS calculated translation initiation score', size=20)
    plt.ylabel('Induction values of UTR-variant', size=20)
    plt.title('Correlation between RBS calculated scores and the actual'
              ' induction levels to which the sequences could grow')
    plt.show()

def ScoreScaler(seqdict):
    """ Not finished. Maybe will not be finished?"""
    # Rank both. For each transcript, get the distance between the rankings.
    # Multiply distance by a value k (0,1) and roof the value. Propose a new
    # list. Iteratively sort values if they have same rank. Better yet! Don't
    # multiply the rank, but multiply the original Halis value. Then you
    # probably avoid duplicates. NOTE you don't think this will lead to
    # anything, because of the mixed perturbation of transcription and
    # translation.
    rbs = []
    rna = []
    for utr, utrDict in seqdict.iteritems():
        rbs.append(utrDict['32starts'])
        rna.append(utrDict['TSrnaDNA'])
    rbsRank = ar(rbs).argsort()
    rnaRank = ar(rna).argsort()
    pos = 0
    for utr, utrDict in seqdict.iteritems():
        utrDict['32start_rank'] = rbsRank[pos]
        utrDict['ISrnaDNA_rank'] = rnaRank[pos]

ReadAndSaveData()
seqdict = ReadData()
#seqdict.pop('P1-44')
#seqdict = RateExtracter(seqdict)
#seqdict = ScoreScaler(seqdict)
#plt.close()
#Correlator(seqdict)
