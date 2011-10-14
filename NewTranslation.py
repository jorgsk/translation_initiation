""" Now that you have the Salis program for calculating relative initiation
rates you don't need your old Translation.py any more *sniff sniff*. This
program reads sequence information into a dictionary, obtains the Salis
translation initiation rates, and then correlates this rate with a measure. """

from __future__ import division

from IPython.Debugger import Tracer
debug = Tracer()

import sys
# Appending the Hsalis calcuator programs to path.
sys.path.append('/home/jorgsk/phdproject/5primeUTR/Hallis_Calculator')
import numpy as np
from numpy import array as ar

import Filereader
import Workhouse
import Energycalc
import MyRBS
import cPickle
import matplotlib.pyplot as plt

plt.ion()
#plt.ioff()

from scipy import stats

# Storage location for cpu-heavy dict
storeAdr = '/home/jorgsk/phdproject/5primeUTR/sequence_data/RBSdictstore'

# Rahmi's dataset begins with +1 and has translation intiation at +32. Feeding
# first 80 elements to RBS calculator.
# Veronika's sets begins with -30 and has translation initiation at +31/32.
# NOTE now Veronika's sequence begins at +1
# Feeding 20:100 to the RBS calculator.
# NOTE RBS does not take temperature into account. wtf. aha.. the latest turner
# lab measurements can not be extrapolated outside 37 degrees.

def ReadAndSaveData():
    """ Reads data and adds RSB info through RSBincorp."""
    #seqs = Filereader.Rahmi104()
    #seqs = Filereader.NikaCombos()
    seqs = Filereader.Fried()

    seqs = RBSincorp(seqs)
    cPickle.dump(seqs, open(storeAdr, 'w'))

    return seqs

def ReadData():
    seqDict = cPickle.load(open(storeAdr, 'r'))
    return seqDict

def RBSincorp(seqs):
    """ Add RBS scores to the sequence dictionary. """
    for seq in seqs:
        # Thus each proposed translation start site must be put in the 0, 1, or
        # 2 group. For each group the sum of the initiation rates should be
        # given. If the proposed start site is after the given one, this should
        # be noted. The question is how to penalize translation from
        # out of reading frame start sites! One way is to find where the stop
        # codons of the two strands are located!
        #
        # Important: for the synthetic operon
        # it doesn't matter if the translation initiation is further up in the
        # gene, as long as the reading frame is maintained. IDEA perhaps use the
        # hsalis on the whole gene??? Do this if the normal approach bears no
        # results.  +- 60 from translation start. Salis suggests minimum +50.

        rbs_start = max(seq.TNstart-160,0) # 0 is beginning of mRNA
        rbs_stop = seq.TNstart+160 # The T7 UTR is 135 bp long

        # the frame that RBS calculates from
        rbs_frame = Workhouse.ReadingFrame(rbs_start)

        sequence = seq.sequence[rbs_start:rbs_stop]
        # seq.RBS is a list of instances of a result class in MyRBS

        RBS = MyRBS.MyRBS(sequence, rbs_frame, rbs_start)
        seq.addRBS(RBS)

    return seqs

def RatePresenter(seqs):
    """ Extracting rates from all reading frames fro the RBS info in MyRBS. """
    # Synthetic operon: name
    # Reading frame : X
    # Assumed translation start
    # Proposed start sites and rates in reading frame 0
    # Proposed start sites and rates in reading frame 1
    # Proposed start sites and rates in reading frame 2
    # Length of peptides in non-productive reading frames:
    for seq in seqs:
        RBSstarts = seq.RBSstarts
        RBSrates  = seq.RBSrates
        print 'Synthetic operon:\t\t{0}'.format(seq.name)
        print 'Assumed translation start:\t{0}'.format(seq.TNstart)
        print 'Reading frame at nt {0}: \t{1}'.format(seq.TNstart, seq.frame+1)
        for index, starts in enumerate(RBSstarts):
            print 'Sites + rates in frame {0}:'.format(index+1)
            for ind, start in enumerate(starts):
                print '\t\t\t\t{0}\t Rate {1:g}'.format(start, RBSrates[index][ind])
        print ''

def Correlator(seqs):
    """ Correlate RBS calculated values with protein_level.

    Correlate with
    1) Rate from assumed start
    2) Rate from assumed start relative to all other rates
    3) Rate from all in same reading frame
    4) Rate from assumed start relative to rates in other reading frames
    5) Rate from all in same reading frame relative to rates in other reading
    6) Rates in other reading frames
    frames
    """

    # NOTE from Friederike's data;
    # for the original data, initiation from ATG # 2 correlates better than
    # initiation from #3 : I expect their sum as well.

    x = [] # Induced
    ydict = {'Assumed start': [],
             'Assumed start relative to all other': [],
             'Reading frame of assumed start': [],
             'Assumed start relative to other reading frames': [],
             'Reading frame of assumed start relative to other reading frame': [],
             'Other reading frames': []}

    y1, y2, y3, y4, y5, y6 = ([], [] ,[] ,[] ,[], [])

    names = []
    for seq in seqs:

        ## skip the newest sequences from Friederike
        #if not seq.name.endswith('_'):
        #if seq.name.endswith('_'):
            #continue

        if seq.name.startswith('T7'):
            continue

        names.append(seq.name)

        x.append(seq.induced_mean)

        y1.append(seq.TNstartRate)

        y2.append(seq.TNstartRate/(seq.frames_sum - seq.TNstartRate))

        y3.append(seq.TNstartFrameRate)

        y4.append(seq.TNstartRate/(seq.frames_sum - seq.TNstartFrameRate +
                                    seq.TNstartRate))

        try:
            y5.append(seq.TNstartFrameRate/(seq.frames_sum - seq.TNstartFrameRate))
        except ZeroDivisionError:
            y5.append(0)

        y6.append(seq.frames_sum - seq.TNstartFrameRate)

    # add them to dicts
    ydict['Assumed start'] = y1
    ydict['Assumed start relative to all other'] = y2
    ydict['Reading frame of assumed start'] = y3
    ydict['Assumed start relative to other reading frames'] = y4
    ydict['Reading frame of assumed start relative to other reading frame'] = y5
    ydict['Other reading frames'] = y6

    #fig, axes = plt.subplots(2)
    fig, axes = plt.subplots()
    nr = -1
    for pnr, (ptitle, ydata) in enumerate(ydict.items()):
        #if not ptitle in ['Reading frame of assumed start', 'Assumed start']:
        #if not ptitle in ['Assumed start']:
        if not ptitle in ['Reading frame of assumed start']:
            continue
        nr += 1

        #ax = axes[nr]
        ax = axes

        ax.set_xlabel('Induced')
        ax.set_ylabel('Translation score')

        ax.scatter(x, ydata)

        # spearman corr
        spman = '{0:.2f}, {1:.4f}'.format(*stats.spearmanr(x, ydata))
        ax.set_title(ptitle + ' ' + spman)

    plt.show()

def SeqScreener(seqs):
    # Screen out i) the 1000 prot_lev outlier and the signal peptide sequnces
    newseqs = []
    for seq in seqs:
        ## remove the sequences with signal sequence
        #if 'omp' in seq.name:
            #continue
        #if 'pelB' in seq.name:
            #continue
        #remove outlier
        if seq.name == 'gm-omp_LIV-2':
            print 'hei!'
            continue
        newseqs.append(seq)

    return newseqs

def main():
    #ReadAndSaveData() # comment out after first run

    seqs = ReadData()
    RatePresenter(seqs)

    Correlator(seqs)

    return seqs

seqs = main()
