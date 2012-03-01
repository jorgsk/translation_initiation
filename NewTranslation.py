""" Now that you have the Salis program for calculating relative initiation
rates you don't need your old Translation.py any more *sniff sniff*. This
program reads sequence information into a dictionary, obtains the Salis
translation initiation rates, and then correlates this rate with a measure. """

from __future__ import division

from IPython.Debugger import Tracer
debug = Tracer()

import sys
# Appending the Hsalis calcuator programs to path. Horrible hack.
sys.path.append('Hallis_Calculator')
import numpy as np
import scipy.stats

import Filereader
import os
import Workhouse
import MyRBS
import math
import cPickle
import matplotlib.pyplot as plt
from operator import attrgetter
from subprocess import Popen, PIPE

#plt.ion()
plt.ioff()

from scipy import stats

# Storage location for cpu-heavy dict

# Rahmi's dataset begins with +1 and has translation intiation at +32. Feeding
# first 80 elements to RBS calculator.
# Veronika's sets begins with -30 and has translation initiation at +31/32.
# NOTE now Veronika's sequence begins at +1
# Feeding 20:100 to the RBS calculator.
# NOTE RBS does not take temperature into account. wtf. aha.. the latest turner
# lab measurements can not be extrapolated outside 37 degrees.

def ReadAndSaveData(storeAdr, dset):
    """ Reads data and adds RSB info through RSBincorp."""
    #seqs = Filereader.Rahmi104()
    #seqs = Filereader.NikaCombos()
    if dset == 'Fried':
        seqs = Filereader.Fried()

    elif dset == 'Growing':
        seqs = Filereader.Growing()

    elif dset == 'celB':
        seqs = Filereader.celB()

    elif dset == 'designer_rna':
        seqs = Filereader.designer()

    seqs = RBSincorp(seqs)

    cPickle.dump(seqs, open(storeAdr + '_' + dset, 'w'))

    return seqs

def ReadData(storeAdr, dset):
    seqDict = cPickle.load(open(storeAdr + '_' + dset, 'r'))
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
    for seq in sorted(seqs, key = attrgetter('name'), reverse=True):
        RBSstarts = seq.RBSstarts
        RBSrates  = seq.RBSrates
        print 'Name:\t\t\t\t{0}'.format(seq.name)
        print 'Induced:\t\t\t{0}'.format(seq.induced_mean)
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

        #if seq.name.endswith('_'):
            #continue

        #if seq.name.startswith('T7'):
            #continue

        # skip the 'old' seqs
        #if seq.name.startswith('wt') or seq.name.startswith('syn'):
            #continue
        #if seq.name[0] in ['A', 'B']:
            #continue

        # skip the 'new' seqs
        if seq.name.startswith('His') or seq.name.startswith('T7')\
           or seq.name.startswith('H39') or seq.name.startswith('LII'):
            continue

        names.append(seq.name)

        x.append(seq.induced_mean)
        #x.append(math.log(2, seq.induced_mean))

        #y1.append(seq.TNstartRate)
        y1.append(math.log(seq.TNstartRate, 2))

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

    # for getting reading frame data on screen
    ofil = open('better_numbers.csv', 'wb')
    for (n, t, i) in zip(names, y3, x):
        ofil.write('\t'.join([str(n), str(t), str(i)])+'\n')
        print n, t, i
    ofil.close()

    # Write the sequences in fasta-format; from 0 to + 50 translation start
    ofile = open('0_to_plus_50_TN_start.faq', 'wb')
    for seq in seqs:
        ofile.write('>{0}\n'.format(seq.name))
        ofile.write('{0}\n'.format(seq.sequence[:seq.TNstart+50]))
    ofile.close()

    fig, axes = plt.subplots()
    nr = -1

    for pnr, (ptitle, ydata) in enumerate(ydict.items()):
        #if not ptitle in ['Reading frame of assumed start', 'Assumed start']:
        if not ptitle in ['Assumed start']:
            continue
        #if not ptitle in ['Reading frame of assumed start']:
            #continue
        nr += 1

        #ax = axes[nr]
        ax = axes

        ax.set_xlabel('Induced')
        ax.set_ylabel('RBS score')

        ax.scatter(x, ydata)

        # spearman corr
        spman = '{0:.2f}, {1:.4f}'.format(*stats.spearmanr(x, ydata))
        ax.set_title(ptitle + ' ' + spman)

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

def probability_reader(seqs, plot, downstream_seqlen=100):
    """
    Make folding-probability plots in the order of induction. 5' hairpin, SD and
    ATG are all highlighted
    """
    here = os.path.dirname(os.path.realpath(__file__))
    fold_dir = os.path.join(here, 'folding_profiles')

    #fig_dir_new = os.path.join(here, 'folding_figures', 'new')
    #fig_dir_old = os.path.join(here, 'folding_figures', 'old')

    fig_dir = os.path.join(here, plot,'dstream_{0}'.format(downstream_seqlen))
    if not os.path.isdir(fig_dir):
        os.makedirs(fig_dir)

    for seq in seqs:

        # make special case for /
        if '/' in seq.name:
            seq.name = seq.name.replace('/', '~')

        if fold_dir:
            # 1 dir for each fold
            seqDir = os.path.join(fold_dir, seq.name)
            if not os.path.isdir(seqDir):
                os.mkdir(seqDir)

        outfilePath = os.path.join(seqDir, seq.name)
        outfileHandle = open(outfilePath, 'wb')

        # Save a fasta file with the sequence from +1 to TN_start + 50
        outfileHandle.write('>{0}\n'.format(seq.name))
        outfileHandle.write(seq.sequence[:seq.TNstart+downstream_seqlen])
        outfileHandle.close()

        # run hybrid-ss with X samples
        cmd = ['hybrid-ss', '--suffix', 'DAT', '--tracebacks', '100',
               outfilePath]

        # change directory to the sec dir
        os.chdir(seqDir)
        p = Popen(cmd, stderr=PIPE)
        p.wait()

        # get the single-stranded probabilities
        ss_array = []
        for line in open('.'.join([seq.name, '37', 'ext']), 'rb'):
            try:
                (nr, pSS, pSS_and_next_SS) = line.split('\t')
            # last line has only 2 values 
            except ValueError:
                (nr, pSS) = line.split('\t')
                pSS.rstrip('\n')

            # skip the first line
            if nr == 'i':
                continue

            # normalize a bit the large and small values
            #if float(pSS) < 0.05:
                #pSS = 0
            #if float(pSS) > 0.95:
                #pSS = 1

            ss_array.append(float(pSS))

        # add the binding probabilities
        seq.binding_probabilities = ss_array
        print sum(ss_array[seq.TNstart-15:seq.TNstart+15]), seq.induced_mean

        # add the total binding energy
        seq.deltaG = open(seq.name+'.37.plot', 'rb').next().split()[-1]

        # don't plot unneceeasartyly
        if not plot:
            continue

        xvals = range(1, len(ss_array)+1)

        fig, ax = plt.subplots()

        # plot the whole strtch
        ax.plot(xvals, ss_array, linewidth=2)

        atg = seq.TNstart
        # plot the first 5 nucleotides, ATG-10 and ATG in green, red, and yellow
        ax.plot(range(1,6), ss_array[:5], linewidth=2, c='g')
        ax.plot(range(atg-15+1, atg+1), ss_array[atg-15:atg], linewidth=2, c='r')
        ax.plot(range(atg+1, atg+3+1), ss_array[atg:atg+3], linewidth=2, c='y')

        ax.set_ylim(-0.2, 1.2)
        ax.set_xlim(-1, len(ss_array)+1)

        # set title
        ax.set_title('{0} {1} {2}'.format(seq.name, seq.TNstartRate,
                                          seq.induced_mean))
        image = os.path.join(seqDir, '{0}.png'.format(seq.name))
        fig.savefig(image, format='png')

        # FOR CELB
        # also save in a common dir

        # FOR FRIED
        #True also save in common folder
        #if seq.name.startswith('His') or seq.name.startswith('T7')\
           #or seq.name.startswith('H39') or seq.name.startswith('LII'):
            #fig_dir = fig_dir_new
        #else:
            #fig_dir = fig_dir_old

        image2 = os.path.join(fig_dir, '{0}.png'.format(seq.name))
        fig.savefig(image2, format='png')


    # Change back to working directory
    os.chdir(here)

    return seqs

def fold_similarity(seqs):
    """
    The H39, LII and T7 without His tag have different induced but same in RBS
    score: demonstrate that the overall change in 5'UTR fold changes more for
    these than for T7 with His. These are only 4 examples.

    Then do the same for the 'old' variants; correlate the similarity of total
    fold change of all variants with the change in induced and RBS.

    A possible conclusion is that 5'UTR fold changes not captured by the RBS are
    important for protein expression.

    To get better differences for fold, use a normalized to random sequences
    parameter. Now
    """

    olds = []
    news = []

    # separate seqs in new and old
    for seq in seqs:
        if seq.name.startswith('His') or seq.name.startswith('T7')\
           or seq.name.startswith('H39') or seq.name.startswith('LII'):
            news.append(seq)
        else:
            olds.append(seq)

    # for the 'new' ones, group in 'wt/syn' pairs and compare their 
    groupings = ['His-UTR', 'H39-UTR', 'LII-10-UTR', 'T7-UTR_without_His', 'T7-UTR']

    foldDist = []
    indDist = []
    rbsDist = []

    for name in groupings:
        wt = [s for s in news if s.name == name+'_wt'][0]
        syn = [s for s in news if s.name == name+'_syn'][0]

        # Get the normalized distances between fold, induced, and RBS score at
        # TN start
        x = wt.binding_probabilities
        y = syn.binding_probabilities

        if len(x) != len(y):
            print('LENGTH ERROR')

        # get the 1-norm fold difference
        dFold = sum(abs(a-b) for (a,b) in zip(x, y))/len(x)
        foldDist.append(dFold)

        # get the difference in induced
        dInd = dist1(wt.induced_mean, syn.induced_mean)
        indDist.append(dInd)

        # get the RBS difference
        dRBS = dist1(wt.TNstartRate, syn.TNstartRate)
        rbsDist.append(dRBS)

        plot_diff(wt, syn, name)

    # print a table
    print('')
    print('Name\nInduced-difference\tFold-difference\tRBS-difference')
    print('')
    for idiff, fdiff, rdiff, nme in zip(indDist, foldDist, rbsDist, groupings):
        print nme
        I = format(idiff, '.4f')
        F = format(fdiff, '.4f')
        R = format(rdiff, '.4f')
        print('\t'.join([I, F, R]))
        print('')
    print('')

    fig, axes = plt.subplots(2, sharex=True, sharey=True)
    axes[0].scatter(foldDist, indDist)
    #axes[0].set_title('Fold vs induced distance')
    axes[0].set_ylabel('Difference in protein expression', size=10)
    axes[0].set_xlabel('Difference in over-all secondary structure', size=10)

    axes[1].scatter(rbsDist, indDist)
    #axes[1].set_title('RBS vs induced distance')
    axes[1].set_ylabel('Difference in protein expression', size=10)
    axes[1].set_xlabel('Difference in RBS-rates', size=10)
    fig.subplots_adjust(wspace=2)

    fig.suptitle('Fold of whole 5UTR explains difference in protein expression '
                 'better than RBS calculator', size=13)

    spea = scipy.stats.spearmanr(foldDist, indDist)
    pear = scipy.stats.pearsonr(foldDist, indDist)
    both = spea + pear

    headr = 'Spearman: ({0:.2f}, {1:.4f})  Pearson: ({2:.2f}, {3:.4f})'.format(*both)
    axes[0].set_title(headr)
    print headr

    spea = scipy.stats.pearsonr(rbsDist, indDist)
    pear = scipy.stats.spearmanr(rbsDist, indDist)
    both = spea + pear

    headr ='Spearman: ({0:.2f}, {1:.4f})  Pearson: ({2:.2f}, {3:.4f})'.format(*both)
    axes[1].set_title(headr)
    print headr

    savedir = '/home/jorgsk/phdproject/5UTR/FriedrikeFolds/report/images'
    image = os.path.join(savedir, 'difference_correlation.png'.format(seq.name))
    fig.savefig(image, format='png', dpi=200)


def plot_diff(wt, syn, name):

    savedir = '/home/jorgsk/phdproject/5UTR/FriedrikeFolds/report/images'

    order = [('wt', wt), ('syn', syn)]

    fig, axes = plt.subplots(2, sharex=True, sharey=True)

    for ind, ax in enumerate(axes):

        # first wt, then syn
        wtORsyn, seq = order[ind]

        x = seq.binding_probabilities

        xvals = range(1, len(x)+1)

        # plot the whole strtch
        ax.plot(xvals, x, linewidth=2)

        atg = seq.TNstart
        # plot the first 5 nucleotides, ATG-10 and ATG in green, red, and yellow
        ax.plot(range(1,6), x[:5], linewidth=2, c='g')
        ax.plot(range(atg-10+1, atg+1), x[atg-10:atg], linewidth=2, c='r')
        ax.plot(range(atg+1, atg+3+1), x[atg:atg+3], linewidth=2, c='y')

        ax.set_ylim(-0.2, 1.2)
        ax.set_xlim(-1, len(x)+1)

        # set title
        ax.set_title(wtORsyn)

        if ind == 1:
            ax.set_xlabel('nt position from transcription +1')
        ax.set_ylabel('Probability of nt single stranded', size=10)

    if name == 'T7-UTR':
        name = name + ' (with HIS)'
    fig.suptitle(name, size=20)
    image = os.path.join(savedir, '{0}.png'.format(seq.name))
    fig.savefig(image, format='png', dpi=200)

def dist1(x,y):
    return np.abs(x-y)/norm1(x,y)

def dist2(x,y):
    return np.abs(x-y)/norm2(x,y)

def norm1(x,y):
    """
    Assumes positive x and y
    """
    return (x+y)

def norm2(x,y):
    """
    Euclidian norm
    """
    return np.sqrt(x*x + y*y)

def rbs_correlation(seqs):
    """
    Do the RBS calculator score correlate with induction values?
    """

    savedir = '/home/jorgsk/phdproject/5UTR/FriedrikeFolds/report/images'

    olds = []
    news = []
    new_no_T7 = []

    # separate seqs in new and old
    for seq in seqs:
        if seq.name.startswith('His') or seq.name.startswith('T7')\
           or seq.name.startswith('H39') or seq.name.startswith('LII'):
            news.append(seq)

            if seq.name.startswith('T7') and 'without' not in seq.name:
                pass
            else:
                new_no_T7.append(seq)

        else:
            olds.append(seq)

    for oldnew, seqs in [('old', olds), ('new', news), ('not7', new_no_T7)]:

        x = [s.induced_mean for s in seqs]
        y = [s.TNstartRate for s in seqs]

        fig, ax = plt.subplots()

        ax.set_xlabel('Induced')
        ax.set_ylabel('RBS score')

        ax.scatter(x, y)

        # spearman corr
        spman = 'Spearman: {0:.2f}, {1:.4f}'.format(*stats.spearmanr(x, y))
        print oldnew, spman

        ax.set_title(spman, size=20)

        if oldnew == 'old':
            fig.suptitle('CDS variants')
            imagePath = os.path.join(savedir, 'CDS_vs_RBS.png'.format(seq.name))
            fig.savefig(imagePath, format='png', dpi=200)

        if oldnew == 'new':
            fig.suptitle('UTR variants')
            imagePath = os.path.join(savedir, 'UTR_vs_RBS.png'.format(seq.name))
            fig.savefig(imagePath, format='png', dpi=200)

        if oldnew == 'not7':
            fig.suptitle('UTR variants (minus T7His)')
            imagePath = os.path.join(savedir, 'UTR_vs_RBS_no_T7.png'.format(seq.name))
            fig.savefig(imagePath, format='png', dpi=200)

def wt_to_syns(seqs):
    seq_dict = dict((seq.name, seq) for seq in seqs)

    single_syns = ['syn6', 'syn13', 'syn15', 'syn27', 'syn30', 'syn36', 'syn42',
                   'syn6-30', 'syn[1-13]', 'syn[1-42]', 'syn[1-15]', 'syn',
                   'syn[1-30]', 'syn[1-36]', 'A-D7', 'B-D5', 'B-B3', 'wt[1-42]',
                   'wt[1-50]', 'wt[1-115]', 'wt36~42']

    #single_syns = ['syn6', 'syn13', 'syn15', 'syn27', 'syn30', 'syn36', 'syn42']

    #single_syns = seq_dict.keys()

    savedir = '/home/jorgsk/phdproject/5UTR/FriedrikeFolds/report/images'
    # For each ssyn, print wt-fold in -- and compare the distance in fold with
    # the distance in expression

    foldDist = []
    indDist = []
    rbsDist = []

    wt = seq_dict['wt']
    wt_probs = wt.binding_probabilities
    #wt_sd = wt_probs[atg-12:atg-2]

    for ssyn in single_syns:

        seq = seq_dict[ssyn]

        x = seq.binding_probabilities

        # get the 1-norm fold difference
        dFold = sum(abs(a-b) for (a,b) in zip(x, wt_probs))/len(x)
        foldDist.append(dFold)

        # get the difference in induced
        dInd = dist1(wt.induced_mean, seq.induced_mean)
        indDist.append(dInd)

        # get the RBS difference
        dRBS = dist1(wt.TNstartRate, seq.TNstartRate)
        rbsDist.append(dRBS)

        ###### plot both ssyn and wt (in --s)
        fig, ax = plt.subplots()
        xvals = range(1, len(x)+1)

        for vector, stripes, thikn in [(x, '-', 4), (wt_probs, '--', 2)]:

            if len(vector) != len(x):
                debug()

            # plot the whole strtch
            ax.plot(xvals, vector, linewidth=thikn, linestyle=stripes, c='b')

            atg = seq.TNstart
            sd_beg = max(0, atg-12+1)
            ax.plot(range(sd_beg+1, atg-2+1), vector[sd_beg:atg-2],
                    linewidth=thikn, linestyle=stripes , c='r')
            ax.plot(range(atg+1, atg+3+1), vector[atg:atg+3],
                    linewidth=thikn, linestyle=stripes, c='y')

        ax.set_ylim(-0.2, 1.2)
        ax.set_xlim(-1, len(x)+1)
        ticks = range(0,65,5)
        ticks[0] = 1
        ax.set_xticks(ticks)

        # set title
        ax.set_title(seq.name +'induced: ' + str(seq.induced_mean))

        ax.set_xlabel('nt position from transcription +1')
        ax.set_ylabel('Probability of nt single stranded', size=10)

        image = os.path.join(savedir, '{0}_and_wt.png'.format(seq.name))
        fig.set_size_inches(16,6)
        fig.savefig(image, format='png', dpi=200)


    #### plot the difference in RBS and the difference in fold
    # print a table
    print('')
    print('Name\nInduced-difference\tFold-difference\tRBS-difference')
    print('')
    for idiff, fdiff, rdiff, nme in zip(indDist, foldDist, rbsDist,
                                        single_syns):
        print nme
        I = format(idiff, '.4f')
        F = format(fdiff, '.4f')
        R = format(rdiff, '.4f')
        print('\t'.join([I, F, R]))
        print('')
    print('')

    fig, axes = plt.subplots(2, sharex=True, sharey=True)
    axes[0].scatter(foldDist, indDist)
    #axes[0].set_title('Fold vs induced distance')
    axes[0].set_ylabel('Difference in protein expression', size=10)
    axes[0].set_xlabel('Difference in over-all secondary structure', size=10)

    axes[1].scatter(rbsDist, indDist)
    #axes[1].set_title('RBS vs induced distance')
    axes[1].set_ylabel('Difference in protein expression', size=10)
    axes[1].set_xlabel('Difference in RBS-rates', size=10)
    fig.subplots_adjust(wspace=2)

    fig.suptitle('Fold of whole 5UTR explains difference in protein expression '
                 'better than RBS calculator', size=13)

    spea = scipy.stats.spearmanr(foldDist, indDist)
    pear = scipy.stats.pearsonr(foldDist, indDist)
    both = spea + pear

    headr = 'Spearman: ({0:.2f}, {1:.4f})  Pearson: ({2:.2f}, {3:.4f})'.format(*both)
    axes[0].set_title(headr)
    print headr, 'fold'

    spea = scipy.stats.pearsonr(rbsDist, indDist)
    pear = scipy.stats.spearmanr(rbsDist, indDist)
    both = spea + pear

    headr ='Spearman: ({0:.2f}, {1:.4f})  Pearson: ({2:.2f}, {3:.4f})'.format(*both)
    axes[1].set_title(headr)
    print headr, 'rbs'

    savedir = '/home/jorgsk/phdproject/5UTR/FriedrikeFolds/report/images'
    image = os.path.join(savedir, 'wt_singleSyn_correlation.png'.format(seq.name))
    fig.savefig(image, format='png', dpi=200)

def main():
    #dset = 'Fried'
    dset = 'celB'
    #plot = 'celb/figures_first_batch'
    #dset = 'designer_rna'
    plot = 'celb/figures_design'
    #dset = 'Growing'
    #dset = 'Growing_syn' # TODO, a non-His based extension

    # where the pickle is stored
    storeAdr = 'sequence_data/RBSdictstore'.format(dset)

    #ReadAndSaveData(storeAdr, dset)

    seqs = ReadData(storeAdr, dset)


    #RatePresenter(seqs)

    # NOTE TO SELF seqs hold instances of DNAClasses classes.
    #Correlator(seqs)

    # get probability of fold + 
    seqs = probability_reader(seqs, plot, downstream_seqlen=51)
    debug()

    # analysis part:
    # Correlate the difference in induced with the difference in fold
    #fold_similarity(seqs) # for 'new' only

    # Make RBS plots for old and new, showing some correlation for 'old' set,
    # buit no correlation for 'new' set.
    #rbs_correlation(seqs)

    # Fried compare the single-bp syn6, syn13, syn15, syn16, syn27, syn30,
    # syn36, and syn42, and syn6-30
    #wt_to_syns(seqs)

    # Growing: compare all growth variants to the WT for TN+50, TN+75, and
    # TN+100
    # T7 and His vs varying extensions

    return seqs

seqs = main()

# How to design sequences for testing your ideas.

# 1) If you decrease His, you can show that the fold becomes more and more
# different, and that this correlates with the difference in expresssion. The
# RBS calculator predicts identical initiation rates, but maybe the wt-syn fold
# is sufficiently different that you can argue that you are taking into account
# more than the RBS?
#
# 2) Alternatively you can create a whole new sequence. Take into account i)
# codon bias ii) restriction sites iii) His-like amino acid?

# Story: unintended consequences of the His-tag?

# 3) You remember (von hippel?) studies of the ribosome docking mechanism. Ideas
# about secondary structure. Maybe something more recent has come out?
