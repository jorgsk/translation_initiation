"""
Apparently, this is the old version of Translation. The new version is only
using Hsalis calculator. This version used only 5' folding.
"""

import csv
import time

import os
import numpy
import rpy
import matplotlib.pyplot as plt
import random
import timeit
import copy
import itertools
from Bio.Seq import Seq
from scipy import stats
from operator import itemgetter
from decimal import *

science = 'science_standard'
biotech = 'full_sequences_standardModified'
Fried = 'xylSstandard'
#filename = biotech #selecting dataset (remember to change sliding window too!)
filename = Fried #selecting dataset (remember to change sliding window too!)
homedir = os.getcwd()

def statistics(useful, seqNrs,system):
    # Total nr. of changes = tchng, as decimal for approximation
    # expected changes at one position is 22*0.07= 1.54
    if seqNrs < 20: print 'Warning: less than 20 sequences'
    tch = Decimal(sum([sum(useful[:,i]) for i in range(5)]))
    nuclnr = len(useful)
    # Making numfreq in percentages of total change
    numFreq = [(Decimal(sum(useful[pos,:]))/tch).quantize(Decimal('0.01'),
                                                    rounding=ROUND_UP) for
               pos in range(nuclnr)]
    numFreq = [str(itm) for itm in numFreq]
    sigchange = numpy.zeros((nuclnr,5))
    sigchange = [[0,0,0,0,0] for nr in range(nuclnr)]
    if system == 'psrandomz':
        p = 0.07
        for j in range(nuclnr):
            for k in range(5):
                texas = rpy.r.binom_test(int(useful[j,k]),seqNrs,p,'t',0.97)
                truep = texas['estimate']['probability of success']
                pvalu = texas['p.value']
                confInt = texas['conf.int']
                if pvalu < 0.03 and truep > confInt[0] and truep < confInt[1]:
                    sigchange[j][k] = [pvalu]
                    #TODO: 2 on paper
    if system == 'randomz':
        p = 0.25
        for j in range(nuclnr):
            for k in range(5):
                texas = rpy.r.binom_test(int(useful[j,k]),seqNrs,p,'t',0.97)
                truep = texas['estimate']['probability of success']
                confInt = texas['conf.int']
                pvalu = texas['p.value']
                if pvalu < 0.03 and truep > confInt[0] and truep < confInt[1]:
                    sigchange[j][k] = [pvalu]
    # Take sigchange and recalculate the probabilities and pvalues and present
    # it nicely together with the expected number of changes.
    # Make a dictionary with the four bases as keys. Their elements are their
    # significant positions

    return sigchange

def TopDeltaG(seq, fro, to):
    # Choses the mean or median of the highest 10% binding energies 
    c = energyCalc(seq, fro, to)[1:]
    en = sorted(c[0],reverse=True)
#    crop = int(numpy.floor(len(en)*0.01)) # 3% of length of c[0]
#    meanie = numpy.mean(en[0:crop])
#    meadie = numpy.median(en[0:crop])
    meadie = en[0] # just taking the top...
    return meadie
    # choose one

def newscrute(randSeqs,limz):
    List = [['Ran'+str(ind), 0, 0, row] for ind, row in enumerate(randSeqs)]
    # list is reduced as criteria are not met!!
    for row in limz:
        strt = row[0][0]
        stpp = row[0][1]
        limt = row[1]
        Enrgz = energyCalc(List,strt,stpp)[1]
        newList = [seq for dni, seq in enumerate(List) if Enrgz[dni] > limt]
        List = newList
        if List ==[]: break
    return List

def scrutinizer(randSeq, fullIndxs, fullSeq, allintrvls='no'):
    # Removing duplicates
    l1 = len(randSeq)
    randSeq = list(set(randSeq))
    print l1 -len(randSeq)
    Sig, unused1, unused2 = sigFind(fullIndxs)
    # incorporate weighing according to coefficients, so that the intervals with
    # high coefficients are selected for higher than the low-intervals
    # maintain statistics of how many are selected out 
    if allintrvls=='no':
        intrvls = [[Sig[1][ind][0],Sig[1][ind][1]] for ind in range(len(Sig[1]))
                   if Sig[2][ind] > 0.45]
    # Selects from all intervals
    if allintrvls=='yes':
        intrvls = [[Sig[1][ind][0],Sig[1][ind][1]] for ind in
                   range(len(Sig[1]))]
    # Preparing the random sequences for energyCalc 
    randSeq = [('Ran'+str(ind), 0, 0, row) for ind, row in enumerate(randSeq)]
    limz = [[row, TopDeltaG(fullSeq,row[0],row[1])] for row in intrvls]
    losses = []
    # list is reduced as criteria are not met!!
    for row in limz:
        strt = row[0][0]
        stpp = row[0][1]
        limt = row[1]
        t1 = time.time()
        Enrgz = energyCalc(randSeq,strt,stpp)[1]
        print time.time() - t1
        newList = [seq for dni, seq in enumerate(randSeq) if Enrgz[dni] > limt]
        losses.append(len(randSeq) - len(newList))
        randSeq = newList
    # if-test to sort according to highest free energy and then
    # return just the UTR part of the sequences if list is nonempty
    modlist = []
    if randSeq != []:
        # Make list into dictionary so you can pluck values later
        tistel = dict()
        for uti in range(len(randSeq)):
            tistel[randSeq[uti][0]] = randSeq[uti][3]
        # Sorting list after highest free energy for [1,49]
        reList = energyCalc(randSeq,0,48)
        reList = [[reList[0][ind],reList[1][ind]] for ind in range(len(reList[0]))]
        reList = sorted(reList,key=itemgetter(1),reverse = True)
        randSeq = [[row[0], tistel[row[0]][7:29], row[1]] for row in reList]
        modlist = [[row[0], 0, 0, 'AACAUGU'+ row[1]] for row in randSeq]
    # List is only UTR + SG, but modlist contains AACAUGU-start as well.
    return randSeq, losses, modlist

def randomizer(n=10):
    # There are 11 nucleotides before SG and 7 after
    # Should the last A be mutated? It is not mutated in any if the samples. 

    # Totally random:
    begUTR = 'AACAUGU'
    wt = 'ACAAUAAUAAU'+'UCAUGAA' #without SG
    SG = 'GGAG'
    aftrUTR = 'CAU'+'AUGAGUAUUCAACAUUUCCGUGUCGCCCUUAUUCCCUUUUUU'
    takefrom = list('GAUC')
    ran = [''.join([random.choice(takefrom) for x in range(18)]) for y in
             range(n)]
    randomz = [begUTR+ran[ind][0:2]+'C'+ran[ind][3:11]+SG+ran[ind][11:13] + 'U'
       + ran[ind][14:18] + aftrUTR for ind in range(n)]
    #Pseudo-random
    house = [list('AUC'), list('GUC'), list('GAC'),list('GAU')]
    ranSeq = []
    for x in range(n):
        rans = ''
        for ind in range(18):
            if random.random() <=0.21:
                if wt[ind] == 'G': rans = rans + random.choice(house[0])
                if wt[ind] == 'A': rans = rans + random.choice(house[1])
                if wt[ind] == 'U': rans = rans + random.choice(house[2])
                if wt[ind] == 'C': rans = rans + random.choice(house[3])
            else: rans = rans + wt[ind]
        ranSeq.append(rans)
    psRandom = [begUTR+ranSeq[ind][0:2]+'C'+ranSeq[ind][3:11]+SG+ranSeq[ind][11:13]
                + 'U' + ranSeq[ind][14:18]+aftrUTR for ind in range(n)]
    randBino, psRandBino = binomFinder(n)
    return randomz, psRandom

def binomFinder(n):
    return ['ost','tae']
    #TODO
    # Should return how many possible sequences can be made with the randomizer.
    # Should also return how many can be returned by the company.

# sigFind returns two dataset from an indxs dictionary: a list of
# significant [Index,'(a,b)', <value>]-pairs, and a list of the nonsignificant pairs.
def sigFind(indxs,pval=0.05):
    spearPvals = dictSort(indxs['Spearman correlation p-value'], 0, rev=False)
    rnaTuples = dictSort(indxs['Spearman correlation'], 0, rev=False)
    Sig = [[ind,(row[0][0],row[0][1]), row[1]] for ind, row
           in enumerate(rnaTuples) if abs(spearPvals[ind][1]) <
           pval]
    NonSig = [[ind,(row[0][0],row[0][1]), row[1]] for ind, row
              in enumerate(rnaTuples) if abs(spearPvals[ind][1])
              >= pval]
    Both = [[ind,(row[0][0],row[0][1])] for ind, row in enumerate(rnaTuples)]
    # Interchange rows and columns
    Sigg = [[row[i] for row in Sig] for i in [0,1,2]]
    NonSigg  = [[row[i] for row in NonSig] for i in [0,1,2]]
    Bothh = [[row[i] for row in Both] for i in [0,1]]
    return [Sigg, NonSigg, Bothh]

# Plot all top results from sliding window analysis. TODO plot surrounding 
# regions of top values from 'data' or do statistics to make sure the value is not an outlier.
def corrPlot(data, indxs,seq, typ):
    tom = 5
    sig, nonsig, both = sigFind(indxs,0.01)
#    plt.figure(1)
#    plt.title('Rank correlation for intervals on mRNA for ' + typ)
#    xakse = [str(row).replace(' ','') for row in both[1]]
#    plt.xticks(both[0],xakse)
#    plt.scatter(sig[0],sig[2], c = 'r')
#    plt.scatter(nonsig[0],nonsig[2], c = 'b')
#    plt.show()
    for tell in range(len(sig[1])):
        strt = sig[1][tell][0]
        stop = sig[1][tell][1]
        c = energyCalc(seq, strt, stop)[1:]
        x = c[0]; y = c[1]
        corrWindowPlotter(x,y,c,tell,strt,stop,typ)

def newslide(pLow, pHigh):
    lowData, lowIndxs = slidingWin(genlist=pLow)
    highData, highIndxs = slidingWin(genlist=pHigh)
    data = {'lowData':lowData, 'highData':highData}
    indxs = {'lowData':lowIndxs, 'highData':highIndxs}
    return [data, indxs]

def corrWindowPlotter(x,y,c,nr,strt,stop,keyword):
        f = plt.figure(nr+2)
        f.text(.5, .95,'For '+keyword+' nucleotides '+str(strt)+' to '+str(stop),
               horizontalalignment='center')
        plt.subplot(211); plt.scatter(x,y, s=20, marker='o', c='w')
        plt.xticks(range(int(min(c[0])), int(max(c[0]))))
        plt.subplot(212); plt.hist(x)
        plt.xticks(range(int(min(c[0])), int(max(c[0]))))
        plt.savefig('images/Scatter_'+str(strt)+'_'+str(stop))
#        plt.show()
#        print('Mean: '+ str(numpy.mean(c[0]))+ ' Median: ' +
#        str(numpy.median(c[0])))

# datasplit() splits the sequence data according to penicillin and according to
# binding energy, and then does sequence statistics on the four sets.
def datasplit(b, indxs):
    # Splitting according to penicillin growth
    pthresh = 1001
    pHigh = [row for row in b if row[1] > pthresh]
    pLow  = [row for row in b if row[1] < pthresh]
    pLowStats = seq_analysis(pLow)
    pHighStats = seq_analysis(pHigh)
    splitData = {'pLow':pLow,'pHigh':pHigh}
    splitSeq_an = {'pLowStats':seq_analysis(pLow),'pHighStats': seq_analysis(pHigh)}
    splitSeq_an['enHigh']= dict(); splitSeq_an['enLow']= dict()
    splitData['enHigh']= dict(); splitData['enLow']= dict()
    # Then splitting according to energy values 
    sOrted = dictSort(indxs['Spearman correlation'], 0, rev=False)
    fromTo = [row for row in sOrted if abs(row[1]) > 0.5]
    adhoc = 50
    for ind,row in enumerate(fromTo):
        strt = row[0][0]
        stop = row[0][1]
        c = energyCalc(seq, strt, stop)[1:]
        x = c[0]; y = c[1]
        mark = [1]*len(c[0])
        corrWindowPlotter(x,y,c,adhoc,strt,stop,'split')
        ethresh = numpy.median(c[0])
        for ind, energ in enumerate(c[0]):
            if energ < ethresh: mark[ind] =0
        enHigh = [rowen for ind, rowen in enumerate(b) if mark[ind] == 1]
        enLow  = [rowen for ind, rowen in enumerate(b) if mark[ind] == 0]

        splitSeq_an['enHigh'][str(row[0])]=seq_analysis(enHigh,fro=strt,to=stop)
        splitSeq_an['enLow'][str(row[0])]=seq_analysis(enLow,fro=strt,to=stop)
        splitData['enHigh'][str(row[0])] = enHigh
        splitData['enLow'][str(row[0])] = enLow
        adhoc = adhoc +1
    return [splitSeq_an, splitData, fromTo]

# readfile() returns the sequence file as a listed list. The data file it reads
# from is already sorted in the high-first low-last fashion. Also removes
# whitespace from gene names
def readfile(name):
    f = open('sequence_data/Rahmi/'+name, 'rt')
    a = csv.reader(f, delimiter='\t')
    b = [[row[0].replace(' ','_'), int(row[1]), int(row[2]), row[3]] for row in a]
    f.close()
    return b

def energyCalc(genlist, fro, to):
    ename = 'ENERGY_TMP'
    tempFasta = open(homedir +'output_data/temp_output/' + ename, 'wt')
    for row in genlist:
        tempFasta.write('>' + row[0] + '\n' + row[-1][fro:to] + '\n')
    tempFasta.close()
    os.chdir(homedir + 'output_data/temp_output')
#    print len(genlist)
    lovely = [float(row.rstrip()) for row in os.popen('cat ' + ename + ' | '
                                                     ' hybrid-ss-min --stream '
                                                      '--tmin=30 --tmax=30', 'r')]
    c1 = [[genlist[ind][0], float(row), float(genlist[ind][1]), float(genlist[ind][2])] for
          ind, row in enumerate(lovely)] # combining free energy and data levels in one list
    c = [[row[i] for row in c1] for i in [0,1,2,3]]# Interchange rows and columns
    os.chdir(homedir)
    return c

def slidingWin(genlist, strt=0, stp=55, winmin=25, winmax=55):
    d = numpy.ones((4, stp-winmin+2-strt, winmax-winmin+1))*-1
    # k is the range of start positions, t the range of end positions.
    # TODO Fix the mismatch between gene starts at 1 and array at 1. Suggest you
    # adopt your own use in the seq_analysis function
    rpy.r.options(warn=-1) # suppress warning messages from R in the last loop
    for k in range(strt, stp-winmin+2):
        if k+winmax-1 > stp: limit = k+winmax-1 - stp
        else: limit =0
        for t in range(k+winmin-1, k+winmax - limit):
            c = energyCalc(genlist, k, t+1)[1:]
#            gradient, intercept, r_val, reg_pval, st_err = stats.linregress(c[0],c[1])
#            r_sqrd = r_val**2
            rankcorr = spearman(c[0], c[1])
            rank_pval = spearman(c[0], c[1],'pval')
            # In each matrix, start position is row number and adjusted stop is column nr
#            d[0][k-strt][t-k-winmin+1] = r_sqrd
#            d[1][k-strt][t-k-winmin+1] = reg_pval
            d[2][k-strt][t-k-winmin+1] = rankcorr
            d[3][k-strt][t-k-winmin+1] = rank_pval
            # d[statistic][strt + index = from nucleotide][winmin + index = to nucleotide]
            # TODO
            # Idea do statistics the correlations. Variance and such, because the
            # values you get are too uniform. Example: starting from 0, return the max,
            # min, mean, and variance of the correlations.
    indxs = dict()
#    indxs['R squared'] = dict()
#    indxs['Lin reg p-value'] = dict()
    indxs['Spearman correlation'] = dict()
    indxs['Spearman correlation p-value'] = dict()
    # Since each row in each d[i] gives a start-from nucleotide, by doing
    # argmax(axis=1) the column indices along those rows wih max values are
    # given. Thus I return for each start position the end position index that gave
    # highest values.
    indrankcorr = d[2].argmax(axis=1)
    for ind, val in enumerate(indrankcorr):
#        indxs['R squared'][strt+ind,strt+winmin+ind+val-1] = d[0][ind][val]
#        indxs['Lin reg p-value'][strt+ind,strt+winmin+ind+val-1] = d[1][ind][val]
        indxs['Spearman correlation'][strt+ind,strt+winmin+ind+val-1] = d[2][ind][val]
        indxs['Spearman correlation p-value'][strt+ind,strt+winmin+ind+val-1] = d[3][ind][val]

    os.chdir(homedir)
    return [d, indxs]

def dictSort(d, xx, rev=True):
    # Returns the dictionary sorted by keys (xx=0) or values (xx=1). Reverse
    # order (largest first) is default 
    return sorted(d.iteritems(), key=itemgetter(xx), reverse=rev)

# call a sequence analysis program. input: 0 to 20 and all of 5prime
# When 'to' is 'X' it means ON THE REAL SEQUENCE: [8,X]
# If system is 777_79, then you are using pseudorandom sequences
def seq_analysis(seq, fro=7, to=20,system='psrandomz'):
    if to > 29: to = 29
    if fro < 7: fro = 7
    # test if you need to reformat...
    if len(seq[0]) != 4:
           seq = [['sumth',0,0,row] for row in seq]

    # wt is the coding strand = non template = NT. Begins at site 1 (0 in list). The site of first
    # mutation is nt 8 (7 in list). Ends at position 29, three nucleotides before translational start.
    wt = 'AACATGTACAATAATAATGGAGTCATGAA'
    # reduce the wt segment to the piece you wish to compare with 
    wt = wt[fro:to]
    # cc is RNA sequence data for coding strand, bb is [induced, uninduced]
    cc = [Seq(row[-1][fro:to]).back_transcribe().tostring() for row in seq]
    bb = [[row[i] for row in seq] for i in [1,2]]
    # The number of nucleotides in the wt NT strand
    G = wt.count('G'); A = wt.count('A'); T = wt.count('T'); C = wt.count('C');
    # The change from the wt NT strand to the other strands. Doing ranking
    # analysis, so relative change should be equivalent to absolute.
    u = [[row.count('G')-G for row in cc], [row.count('A')-A for row in cc],
         [row.count('T')-T for row in cc], [row.count('C')-C for row in cc],
         [row.count('-')-0 for row in cc]]
    # Number of sequences and number of UTR positions to be analysed
    seqnum = len(cc); utrnum = len(wt)
    x = range(seqnum)
    # Finds the change in Purine = AG. Pyrimidine = TC from the wt
    # G, A, T, C, AG, TC, AT, AC, GT, GC
    D = {'G':[u[0][i] for i in x],
         'A':[u[1][i] for i in x],
         'T':[u[2][i] for i in x],
         'C':[u[3][i] for i in x],
         'purine':[u[0][i] + u[1][i] for i in x],
         'GT':[u[0][i] + u[2][i] for i in x],
         'GC':[u[0][i] + u[3][i] for i in x],
         'AT':[u[1][i] + u[2][i] for i in x],
         'AC':[u[1][i] + u[3][i] for i in x],
         'pyrimidine':[u[2][i] + u[3][i] for i in x],
         'deletion':[u[4][i] for i in x]
        }
    # Evaluates spearman against bb - the order of 'induced' in the original csv file
    # TODO: generate data with KNOWN percentage changes and TEST if your method
    # works or even makes sense.
    nList = ['G','A','T','C','deletion','purine','GT','GC','AT','AC','pyrimidine']
    1/0
    E = dict()
    for item in nList:
        E[item] = [spearman(bb[0], D[item]), spearman(bb[0], D[item],'pval')]
    E = dictSort(E,1)
    # F is [pos][seq][nucleotide]
    F = numpy.zeros((seqnum,utrnum,5),int)
    # H is [nucleotide][seq][pos] (more useful!)
    H = numpy.zeros((5,utrnum,seqnum),int)
    for ind, utr in enumerate(cc):
        for pos in range(len(utr)):
            if utr[pos] != wt[pos]:
                if utr[pos] == 'G': F[ind][pos][0] = 1
                if utr[pos] == 'A': F[ind][pos][1] = 1
                if utr[pos] == 'T': F[ind][pos][2] = 1
                if utr[pos] == 'C': F[ind][pos][3] = 1
                if utr[pos] == '-': F[ind][pos][4] = 1
    for k in range(5):
        for j in range(utrnum):
            for i in range(seqnum):
                H[k][j][i] = F[i][j][k]
    useful = [[sum(temo[pos,:]) for temo in H] for pos in range(utrnum)]
    # useful is [position][#G,#A,#T,#C,#-]], where #A is the change in As
    useful = numpy.array(useful)

    # tic gives correlation at each UTR position between substitution type and expression
    tic = dict()
    for ind, item in enumerate(nList[0:5]):
        tic[item] = [[pos, spearman(bb[0], H[ind][pos]),
                      spearman(bb[0], H[ind][pos],'pvalue')]
                     for
                     pos in range(utrnum) if spearman(bb[0], H[ind][pos], 'pval') < 0.05]

    deviant_stat = statistics(useful,len(cc),system)
    R = {'UTR ' + str((fro,to)):E, 'Changes at each pos': tic, 'significant '
    'changes': deviant_stat}
    return R, useful

def spearman(x, y, optn='rho'):
    if optn == 'rho': val = rpy.r.cor_test(x, y, method="spearman")['estimate']['rho']
    else: val = rpy.r.cor_test(x, y, method="spearman")['p.value']
    return val

def duplicateTester(seq):
    # This function identifies duplicate sequences in the sequence file in
    # "dupes" and returns a new sequence list in 'newSeq' that has replaced the
    # duplicates with mean expression values. Presumes there are no duplicate
    # names: these must be removed from the original file!
    # First investigate duplicate names HL95, HL87, and UP_C, (DONE) and go to the
    # other dictionary you make of seq and make sure that is OK
    #TODO by having 0 to 40 you're blocking for other datasets...
    tim = [row[-1][0:40] for row in seq]
    stein = []
    skipseqs = []
    for ind, row in enumerate(tim):
        n = tim[ind:].count(row) # nr. of seqs from ind and out
        # if sequence is in skipseqs it has already been counted!
        if n > 1 and (seq[ind][0] in skipseqs) == False:
            tam = copy.deepcopy(ind)
            while n > 1:
                i = tim.index(row,tam+1)
                tam = copy.deepcopy(i)
                n = n -1
                stein.append((ind,i))
                skipseqs.append(seq[i][0])
    dupes = [[seq[row[0]][0], (seq[row[0]][1],seq[row[0]][2]), seq[row[1]][0], (seq[row[1]][1],seq[row[1]][2])] for row in stein]
    newSeq = dict()
    for row in seq:
        newSeq[row[0]] = [row[1], row[2], row[-1]]
    for el in range(len(skipseqs)):
        del newSeq[skipseqs[el]]
    newVals = dict()
    cdupes = copy.deepcopy(dupes)
    namelist = [row[0] for row in dupes]
    for inx, row in enumerate(cdupes):
        n = namelist.count(row[0])
        t = copy.deepcopy(n)
        tempInd = [row[1][0], row[3][0]]
        tempUnInd = [row[1][1], row[3][1]]
        # remove the assimilated sequence name from the seq list
        while t > 1:
            tempInd.append(cdupes[inx+t-1][3][0])
            tempUnInd.append(cdupes[inx+t-1][3][1])
            del cdupes[inx+t-1]
            t = t-1
        # Mean
        newInd = sum(tempInd)/float(n+1)
        newUnInd = sum(tempUnInd)/float(n+1)
        # Unbiased standard deviation
        stdInd = numpy.sqrt(sum([(el-newInd)**2 for el in tempInd])/float(n))
        stdUnInd = numpy.sqrt(sum([(el-newUnInd)**2 for el in
                                   tempUnInd])/float(n))
        newVals[row[0]] = [newInd, stdInd, newUnInd, stdUnInd]
    for row in newVals:
        newSeq[row][0] = newVals[row][0]
        newSeq[row][1] = newVals[row][2]
    newSeq = [[row, newSeq[row][0], newSeq[row][1], newSeq[row][2]] for row in
              newSeq]
    # Sort the new sequence file first after uninduced and then after induced
    newSeq = sorted(newSeq,key=itemgetter(2),reverse=True)
    newSeq = sorted(newSeq,key=itemgetter(1),reverse=True)
    return dupes, newSeq, newVals, skipseqs

# Creates all possible sequences of length nr and analyses them against the best
# sequences of the original data set. Sig is the set of significant sequences
# determined by the sigFind() function.
def createAndSelect(Sig, fullSeq):
    nr = 7
    a = 4**nr
    print 'calculating '+ str(a) +' sequences'
    # Finding the sequence ranges where correlation coefficent is greater than
    # 0.45 (TODO justify this value by comparing the results of your set with
    # random DNA sequences -- when is there a difference in significance?)
    intrvls = [(Sig[1][ind][0],Sig[1][ind][1]) for ind in range(len(Sig[1]))
                   if Sig[2][ind] > 0.45]
    # Finding the top binding energy of each significant interval
    limz = [(row, TopDeltaG(fullSeq,row[0],row[1])) for row in intrvls]
    # Generating all possible permutations of GATC of length 'nr'
    begUTR = 'AACAUGU'
    wt = 'ACAAUAAUAAU'+'UCAUGAA' #without SG
    SG = 'GGAG'
    aftrUTR = 'CAU'+'AUGAGUAUUCAACAUUUCCGUGUCGCCCUUAUUCCCUUUUUU'
    x = begUTR
    os.chdir(homedir + 'output_data/temp_output')
    survivers = []
    # for each sequence, calculate the log(Z) score and keep if it is below the
    # limit for all intervals
    #TODO IDEA: if doing 1000 at a time, whenever one sequence is discovered
    #that does not fulfill the criteria -> jump to the next start. How to do
    # that is not obvious: maybe for i in range(1000) for i in range(len(limz))
    # and then do a return in the inner loop?
    y = 1001
    for seq in itertools.product('GATC',repeat=nr):
        if y > 1000:
            y=1
            tempstor = []
        x = ''.join(seq)
        x = begUTR+x[0:11]+SG+x[11:18]+aftrUTR
        tempstor.append(x)
        if y == 1000:
            survs = newscrute(tempstor,limz)
            if survs != []: survivers.append(survs)
        y = y+1
    return survivers

def surviverEquSeq(survivers,seq):
    seqsearch = [row[-1][7:29] for row in seq]
    survs = [row[1] for row in survivers]
    namelist = []
    for count in range(len(survivers)):
        if survs[count] in seqsearch:
            namelist.append(survs(count))
    return namelist

# a function that evaluates a set of sequences against the
# sequence result of the full analysis. This should be done in such a way that
# however the original data-set analysis changes, it should always be able to
# run a fresh data set against the data-set analysis to find out how this set is
# predicted to perform. Next step: include it in the data set, re-run the
# original analysis and see how they actually perform :)

# Short: takes new sequences + expression values and ranks among the existing
# dataset how well these perform according to calculated mRNA energies. Compares
# this to how well it _really_ performs. Should it be flexible enough to be used as
# a test for the mRNA folding energy model?


def shineDalgarno(seq):
    # Compare the sequences from 7 to 28 to the aSD sequence
    # TODO Save each 7to28 sequence in seq to a separate file and run the
    # hybrid-min program, extract the binding energies as well as the
    # nucleotides that bound to the mRNA an save this in a python list.
    os.chdir(homedir + 'output_data/antiSD')
    SDlist = [[row[0], row[3][7:28]] for row in seq]
    for SDrow in SDlist:
        storfil = open(SDrow[0], 'wt')
        storfil.write(SDrow[1])
        storfil.close()
        os.system('hybrid-min '+SDrow[0]+' anti') #calculating binding
        anneal = open(SDrow[0] +'-anti.asc','rt') #finding annealings
        binden = open(SDrow[0] +'-anti.dG','rt') #finding binding energy
        anne = csv.reader(anneal)
        bind = csv.reader(binden)
        for row in anne: print row
        for row in bind: print row
        print('\n')

    os.chdir(homedir)
#    ename = 'ENERGY_TMP'
#    tempFasta = open(homedir +'output_data/temp_output/' + ename, 'wt')
#    for row in genlist:
#        tempFasta.write('>' + row[0] + '\n' + row[-1][fro:to] + '\n')
#    tempFasta.close()
#    os.chdir(homedir + 'output_data/temp_output')
##    print len(genlist)
#    lovely = [float(row.rstrip()) for row in os.popen('cat ' + ename + ' | '
#                                                     ' hybrid-ss-min --stream '
#                                                      '--tmin=30 --tmax=30', 'r')]
#    c1 = [[genlist[ind][0], float(row), float(genlist[ind][1]), float(genlist[ind][2])] for
#          ind, row in enumerate(lovely)] # combining free energy and data levels in one list
#    c = [[row[i] for row in c1] for i in [0,1,2,3]]# Interchange rows and columns
#    os.chdir(homedir)
#    return c


# Create a plotting function with a loop that goes through all the relevant
# plots.

#t1 = time.time()

#dupes, newSeq, newVals, skipseqs = duplicateTester(seq)
#seq = newSeq
#shineDalgarno(seq)
#fullData, fullIndxs = slidingWin(genlist=seq)
#fullSeqStats, fullUseful = seq_analysis(seq)
#splitSeqStats, splitSeq, fromTo = datasplit(seq, fullIndxs)
#splitPenData, splitPenIndxs = newslide(splitSeq['pLow'], splitSeq['pHigh'])
#splitpenStats, splitPenUseful = seq_analysis(splitSeq['pHigh'])
#print time.time()-t1
# Find significant sequences from data-set
#significants, unused1, unused2 = sigFind(fullIndxs)
#survs = createAndSelect(fullIndxs,seq)
#t1 = time.time()

#randomz, psrandomz = randomizer(20)
#surviverSeqs, losses, modsurvs = scrutinizer(psrandomz, fullIndxs, seq,
#                                             allintrvls='no')
#surviverEquSeq(surviverseqs)
#wt = 'ACAAUAAUAAUGGAGUCAUGAA'
#print time.time()-t1
#surviverStats, surviverUseful = seq_analysis(modsurvs,system='psrandomz')
#accidents = surviverEquSeq(surviverSeqs,seq)
#print accidents

#corrPlot(splitPenData, splitPenIndxs['highData'], 'Penicillin High')
#corrPlot(splitPenData, splitPenIndxs['lowData'], 'Penicillin Low')
#corrPlot(fullData, fullIndxs, seq, 'Xyls')

#t = timeit.Timer('seq_analysis(seq)','from newmain import seq_analysis, seq')
#print(t.timeit(number=5)/5)

