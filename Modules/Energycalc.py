""" Calculate energy-levels: RNA/DNA and RNA/RNA """
from Bio.Seq import Seq
import copy
import numpy as np

# RNA-DNA duplexes from 2002 paper
NNRD = {'TT':-0.4,'TG':-1.6,'TC':-1.4,'TA':-1.0,'GT':-1.0,'GG':-1.5,'GC':-1.2,'GA':-0.9,'CT':-1.4,'CG':-2.4,'CC':-2.2,'CA':-1.5,'AT':-0.3,'AG':-0.8,'AC':-1.0,'AA':-0.2,'Initiation':1.0}
# EnLibRNA is for calculating expected random sequence values
EnLibRNA = copy.deepcopy(NNRD)
del(EnLibRNA['Initiation'])

# DNA-DNA duplexes from 2004 paper in the 5->3 direction
NNDD = {'AA':-1.0,'TT':-1.0,'AT':-0.88,'TA':-0.58,'CA':-1.45,'TG':-1.45,'GT':-1.44,'AC':-1.44,'CT':-1.28,'AG':-1.28,'GA':-1.3,'TC':-1.3,'CG':-2.17,'GC':-2.24,'GG':-1.84,'CC':-1.84,'TerminalAT':0.05,'Initiation':1.96,'Symmetry':0.43}
EnLibDNA = copy.deepcopy(NNDD)
del(EnLibDNA['Initiation'])

def RNA_DNAexpected(length):
    """ Returns expected RNA/DNA energy for random sequence of length 'length'. """
    mean = np.mean(EnLibRNA.values())
    # IT's flawed!! The nucleotides are not evenly distributed in the table!
    if length == 8:
        return -7.275
    if length == 9:
        return -8.447
    else:
        print "Length must be 8 or 9".
        1/0
    return mean*(length-1)

def DNA_DNAexpected(length):
    """ Returns expected DNA/DNA energy for random sequence of length 'length'. """
    mean = np.mean(EnLibDNA.values())
    return mean*(length-1)

def RNA_DNAenergy(sequence):
    """ Calculate the DNA/RNA binding energy of 'sequence'. Now skipping
    initiation cost. """
    indiv = list(sequence) # splitting sequence into individual letters
    neigh = [indiv[cnt] + indiv[cnt+1] for cnt in range(len(indiv)-1)]
    energ = sum([NNRD[nei] for nei in neigh])
#    energ = energ + NNRD['Initiation']
    return energ

def DNA_DNAenergy(sequence):
    """ Calculate the DNA/DNA binding energy of 'sequence'. """
    indiv = list(sequence)
    neigh = [indiv[cnt] + indiv[cnt+1] for cnt in range(len(indiv)-1)]
    energ = sum([NNDD[nei] for nei in neigh])
    # Terminal AT penalty
    if neigh[-1][-1] in ['A','T']:
        energ = energ + NNDD['TerminalAT']
    if neigh[0][0] in ['A','T']:
        energ = energ + NNDD['TerminalAT']
    # Symmetry penalty
    if Seq(sequence).complement().tostring() == sequence[::-1]:
        energ = energ + NNDD['Symmetry']
    energ = energ + NNDD['Initiation']
    return energ

def PhysicalRNA(sequence, msat, msatyes,maxlen):
    """ Return RNA/DNA binding energies for ITS-sequence. If msatyes == 'yes':
        considering max size of RNA/DNA hybrid as 'maxlen', adding first up to
        'maxlen', and then follows a window of size 'maxlen' until either 20 or
        msat. (pos1,pos2) = (0,6) means first 6 sequences
        (1,2,...,6) in the real sequence. Adds trailing E[energy] if
        msatyes=='yes'."""
    seqlist = list(sequence)
    enlist = []
    for pos1 in range(len(seqlist)-maxlen+1):
        if pos1 == 0:
            for pos2 in range(pos1+2,pos1+maxlen+1):
                subseq = ''.join([seqlist[tic] for tic in range(pos1,pos2)])
                enlist.append([RNA_DNAenergy(subseq),[pos1, pos2]])
        else:
            subseq = ''.join([seqlist[tic] for tic in range(pos1,pos1+maxlen)])
            enlist.append([RNA_DNAenergy(subseq),[pos1,pos1+maxlen]])
    # Adding expected RNA/DNA energy for maxlen bond +1 after msat values for fair
    # comparison between transcripts of different lengths.
    expEnrgy = RNA_DNAexpected(maxlen)
    if msatyes == 'yes':
        for zeror in range(int(msat)-1,19):
            enlist[zeror][0] = expEnrgy
    return enlist

promoter = 'ATAATAGATTC' # Last 11 nt's of promoter sequence (next is +1)
bubsize = 13
def PhysicalDNA(sequence):
    """ Return DNA/DNA hybrid energy of ITS-sequence using the (-11, 20)-element.
    (pos1,pos2) = (0,13) really means that the sequence 0,1,...,12 is joined. """
    seqlist = list(promoter+sequence)
    enlist = []
    pos1 = 0
    for pos2 in range(pos1+bubsize,len(seqlist)+1):
        subseq = ''.join([seqlist[tic] for tic in range(pos1,pos2)])
        enlist.append((DNA_DNAenergy(subseq),(pos1,pos2)))
    return enlist

