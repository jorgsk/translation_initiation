""" Program to read csv text data and return it to th program the way you want
it. Prior to this reading, the data should be properly saved in a csv file using
Formatfixer.py
"""

from IPython.Debugger import Tracer
debug = Tracer()

import csv
import numpy as np
import os
# My own modules
import Workhouse
import DNAClasses
from Bio import SeqIO
from Bio import Seq

homedir = os.getcwd()

csv.register_dialect('jorgsk', delimiter='\t', quotechar="'",
                     quoting=csv.QUOTE_NONNUMERIC)

def AbortiveP():
    """ Import Hsu's abortive probabilities. Return dictionary with key
    promotername, which again contains a dict with keys 'Mean' and 'Std'.

    The data is not the same as in the 2006 paper.
    """
    fr = open('sequence_data/Hsu/abortiveProbabilities.csv','rt')
    reader = csv.reader(fr, delimiter='\t')
    temp = [line for line in reader]
    # Putting zeros where temp has empty string
    for row in temp[2:]:
        for ind, element in enumerate(row):
            if element == '':
                row[ind]=0
    tamp = np.array(temp)
    tamp = np.transpose(tamp)
    apDict = dict([(promoter,dict()) for promoter in temp[0][1::2]])
    # Even columns in tamp[1:] have the promoter name (even = 0,2,4,..). The
    # next carries the std and is reached by ind+1
    ost = tamp[1:]
    for ind, column in enumerate(ost):
        if not np.mod(ind,2): # This is an even column
            apDict[column[0]]['Mean'] = \
            np.array(Workhouse.StringOrFloat(column[2:-1]))/100
            apDict[column[0]]['Std'] = \
            np.array(Workhouse.StringOrFloat(ost[ind+1][2:-1]))/100
            apDict[column[0]]['FL+std'] = (float(column[-1]),
                                           float(ost[ind+1][-1]))
    return apDict

# Labels of CSV file ['Name','Sequence','PY','PYst','RPY','RPYst','RIF','RIFst','APR','MSAT','R']
def PYHsu(filename):
    """ Read Hsu csv-file with PY etc data. """
    f = open(homedir+'/sequence_data'+filename, 'rt')
    a = csv.reader(f, delimiter='\t')
    b = [[Workhouse.StringOrFloat(v) for v in row] for row in a]
    f.close()

    return b

def Rahmi104():
    """ Read Rahmi's 104 sequences and return as dictionary. """
#    Read the fixed Rahmi csv data and return a dict. 
    f = open(homedir+'/sequence_data/Rahmi/full_sequences_Fixed_for_all_FINAL','rt')
    # In FINAL I have changed some of the '-' values for induced/uninduced
    a = csv.reader(f, delimiter='\t')
    dicie = dict()
    for line in a:
        dicie[line[0]] = dict()
        dicie[line[0]]['Induced'] = float(line[1])
        dicie[line[0]]['Uninduced'] = float(line[2])
        sequence = line[3].replace('-','') # removing nucleotide deletion markers
        dicie[line[0]]['Sequence'] = sequence
    return dicie

def MinusTen():
    """ Reading the minus 10 elements of E. coli K12 as found on regulonDB. """
    fread = open('sequence_data/minusTens.csv','rt')
    min10s = []
    for row in fread.readlines():
        min10s.append(row.rstrip())
    return min10s

def NikaCombos(nafold=False):
    """ Read the promoter, UTRs, and genes (some with signal sequence) from a
    Fasta file into a dictionary. Then combine these elements to form all the
    combinations that Veronika is investigating in her experiments.

    Veronika's sequences have a -promoter of length 30, so the +1 transcription
    nucleotide is index 30. The translation start varies for the different ones,
    so this I don't know. pKO1 is not a UTR. It represents the wt+bla for this
    set of experiments. It should only be used to analyze the experiments
    themselves. Thus entries with f.ex 'celB_bla' should be discarded in this
    analysis. """

    # Reading the genes, utrs, and pKO1
    genes_utrs = open('sequence_data/Veronika/genes_utrs','rU')
    allseqs = dict()
    for record in SeqIO.parse(genes_utrs, "fasta"):
        allseqs[record.name] = str(record.seq)
    genes_utrs.close()

    # Reading the induced values for the different experiments.
    synthesis = open('sequence_data/Veronika/my_synthesis.csv','rU')
    synread = csv.reader(synthesis, dialect='jorgsk')
    header = synread.next()
    syndict = dict()
    for line in synread:
        rawline = line[1:]
        # excluding non-numbers from raw-line
        refline = [elem for elem in rawline if type(elem) in [int, float]]
        syndict[line[0]] = refline
    # NOTE I'm taking bla out of the output. If needed later, remember that it's
    # easies to include bla as a UTR since it occurs not with other utrs. Then
    # do a quickfix.
    genes = ['celB', 'gm', 'gm-omp', 'infa', 'infa-pelB', 'synluc']

    UTRs = ['LII-10', 'LII-12', 'LIV-1', 'LIV-2', 'LV-1', 'LV-2', 'wt']
    # Translation start relative to transcription start (python-index).
    TNstart = {'wt':32,'LIV-1':33,'LIV-2':33,'LV-1':32,'LV-2':32,'LII-10':32,'LII-12':30}
    # SD start and stop. Hand curated, like TNstart.
    SDstart = {'wt':16,'LIV-1':19,'LIV-2':19,'LV-1':16,'LV-2':17,'LII-10':16,'LII-12':16}
    SDstop = {'wt':21,'LIV-1':23,'LIV-2':22,'LV-1':21,'LV-2':21,'LII-10':21,'LII-12':21}

    promoter = allseqs['Promoter']
    TNobjects = []
    Na_arguments = []

    for gene in genes:
        for utr in UTRs:
            # NOTE removing the promoter part for now
            sequence = allseqs[utr] + allseqs[gene]
            name = gene+'_'+utr
            induced = syndict[name]
            TN_start = TNstart[utr]
            SD_start = SDstart[utr]
            SD_stop = SDstop[utr]
            tnobject = DNAClasses.TNobject(gene, utr, sequence, TN_start, induced)
            TNobjects.append(tnobject)
            # for Nafold objects I need the sequence from transcription start and on
            narg = [sequence, name, SD_start, SD_stop, TN_start, induced]
            # just returning the arguments, because sample data is needed too.
            Na_arguments.append(narg)

    if nafold:
        return Na_arguments
    else:
        return TNobjects

def Fried():
    mutantsf = os.path.join(homedir, 'sequence_data', 'Fried',
                            'utr_variants_new.csv')
    cdsf = os.path.join(homedir, 'sequence_data', 'Fried',
                        'xylS1000WT_for_adding_to_csv_sequence')

    cds = [l for l in open(cdsf, 'rb')][0].rsplit()[0]

    utr = 'ACCGTGAACC'

    # first do the original data
    seq_objects = []
    for line in open(mutantsf, 'rb'):
        splitline = line.split()
        seq_id = splitline[0].strip("\'")
        (induced, noninduced) = splitline[-2:]
        sequence = splitline[1:-2]
        seq = ''.join([s.strip("\'") for s in sequence])
        full_seq = utr + seq + cds[:100]
        # What is assumed TN start?
        TN_start = 'auto'
        induced = int(induced)
        d = ''

        tnObject = DNAClasses.TNobject(seq_id, d, full_seq, TN_start, induced)
        seq_objects.append(tnObject)

    # Then do the new data
    new_data = os.path.join(homedir, 'sequence_data', 'Fried',
                        'new_sequences')

    seqdict = {}
    for l in open(new_data, 'rb'):
        name, utr = l.split()
        seqdict[name] = utr

    for name, utr in seqdict.items():
        # don't store the wild type
        if name.startswith('xylS') or name.startswith('naturlig'):
            continue

        induced = -100
        TN_start = 'auto' # auto = first start codon

        full_seq_wt = utr + seqdict['xylSwt']
        tnObject = DNAClasses.TNobject(name, 'wt', full_seq_wt, TN_start, induced)
        seq_objects.append(tnObject)

        full_seq_syn = utr + seqdict['xylSsyn']
        tnObject = DNAClasses.TNobject(name, 'syn', full_seq_syn, TN_start, induced)
        seq_objects.append(tnObject)

    return seq_objects

# Only run if used as standalone program
if __name__ == '__main__':
    Fried()

