"""
Purpose: make synonymous substitutions to show 1) you can/can't get the original
fold back using all codon substituttions and 2) you can't get the original fold
back with non-synonymous substitutions outside the coding regions

Consider the first 15 codons (45 nt)
"""

from Bio import SeqIO
from pandas import *
import random
import os
from subprocess import Popen, PIPE
import cPickle as pickle

from IPython.Debugger import Tracer
debug = Tracer()

#from matplotlib import pyplot as plt

def probability_array(seq, fold_until, sample_name):
    """
    Write sequence to file, run hybrid-ss, and return the probability array of
    the sequence.
    """

    here = os.path.dirname(os.path.realpath(__file__))
    sampled_dir = os.path.join(here, 'sequence_data/Fried/sampled_syns')

    outfilePath = os.path.join(sampled_dir, '{0}'.format(sample_name))
    outfileHandle = open(outfilePath, 'wb')

    outfileHandle.write('>{0}\n'.format(sample_name))
    outfileHandle.write(seq[:fold_until])
    outfileHandle.close()

    # run hybrid-ss with 100 samples
    cmd = ['hybrid-ss', '--suffix', 'DAT', '--tracebacks', '200',
           sample_name]

    # change directory to the sec dir
    os.chdir(sampled_dir)
    p = Popen(cmd, stderr=PIPE)
    p.wait()

    # get the single-stranded probabilities
    ss_array = []
    for line in open('.'.join([sample_name, '37', 'ext']), 'rb'):

        try:
            (nr, pSS, pSS_and_next_SS) = line.split('\t')
        # exception is that the last line has only 2 values 
        except ValueError:
            (nr, pSS) = line.split('\t')
            pSS.rstrip('\n')

        # don't append for the first line
        if nr == 'i':
            continue

        ss_array.append(float(pSS))

    # go back to whence you came!
    os.chdir(here)

    return ss_array

def codon_structures(threshold):
    """
    Return data structures useful for working with codons
    """

    # Get a codon table so you know what codons you can substitute with
    from Bio.Data import CodonTable
    bacterial_table = CodonTable.unambiguous_dna_by_id[11]

    # the codon <--> dna mappings
    letter_codon = bacterial_table.back_table
    codon_letter = bacterial_table.forward_table

    # For each codon, list all possible versions
    # Read the codon frequency table
    codon_names = {}
    codon_freqs = {}

    for line in open('codon_freqs', 'rb'):
        (name, triplet, d, d, rel_freq) = line.split()

        codon_names[triplet] = name
        codon_freqs[triplet] = float(rel_freq) # frequency for this amino acid

    # here you set the limit for hwo long the wt and syn should be. Rather, you
    # should read them from file and thus set the limit. Always set the limit up to
    # a codon triplet.

    # Experiment: make a data table of the codons!
    codonDict = {'letter': Series(codon_letter),
                 'name': Series(codon_names),
                 'freq': Series(codon_freqs)}

    codonFrame = DataFrame(codonDict)

    # make a subsitution table (each codon points to all synonymous ones)
    codon_subst = {}

    # only count codons used more than 10% of the time to avoid rare ones
    for codon, letter in codon_letter.items():
        codon_subst[codon] = []
        for c, l in codon_letter.items():
            if letter == l and codon_freqs[c] > threshold:
                codon_subst[codon].append(c)

    return codonFrame, codon_subst

def get_cindex(seqs, tn_start):
    """
    Return a list of indexes where wt and syn are the same; start from tn_start
    as provided.
    """

    wt_translated = seqs['wt'][10:]
    syn_translated = seqs['syn'][10:]

    codon_rangeWt = range(0, len(wt_translated),3)
    codon_rangeSyn = range(0, len(syn_translated),3)

    wt_codons = [wt_translated[pos:pos+3] for pos in codon_rangeWt]
    syn_codons = [syn_translated[pos:pos+3] for pos in codon_rangeSyn]

    change_indexes = []
    #combs = 0

    for sCodonNr, sCodon in enumerate(syn_codons):
        if sCodon != wt_codons[sCodonNr]:
            change_indexes.append(sCodonNr)
            #print sCodon, sCodonNr, wt_codons[sCodonNr]
            #print codon_subst[wt_codons[sCodonNr]]
            #print [codonFrame.ix[c]['letter'] for c in codon_subst[wt_codons[sCodonNr]]]
            #print [format(codon_freqs[c], '.2f') for c in codon_subst[wt_codons[sCodonNr]]]
            #combs *= len(codon_subst[wt_codons[sCodonNr]])

    #return change_indexes, combs
    return change_indexes, syn_codons

def codon_randomizer(syn_codons, change_indices, codon_subst):
    """
    Randomize codons in seq at positions in codon_index. Allowed substitutions
    given in codon_subst.
    """
    codons = []
    for cInd, cod in enumerate(syn_codons):
        # If you can change this one put random codon
        if cInd in change_indices:
            rCod = codon_subst[cod][random.randrange(0, len(codon_subst[cod]))]

            codons.append(rCod)
        # If not, include the one in syn already
        else:
            codons.append(cod)

    return codons

def randomize_syn(utr, wt_probs, rands, syn_codons, change_indices, codon_subst):

    keepers = []
    distances = []

    for sample_nr in range(rands):

        # randomly change 'syn' seqs
        rand_codons = codon_randomizer(syn_codons, change_indices, codon_subst)
        # add the utr because only CDS has been changed
        rand_seq = utr + ''.join(rand_codons)

        rand_probs = probability_array(seq=rand_seq, fold_until=60,
                                       sample_name=str(sample_nr))

        x = wt_probs
        y = rand_probs
        dFold = sum(abs(a-b) for (a,b) in zip(x, y))/len(x)

        if dFold < 0.05:
            keepers.append(y)

        distances.append((rand_seq, dFold))

    return keepers, distances

def main():
    # Get wt and syn seqs
    wtsyn = os.path.join(here, 'sequence_data/Fried/full_syn_wt.fasta')
    seqs_types = SeqIO.to_dict(SeqIO.parse(wtsyn, 'fasta'))
    # Get rid of the extra information
    seqs = dict((name, s.seq.tostring()) for (name, s) in seqs_types.items())

    # get the probabilitiesof WT
    wt_probs = probability_array(seq=seqs['wt'], fold_until=60, sample_name='wt')

    # get codon datstructures. codon_subst only contains codons used more than
    # 10% of the time.
    codonFrame, codon_subst = codon_structures(threshold=0.1)

    # get the codon-indexes at which syn and wt are the same; these are the
    # codons you may optimize
    change_indices, syn_codons = get_cindex(seqs, tn_start=10)

    utr = seqs['wt'][:10]

    r = 100000

    pickfile = os.path.join(here, 'sequence_data/pickled/distances_{0}'.format(r))

    # Only calculate if the file doesn't exist
    if not os.path.isfile(pickfile):

        keepers, distances = randomize_syn(utr, wt_probs, r, syn_codons,
                                           change_indices, codon_subst)
        pickle.dump((keepers, distances), open(pickfile, 'wb'), protocol=2)

    else:
        keepers, distances = pickle.load(open(pickfile, 'rb'))

    limits = [0.005, 0.01, 0.05, 0.1]

    counter = dict((lim, 0) for lim in limits)

    for val in distances:
        for l in limits:
            if val < l:
                counter[l] += 1

    print counter

    # output the mean, std and % that is 0.5%, 1%, 5% and 10% similar
    # output the percentage that is < 1% similar, 5% similar and 10% similar;
    # note 0.5% similar. Note that 50% similar is the same as random.

    # For those that are 1000*

    # note you can process 60000 seqs per hour. That means 1k hours do to all 60
    # mill possible combinations. That's 40 days.

    # RESULT 
    # Recall that 50% difference is 100% different
    # 95% of sequences are more than 10% different in fold
    # (55/10000) == 0.5% of sequences are 10% similar
    # (5/10000) == 0.05% of sequences are 1% similar
    # (2/10000) == 0.02% of sequences are 0.5% similar
    # (1/10000) == 0.01% of sequences fold the same

    # Thus it is possible to recover the fold down to +60 by changing the
    # remaining codons.


    debug()

here = os.path.dirname(os.path.realpath(__file__))
if __name__ == '__main__':
    main()


# OK, by now you've got the new sequences. Now you need to fold them and compare
# them to the WT fold; we're looking for something within 5/10 % of the fold.
# However, I don't see the reason for this because I don't find a correlation
# between the fold and the expresssion for the CDS variants. Clearly, it's only
# in certain places that the fold matters. Maybe I should simply correlate with
# the +30 fold? +40? Nobody knows.

# You are seeing that the probability structures change quite a bit with the
# sampling the sampling rate.
