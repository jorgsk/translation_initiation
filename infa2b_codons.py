"""
Make synonymous mutations to infa2b and do some folding calculations

There are several approaches you can make

1) Average probability of not being in fold
2) Fold energy of first 30 nucleotides
3) ???

Ideally you'd have a way of fishing out deplorable hairpins
"""

import SequencePred

from subprocess import Popen, PIPE
from matplotlib import pyplot as plt

plt.ion()

from IPython.Debugger import Tracer
debug = Tracer()

class SequenceCandidate(object):
    """
    For 5' region candidates of infa2b.
    """

    # know about the codons so you can calculate the codon usage frequency
    _codonFrame, _codon_subst = SequencePred.codon_structures(threshold=0.05)

    def __init__(self, sequence, name, codon_freedom):
        self.seq = sequence
        self.name = name
        self.TNS = 32 # translation start site
        #self.TNS = 27 # translation start site
        # freedom = how many codons were changed after ATG
        self.codon_freedom = codon_freedom

        # the codon score is automaticall calculated as the average of the
        # codage usage bias from codon 1 to codon_freedom
        self.codon_score = self._get_codon_score()

    def _get_codon_score(self):
        """
        The codon score is the sum of the usage frequency of the codons that
        have been changed (the free codons)
        """
        # get the free codons only
        _free_codons = []

        for i in range(self.TNS+3, self.TNS+3 + self.codon_freedom*3, 3):
            _free_codons.append(self.seq[i:i+3])

        return sum([self._codonFrame.ix[c]['freq'] for c in _free_codons])


    def get_min_energy(self, TNplus=30):
        """
        Return the minimum folding energy. TNplus is the number of basepairs to
        fold upstream the translation start site.
        """

        cmd = ['hybrid-ss-min', '--quiet', self.seq[:self.TNS + TNplus]]

        # Return the minimum folding energy of the whole structure

        return float(Popen(cmd, stdout=PIPE).stdout.next().rstrip())

    def fold_probabilities(self):
        """
        Return an array of folding probabilities given an ensemble calculation
        """

    def ensemble_average_energy(self):
        """
        Sample the stochastic ensemble of structures and calculate the average
        energy
        """

    def relevant_hairpin_average_energy(self):
        """
        With the RNAStructure package you have the possibility to evaluate the
        energies of the specific hairpins.

        You don't fully understand everything that's got to do with bulges and
        etc, but all you'd need to do is sum all hairpin energies within +/- 15
        to 20 of the TN start site.

        One possibility is that the hairpin sum doesn't matter but the max of
        each hairpin influences the process.

        Damn, now you're deep in the mire aren't you.

        How to weight the different cases though?
        """

def codon_combos(max_lens):
    """
    Return the indices of condon combinations. Max_lens is a list of the
    maximum number of synonymous codon changes at each position in
    range(nr_codons)

    Check if all the permuted indices are less or equal to the max_lens indices
    For example if you have the following possibilityies [[ARG, VAR], [TAR],
    [SAR, MAR, HAR]] then max_lens should be [1, 0, 2], assuming 0-based indices

    This is a new version
    """
    ranges = [xrange(x) for x in max_lens]

    from itertools import product

    for perm in product(*ranges):
        yield perm

def get_fullseqs(UTR5, infa2b, maxlens, mutable_codons,
                 codon_subst, codon_freedom):
    """
    Use maxlens to get all possible codon combinations. Then construct all
    possible sequence combinations using the codon table.
    """
    seqs = []

    for indx in codon_combos(maxlens):
        this_seq = ''
        for cod_pos, syn_indx in enumerate(indx):
            codon = mutable_codons[cod_pos]
            this_seq += codon_subst[codon][syn_indx]

        seqs.append(this_seq)

    full_seqs = []

    for seq_nr, seq in enumerate(seqs):
        full_seq = UTR5 + 'ATG' + seq + infa2b[len(seq):]

        seqObj = SequenceCandidate(full_seq, 'seq_'+str(seq_nr), codon_freedom)

        # do some calculations on the seqobj
        seqObj.minEn = seqObj.get_min_energy()
        full_seqs.append(seqObj)

    return full_seqs

def sd_binding(seq):
    """
    Scan all hexamers for a match with the anti-shine dalgarno
    AGGAGG -> look for match to TCCTCC. This is rna-rna data. need to get some
    look-up tables.

    This ignores cost of initiation which is large.
    """
#Nearest-Neighbor Model, 1 M NaCl, pH 7a ∆G°37(kcal/mol)
    indiv = list(seq)
    neigh = [indiv[cnt] + indiv[cnt+1] for cnt in range(len(indiv)-1)]

    rna_rna = {
        'AA' : -0.93,
        'AU' : -1.10,
        'UA' : -1.33,
        'CU' : -2.08,
        'CA' : -2.11,
        'GU' : -2.24,
        'GA' : -2.35,
        'CG' : -2.36,
        'GG' : -3.26,
        'GC' : -3.42}

    energ = sum([rna_rna[nei] for nei in neigh])

    return energ
# Ignoring symmetty and stuff
#terminal AU            4.09 (0.22) 
#symmetry            0.45d (0.04)
            #0.43 0
#(self-complementary)
 #Numbers in parentheses are uncertainties for parameters. b Calcu-


def main():

#1 Get the sequecnes
    UTR5 = 'AACAGAAACAATAATAATGGAGTCATGAACAT' # wt UTR = 32 bp
    # NB! gene is without ATG
    #infa2b = 'GCTTGCGATCTGCCGCAGACCCATAGCCTGGGTAGCCGTCGCACCCTGATGCTGCTGGCACAGATG'
    # UPDATE: new, 5nt shorter 5UTR and infa3b without GCT as first codon (it
    # was for cloning purposes when used together with constructs)
    #infa2b = 'TGCGATCTGCCGCAGACCCATAGCCTGGGTAGCCGTCGCACCCTGATGCTGCTGGCACAGATG'
    # update2: this is the non-optimized infa2b
    infa2b = 'TGTGATCTGCCTCAAACCCACAGCCTGGGTAGCAGGAGGACCTTGATGCTCCTGGC'

    #2 Get codon tools of codons with more than 5% frequencyt
    codonFrame, codon_subst = SequencePred.codon_structures(threshold=0.05)

    #3 Get the codons You can change codon 1->8 (not counting ATG)
    myCodons = [infa2b[i:i+3] for i in range(0, len(infa2b), 3)]

    codon_freedom = 8

    # Run with a subset first to get the hang of it
    mutable_codons = myCodons[:codon_freedom]
    #mutable_codons = myCodons[1:16]

    #3 Generate all unique combinations of the codons with itertools
    maxlens = [len(codon_subst[c]) for c in mutable_codons]
    # first generate all possible permutations up to 6
    #Then filter those based on the max for each codon

    wtObj = SequenceCandidate(UTR5+'ATG'+infa2b, 'WT', codon_freedom)
    wtEn, wtScore = wtObj.get_min_energy(), wtObj.codon_score

    # the full UTR5 + ATG + synonymous mutation variant
    # full_seqs are objects that can be worked on
    seqObjects = get_fullseqs(UTR5, infa2b, maxlens,
                             mutable_codons, codon_subst, codon_freedom)

    ens = [sO.minEn for sO in seqObjects]
    rare_cods = [sO.codon_score for sO in seqObjects]

    #plt.hist(ens, bins=40)
    #plt.hist(rare_cods, bins=40)
    plt.scatter(ens, rare_cods)
    plt.scatter([wtEn], [wtScore], color='r')

    plt.xlabel('Energy')
    plt.ylabel('Codon Rarity Index')

    plt.title('WT infa2b has non-rare codons and tight fold compared to'\
             ' all possible synonymous mutation variants')

    fig, ax = plt.subplots()
    ax.hist(ens, bins=30)

    debug()

if __name__ == '__main__':
    main()


# Now you have the library of all possible sequences!
# With some luck you're done today with this work.

# What folding should I consider? What measures?

# The fold of the first 62 nucleotides (72, 82 -- make sure it's robust to
# length) and second the codon bias, which can be a sum (just rank the two and
# find the minimum rank)

#4 Calculate the folding energy of original sequence, 30 nt into CDS
# TODO the last way was good for Fried's seqs. Maybe this not so good in general.
#wt_probs = SequencePred.probability_array(seq=gene, fold_until=62, sample_name='wt')

#5 Compare the folding energy distribution to the original one
