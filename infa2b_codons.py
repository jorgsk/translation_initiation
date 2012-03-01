"""
Make synonymous mutations to infa2b and do some folding calculations

There are several approaches you can make

1) Average probability of not being in fold
2) Fold energy of first 30 nucleotides
3) ???

Ideally you'd have a way of fishing out deplorable hairpins
"""
import SequencePred

from IPython.Debugger import Tracer
debug = Tracer()

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

#1 Get the infa2b sequence
UTR5 = 'AACAGAAACAATAATAATGGAGTCATGAACAT' # wt UTR = 32 bp
infa2b = 'ATG' + 'GCTTGCGATCTGCCGCAGACCCATAGCCTGGGTAGCCGTCGCACCCTGATGCTGCTGGCACAGATG'

gene = UTR5 + infa2b

#2 Get codon tools of codons with more than 5% frequencyt
codonFrame, codon_subst = SequencePred.codon_structures(threshold=0.05)

#3 Get the codons You can change codon 2->8
myCodons = [infa2b[i:i+3] for i in range(0, len(infa2b), 3)]

# Run with a subset first to get the hang of it
#mutable_codons = myCodons[1:9]
mutable_codons = myCodons[1:3]

#3 Generate all unique combinations of the codons with itertools
maxlens = [len(codon_subst[c]) for c in mutable_codons]
# first generate all possible permutations up to 6
#Then filter those based on the max for each codon

mutant_indices = [i for i in codon_combos(maxlens)]

seqs = []

for indx in mutant_indices:
    this_seq = ''
    for cod_pos, syn_indx in enumerate(indx):
        codon = mutable_codons[cod_pos]
        this_seq += codon_subst[codon][syn_indx]

    seqs.append(this_seq)

# Now you have the library of all possible sequences!
# The next is simply to insert the new codons in genes and do da fold!
# With some luck you're done today with this work.

#4 Calculate the folding energy of original sequence, 30 nt into CDS
wt_probs = SequencePred.probability_array(seq=gene, fold_until=62, sample_name='wt')

#5 Compare the folding energy distribution to the original one
