"""
Make synonymous mutations to infa2b and do some folding calculations

There are several approaches you can make

1) Average probability of not being in fold
2) Fold energy of first 30 nucleotides
3) ???

Ideally you'd have a way of fishing out deplorable hairpins
"""

import SequencePred
import os
import numpy as np

from operator import attrgetter
from itertools import izip

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

    def __init__(self, sequence, name, codon_freedom, TNstart):
        self.seq = sequence
        self.name = name
        self.TNS = TNstart

        # how many codons should be changed after ATG
        self.codon_freedom = codon_freedom

        # the codon score is automatically calculated as the average of the
        # codage usage bias from codon 1 to codon_freedom
        self.codon_score = self._get_codon_score()

        # The amino acid sequence
        self.peptide_seq = self.get_peptide_seq()

        # prepare for energies
        self.minEns = []

    def __str__(self):

        return str((self.name, str(self.minEns)))

    def get_peptide_seq(self):
        """
        """
        # get the free codons only
        _free_codons = []

        for i in range(self.TNS+3, self.TNS+3 + self.codon_freedom*3, 3):
            _free_codons.append(self.seq[i:i+3])

        return ''.join([self._codonFrame.ix[c]['letter'] for c in _free_codons])

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


    def get_min_energy(self, temp, TNplus):
        """
        Return the minimum folding energy. TNplus is the number of basepairs to
        fold upstream the translation start site.
        """

        cmd = ['hybrid-ss-min', '-t', str(temp), '-T', str(temp),
               '--quiet', self.seq[:self.TNS + TNplus]]

        # Return the minimum folding energy of the whole structure

        return float(Popen(cmd, stdout=PIPE).stdout.next().rstrip())

    def get_detailed_energy(self, TNplus):
        """
        Calculate the minimum folding energy and get the energy of the
        nucleoltide bonds from +/- 10, 15 and 20 nucleotides around translation
        start.

        The idea is that only those chemical bonds matter. Parse the output
        files and
        """

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
                 codon_subst, codon_freedom, temp):
    """
    Use maxlens to get all possible codon combinations. Then construct all
    possible sequence combinations using the codon table.
    """
    seqs = []

    # 0 create all sequences from the codon combos
    for indx in codon_combos(maxlens):
        this_seq = ''
        for cod_pos, syn_indx in enumerate(indx):
            codon = mutable_codons[cod_pos]
            this_seq += codon_subst[codon][syn_indx]

        seqs.append(this_seq)

    full_seqs = {}

    # 1 create the objects
    for seq_nr, seq in enumerate(seqs):
        full_seq = UTR5 + 'ATG' + seq + infa2b[len(seq):]

        TNstart = 32
        seqObj = SequenceCandidate(full_seq, 'seq_'+str(seq_nr), codon_freedom,
                                  TNstart)

        full_seqs[seqObj.name] = seqObj

    # calculate folding energies for seqs of len 20, 30, 40, 50 downstream
    # translation start and calculate the folding energy in each case
    seq_names = [s.name for s in full_seqs.values()]

    # 2 calculate seq ens for various lengths
    for TNplus in [30]:
    #for TNplus in [20, 30, 40]:

        temp_file = 'sequence_data/temp_files/temp_ens_{0}'.format(TNplus)

        with open(temp_file, 'wb') as handle:

            # Oops: now we're not writing sequences in the 0 to N order any
            # more, but loop thorugh the dict sorted to fix that
            #debug()
            for seq_name, seqObj in full_seqs.items():

                handle.write('>{0}\n{1}\n'.format(seqObj.name,
                                                  seqObj.seq[:seqObj.TNS + TNplus]))

        # supply the seq_obj names so you get the right energy for the right
        # object!
        ens = get_ens(temp_file, temp, seq_names)

        for name, en in ens.iteritems():
            full_seqs[name].minEns.append((TNplus, float(en)))

    # 3 average for the 3 energies
    for seq_name, seq_obj in full_seqs.items():

        #just_ens = [t[1] for t in seq_obj.minEns]
        just_ens = [t[1] for t in seq_obj.minEns]

        seq_obj.meanEn = np.mean(just_ens)
        seq_obj.stdEn = np.std(just_ens)

        seq_obj.minEn30 = just_ens[0]

    return full_seqs

def sd_binding(seq):
    """
    Scan all hexamers for a match with the anti-shine dalgarno
    AGGAGG -> look for match to TCCTCC. This is rna-rna data. need to get some
    look-up tables.

    This ignores cost of initiation which is large.
    """
    #Nearest-Neighbor Model, 1 M NaCl, pH 7a deltaG 37(kcal/mol)
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

def get_ens(seq_path, temp, seq_names):
    """
    Rely on the fact that energies were written like this: 0,1,...,N-1,N
    """

    out = {}

    here = os.path.dirname(os.path.realpath(__file__))

    os.chdir(os.path.dirname(seq_path))

    new_path = os.path.join(here, seq_path)

    cmd = ['hybrid-ss-min', '-t', str(temp), '-T', str(temp), '-E', new_path]

    p = Popen(cmd, stdout=open('/dev/null', 'wb'))
    p.wait()

    handle = open(new_path + '.dG', 'rb')
    headr = handle.next()

    for seq_name, line in izip(seq_names, handle):

        out[seq_name] = line.split()[1]

    os.chdir(here)

    return out

def color_points(codon_freedom, TNstart, temp):
    """
    Return a set of folding energies for a set of sequences along with colors
    for plotting and a name.

    dict[name] = (energy, codon_score)
    """

    tested_files = 'celb/tested_seqs/thusfarly_testedseqs.txt'

    en_set = {}

    for line in open(tested_files, 'rb'):

        # only some have a specified color
        try:
            name, seq, color = line.split('\t')
            color = color.rstrip()
        except ValueError:
            name, seq = line.split()
            color = 'c'

        seqobj = SequenceCandidate(seq, name, codon_freedom, TNstart)
        energy = seqobj.get_min_energy(temp, TNplus=30)
        codonS = seqobj.codon_score

        en_set[name] = (float(energy), float(codonS), color)

    return en_set

def get_more_energies_for_paper():
    """
    Get energies and codon score for pelB random variants. Write to file.
    """

    # folding temperature
    temp = 30

    codon_freedom = 8
    TNstart = 32

    in_path = 'random_celb'

    out_handle = open('random_celb_output.txt', 'wb')

    seq_objects = []
    for line in open(in_path):
        name, seq = line.split()

        obj = SequenceCandidate(seq.upper(), name, codon_freedom, TNstart)
        en, cScore = obj.get_min_energy(temp, TNplus=30), obj.codon_score

        seq_objects.append((name, en, cScore))
        debug()

        out_handle.write('\t'.join([name, str(en), str(cScore)]) + '\n')

    out_handle.close()

    # make objects


def get_energies_for_paper():
    """
    Get energies and codon score for celB sequence variants. Write to file.
    """

    out_file = open('celb/energy_and_codonscore/energy_codon_score.txt', 'wb')

    out_file.write('Name\tMin_energy\tCodon_score\n')

    temp = 30
    TNstart = 32

    # the 5UTR plus atg
    UTR_atg = 'AACAGAAACAATAATAATGGAGTCATGAACATATG'

    # the first 82 nt of celB (we're only going to go up to 69
    celB_82 = 'CCCAGCATAAGCCCATTTGCCGGCAAGCCGGTCGATCCGGACCGTCTTGTCAATATCGACGCCCT'

    # infa2b
    infa2b = 'TGCGATCTGCCGCAGACCCATAGCCTGGGTAGCCGTCGCACCCTGATGCTGCTGGCACAGATG'

    # the variants have different lenghts of celB
    celb_lens = [3, 5, 8, 10, 15, 20, 23, 25, 30, 38, 69]

    codon_freedom = 8

    for cBlen in celb_lens:
        name = 'celB' +'_{0}'.format(cBlen)

        obj = SequenceCandidate(UTR_atg + celB_82[:cBlen] + infa2b, 'name',
                                codon_freedom, TNstart)


        en_score = obj.get_min_energy(temp, TNplus=30)
        codon_score = obj.codon_score

        if cBlen == 23:
            print cBlen
            print obj.seq
            print en_score
            print codon_score

        continue

        out_file.write('\t'.join([name, str(en_score), str(codon_score)]) + '\n')

    out_file.close()

def write_energies_again(tested_data, wtObj, pelbObj):

    out_file = open('celb/energy_and_codonscore/energy_codon_score_THE_REST.txt', 'wb')

    out_file.write('Name\tMin_energy\tCodon_score\n')

    temp=30

    for name, (en_score, codon_score, color) in tested_data.items():

        out_file.write('\t'.join([name, str(en_score), str(codon_score)]) + '\n')

    out_file.close()


def main():

    # some more energies
    get_energies_for_paper()
    # random cel b energies and codon score
    #get_more_energies_for_paper()
    debug()

    #1 Get the sequecnes
    UTR5 = 'AACAGAAACAATAATAATGGAGTCATGAACAT' # wt UTR = 32 bp
    # NB! gene is without ATG
    # NB2: This is the 'cloning' version
    #infa2b = 'GCTTGCGATCTGCCGCAGACCCATAGCCTGGGTAGCCGTCGCACCCTGATGCTGCTGGCACAGATG'
    # UPDATE: new, 5nt shorter 5UTR and infa3b without GCT as first codon (it
    # was for cloning purposes when used together with constructs)
    infa2b = 'TGCGATCTGCCGCAGACCCATAGCCTGGGTAGCCGTCGCACCCTGATGCTGCTGGCACAGATG'

    #2 Get codon tools of codons with more than 5% frequency
    codonFrame, codon_subst = SequencePred.codon_structures(threshold=0.05)

    #3 Get the codons You can change codon 1->8 (not counting ATG)
    myCodons = [infa2b[i:i+3] for i in range(0, len(infa2b), 3)]

    codon_freedom = 8

    # folding temperature
    temp = 30

    # Run with a subset first to get the hang of it
    mutable_codons = myCodons[:codon_freedom]
    #mutable_codons = myCodons[1:16]

    #3 Generate all unique combinations of the codons with itertools
    maxlens = [len(codon_subst[c]) for c in mutable_codons]
    # first generate all possible permutations up to 6
    #Then filter those based on the max for each codon

    # wild type ifn-a2b
    TNstart = 32
    wtObj = SequenceCandidate(UTR5+'ATG'+infa2b, 'WT', codon_freedom, TNstart)
    wtEn, wtScore = wtObj.get_min_energy(temp, TNplus=30), wtObj.codon_score

    pelB = 'AAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGTGCGATCTGCCG'
    # the pelB variant
    TNstart = 32
    pelbObj = SequenceCandidate(UTR5+'ATG'+pelB, 'WT', codon_freedom, TNstart)
    pelbEn, pelbScore = pelbObj.get_min_energy(temp, TNplus=30), pelbObj.codon_score

    #pepseq = wtObj.get_peptide_seq()

    # the full UTR5 + ATG + synonymous mutation variant
    # full_seqs are objects that can be worked on
    # NOTE now you have 3 energies instead of 1
    seqObjects = get_fullseqs(UTR5, infa2b, maxlens, mutable_codons,
                              codon_subst, codon_freedom, temp)

    ens = [sO.minEn30 for sO in seqObjects.values()]
    rare_cods = [sO.codon_score for sO in seqObjects.values()]

    # to check that you get the same peptide sequence
    #pepseqs = set([d.get_peptide_seq() for d in seqObjects.values()])

    # Get some color points for the extra tested data
    #tested_data = color_points(codon_freedom, TNstart, temp)

    # Write more energy+codon_usage for the paper
    #write_energies_again(tested_data, wtObj, pelbObj)

    #keepers = ['AAA1/2', 'AAA2/2', 'HP1', 'HP2', 'HP3', 'NativeCodons']
    # screen for the ones you want
    #for k, v in tested_data.items():
        #if k not in keepers:
            #tested_data.pop(k)

    # get the data that was selected based on energy and codon usage
    tested2 = 'celb/sequenced_stuff_for_plotting.txt'
    tested_data2 = get_tested(tested2, codon_freedom, TNstart, UTR5, infa2b)

    plot_results(ens, rare_cods, wtEn, wtScore, extra_data=tested_data2)

    # sequence_selecter
    #seq_select(seqObjects, codon_freedom)

def get_tested(path_to_seqs, codon_freedom, TNstart, UTR5, infa2b):

    en_set = {}
    temp = 30

    for line in open(path_to_seqs, 'rb'):
        if line.startswith('S'):
            name, seq, color = line.split()

        # add UTR and downstream gene
        full_seq = UTR5 + seq + infa2b[24:]

        seqobj = SequenceCandidate(full_seq, name, codon_freedom, TNstart)
        energy = seqobj.get_min_energy(temp, TNplus=30)
        codonS = seqobj.codon_score

        en_set[name] = (float(energy), float(codonS), color)

    return en_set


def seq_select(seqObjs, codon_freedom):
    """
    Select some of the sequences for expression

    New: select first on minimal mean energy, second on minimal std
    """

    forbidden1 = 'CATATG' #NdeI
    forbidden2 = 'CCATGG' #NcoI

    atr = 'codon_score'

    mins = sorted(seqObjs.values(), key=attrgetter(atr))

    counter = 0
    for a in mins:
        checkseq = a.seq[a.TNS:]
        if forbidden1 in checkseq or forbidden2 in checkseq:
            counter += 1

    # write last 50 to fasta file with name with folding energy and codon value
    outhandle_low = open('infa2b/for_testing_low_{0}.fasta'.format(atr), 'wb')
    for seq_nr, low in enumerate(mins[-100:]):
        #checkseq = low.seq[low.TNS:low.TNS + codon_freedom*3 + 3]
        checkseq = low.seq[:low.TNS + codon_freedom*3 + 3]
        #if forbidden1 in checkseq or forbidden2 in checkseq:
        if forbidden1 in checkseq[low.TNS:] or forbidden2 in checkseq[low.TNS:]:
            continue
        else:
            #ens = format(low.meanEn, '.3f') + '(' + ','.join([str(t[1]) for t in low.minEns]) +')'
            ens = str(low.minEn30)
            title = '>Seq_{0}-Ens:{1},-Co:{2}'.format(seq_nr, ens, low.codon_score)
            outhandle_low.write('{0}\n{1}\n'.format(title, checkseq))

    outhandle_low.close()

    outhandle_high = open('infa2b/for_testing_high_{0}.fasta'.format(atr), 'wb')

    for seq_nr, low in enumerate(mins[:100]):
        #checkseq = low.seq[low.TNS:low.TNS + codon_freedom*3 + 3]
        checkseq = low.seq[:low.TNS + codon_freedom*3 + 3]
        #if forbidden1 in checkseq or forbidden2 in checkseq:
        if forbidden1 in checkseq[low.TNS:] or forbidden2 in checkseq[low.TNS:]:
            continue
        else:
            #ens = format(low.meanEn, '.3f') + '(' + ','.join([str(t[1]) for t in low.minEns]) +')'
            ens = str(low.minEn30)
            title = '>Seq_{0}-Ens:{1},-Co:{2}'.format(seq_nr, ens, low.codon_score)
            outhandle_high.write('{0}\n{1}\n'.format(title, checkseq))

    outhandle_high.close()

def plot_results(ens, rare_cods, wtEn, wtScore, extra_data=False):

    fig0, ax0 = plt.subplots()
    ax0.scatter(ens, rare_cods)

    # Add the WT with a red dot
    ax0.scatter([wtEn], [wtScore], color='DarkRed')

    #ax0.annotate('WT',
            #xy=(wtEn, wtScore), xytext=(-20,20),
            #textcoords='offset points', ha='right', va='bottom',
            #bbox=dict(boxstyle='round, pad=0.5', fc='DarkRed', alpha=0.9),
            #arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'),
                #size=12)

    # add tested data, with little labels
    for name, (energy, codon_score, colr) in extra_data.items():
        # add the point
        ax0.scatter([energy], [codon_score], color=colr)

        #add the label
        ax0.annotate(name,
                xy=(energy, codon_score), xytext=(-20,-20),
                textcoords='offset points', ha='right', va='bottom',
                bbox=dict(boxstyle='round, pad=0.5', fc=colr, alpha=0.9),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'),
                    size=12)

    #ax0.set_xlabel('Folding energy of ribosome binding site (kcal/mol)', size=22)
    #ax0.set_ylabel('Codon usage index', size=22)

    plt.setp(ax0.get_xticklabels(), fontsize=14)

    ax0.set_ylim(1.75, 4.5)
    ax0.set_xlim(-20, -2)

    fig0.set_dpi(300)

    fig0.set_size_inches(18,12)

    #fig, ax = plt.subplots()
    #ax.hist(ens, bins=50)
    fig_dir = '/home/jorgsk/phdproject/5UTR/MutationCorrelation/celb/scatter_fig'

    for formt in ['pdf', 'eps', 'png', 'tiff']:

        name = 'for_poster.' + formt

        fig0.savefig(os.path.join(fig_dir, name), transparent=True,
                     format=formt, dpi=300)


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
