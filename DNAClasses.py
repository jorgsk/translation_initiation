# Ensure 4/5 is float and not 0
from __future__ import division

from IPython.Debugger import Tracer
debug = Tracer()

import Workhouse
import numpy as np

class DNAobject(object):
    """ Ancestor class for translation objects"""

    def __init__(self, name, sequence, TNstart):
        self.name = name
        self.sequence = sequence
        # initiation should be the starting index of the +1 nt

        # if TN start is auto, make TN start the first ATG after 0
        if TNstart == 'auto':
            try:
                self.TNstart = sequence.index('ATG')
                if self.TNstart == 0:
                    print 'Found ATG at 0, will ignore. Going for 2nd ATG'
                    self.TNstart = sequence[1:].index('ATG') + 1
            except ValueError:
                print 'ATG not found for {0}, check your sequences'.format(name)

        else:
            self.TNstart = TNstart
        # reading frame is 0, 1, or 2 depending on if (x-0)/3, (x-1)/3, or
        # (x-2)/3, respectively, is a whole number, where x is translation
        # start.
        self.frame = Workhouse.ReadingFrame(self.TNstart)

    def NucCount(self, start, stop):
        """ Count the types of nucleotides from start to stop. """
        self.Gs = self.sequence[start,stop].count('G')
        self.As = self.sequence[start,stop].count('A')
        self.Ts = self.sequence[start,stop].count('T')
        self.Cs = self.sequence[start,stop].count('C')

    def __repr__(self):
        """ Instead of seqs[0] returning DNAclass object bla bla bla, it will
        not return the name! :) So much nicer to work with! """
        return self.name

class TNobject(DNAobject):
    """ Class for translation objects with attributes that are relevant for
    translation initiation. Has DNAobject as ancestor. """

    def __init__(self, gene, utr, sequence, TNstart, induced):
        if utr != '':
            name = gene+'_'+utr
        else:
            name = gene

        DNAobject.__init__(self, name, sequence, TNstart)
        self.induced_list = induced # can be list of ints or int
        self.induced_mean = np.mean(induced)

        self.gene = gene
        self.utr = utr

    def addRBS(self, RBS):
        """ RBS brings proposed initiation sites and their rates. I order this
        into two lists and I sum the rates for each reading frame in framesums.
        """
        rates = [[],[],[]] # indices == reading frames
        starts = [[],[],[]] # indices == reading frames
        for result in RBS:
            # if there is more than 1 init rate, I choose only the first one??
            # What if there are more sites? Shit man. What was my reasoning?
            # ANSWER: in each result there is only one start site and one rate.
            # But there are many results! that's why you loop through them.
            if type(result.init_rate) is not float:
                if len(result.init_rate) == 1:
                    result.init_rate = result.init_rate[0]
                if len(result.start_site) == 1:
                    result.start_site = result.start_site[0]

            # start_site has been corrected to be correct for the original DNA
            # sequence
            if result.frame == 0:
                rates[0].append(result.init_rate)
                starts[0].append(result.start_site)

            if result.frame == 1:
                rates[1].append(result.init_rate)
                starts[1].append(result.start_site)

            if result.frame == 2:
                rates[2].append(result.init_rate)
                starts[2].append(result.start_site)

        self.RBSrates = rates
        self.RBSstarts = starts
        # the index in rates of what WE BELEIVE is the TNstart

        # rates will be eg. [[6,45,50],[16],[32]] and your TNstart is 32, and
        # your self.frame is 2. Thus you get index 0.
        myrateindex = starts[self.frame].index(self.TNstart)

        # myrate is the assumed rate of production from self.TNstart
        self.TNstartRate = rates[self.frame][myrateindex]

        framesums = [sum(rates[n]) for n in range(3)]

        #all_framesum is the sum of all rates
        self.frames_sum = sum(framesums)

        #my_framesum is the sum of rates in mystart's readingframe
        self.TNstartFrameRate = framesums[self.frame]

class TSobject(object):
    """ Class for transcription objects with attributes that are relevant for
    transcription initiation. Has DNAobject as ancestor. """

    def __init__(self, name, sequence, dataset, TNstart, PY):
        DNAobject.__init__(self,name,sequence,TNstart)

        self.PY = PY
        self.rna_dna1_10 = Energycalc.RNA_DNAenergy(self.sequence[:10])
        self.rna_dna1_15 = Energycalc.RNA_DNAenergy(self.sequence[:15])
        self.rna_dna1_20 = Energycalc.RNA_DNAenergy(self.sequence[:])
        # Redlisted sequences 
        self.redlist = []
        self.redlisted = False
        # Add self.promoter?
