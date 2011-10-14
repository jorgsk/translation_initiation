#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#This file is part of the Ribosome Binding Site Calculator.

#The Ribosome Binding Site Calculator is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#The Ribosome Binding Site Calculator is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#mYou should have received a copy of the GNU General Public License
#along with Ribosome Binding Site Calculator.  If not, see <http://www.gnu.org/licenses/>.

#Copyright 2008-2009 Howard Salis

# This version by Jørgen reads from a specified file instead from command line.

from RBS_Calculator import RBS_Calculator
import sys, math

class Result(object):
    """ Class of result objects that have a relative initiation site, relative
    translation rate, and kinetic parameter I don't understand. """
    def __init__(self, RBSstart_pos, expr, ks, subDNA_frame, subDNA_start):
        # startframe is of the mRNA relative to the full length DNA.
        # relative_frame is the frame of translation start within the mRNA.
        # Use UTRstart to get return values relative to transcription +1

        # start_site is relative to full length DNA by adding start of
        # invesitagted subsequence to where RBS estimates that translation
        # starts on this subsequence
        self.start_site = subDNA_start + RBSstart_pos

        # Finding the frame position that RBSstart is in
        RBS_frame = [(RBSstart_pos-n)%3 for n in range(3)]
        # scale TLstart frames found by RBS by rbs_frame
        relative_frame = RBS_frame.index(0)
        # relative to full length frame
        self.frame = (subDNA_frame+relative_frame)%3

        self.init_rate = expr
        self.kinetic = ks

def MyRBS(sequence, subDNA_frame, subDNA_start, start=0):
    """ Recieve sequences and return starting points and relative initiation
    frequencies and a kinetic folding parameter. I don't know how to use the
    kinetic folding parameter. """

    # You can force the calculator to start at a specific start codon
    if start == 0:
        start_range = [0, len(sequence)]
    else:
        start_range = [int(start), int(start)+1]

    name = "no name"

    #Create instance of RBS Calculator
    calcObj = RBS_Calculator(sequence, start_range, name)
    calcObj.calc_dG()

    dG_total_list = calcObj.dG_total_list[:]
    start_pos_list = calcObj.start_pos_list[:]
    kinetic_score_list = calcObj.kinetic_score_list[:]

    expr_list = []
    for dG in dG_total_list:
        expr_list.append(calcObj.K * math.exp(-dG/calcObj.RT_eff))

    returnlist = []
    for (expr, RBSstart_pos, ks) in zip(expr_list, start_pos_list, kinetic_score_list):
        result_object = Result(RBSstart_pos, expr, ks, subDNA_frame, subDNA_start)
        returnlist.append(result_object)
    # Return the number of start sites, and a list of start positions, relative
    # initiation rates, and the kinetic parameter.
    return returnlist
