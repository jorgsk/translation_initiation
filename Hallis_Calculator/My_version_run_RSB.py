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

#You should have received a copy of the GNU General Public License
#along with Ribosome Binding Site Calculator.  If not, see <http://www.gnu.org/licenses/>.

#Copyright 2008-2009 Howard Salis

# This version by Jørgen reads from a specified file instead from command line.

from RBS_Calculator import RBS_Calculator
import sys, math

def MyRBS(sequence, start=0):
    """ Recieve sequences and return starting points and relative initiation
    frequencies and a kinetic folding parameter """

    if start == 0:
        start_range = [0, len(seq)]
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

    startsitenr = len(expr_list))
    returnlist = []
    for (expr,start_pos,ks) in zip(expr_list,start_pos_list,kinetic_score_list):
        returnlist.append([start_pos, expr, ks])
    return startsitenr, returnlist
