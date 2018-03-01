#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:30:20 2018

@author: sadra
"""
from ana_STL import STL
from ana_STL import STL_computation

F=STL_computation(2,10)

p1=F.add_predicate(2,">",0.15)
p2=F.add_predicate(2,"<",0.6)

phi_1=F.Conj([p1,p2])
phi_3=F.G(range(2,11),F.Conj([p1,p2]))

F.best_signal(phi_3)
#F.worst_signal(phi_3)