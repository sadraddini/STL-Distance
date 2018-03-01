#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 20:19:58 2018

@author: sadra
"""

from ana_STL import STL
from ana_STL import STL_computation

F=STL_computation(2,1)

p1=F.add_predicate(1,">",0.2)
p2=F.add_predicate(1,"<",0.6)

p3=F.add_predicate(1,">",0.25)
p4=F.add_predicate(1,"<",0.67)

phi_1=F.Conj([p1,p2])
phi_2=F.Conj([p3,p4])
r=F.directed_distance(phi_1,phi_2)
print(r)