#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 20:19:58 2018

@author: sadra
"""

from ana_STL import STL_computation
from ana_STL import directed_distance

F=STL_computation(1,10)

p1=F.add_predicate(1,">",0.1)
p2=F.add_predicate(1,"<",0.3)

p3=F.add_predicate(1,">",0.1)
p4=F.add_predicate(1,"<",0.3)

f1=F.F(range(11),F.Conj([p1,p2]))
f2=F.F(range(8),F.Conj([p3,p4]))
r=directed_distance(F,f1,f2)
print(r)
