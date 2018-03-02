#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 20:19:58 2018

@author: sadra
"""

from ana_STL import STL_computation
from ana_STL import directed_distance

F=STL_computation(2,10)

p1=F.add_predicate(1,">",0.2)
p2=F.add_predicate(1,"<",0.6)

p3=F.add_predicate(1,">",0.2)
p4=F.add_predicate(1,"<",0.6)

f1=F.G(range(10),F.Conj([p1,p2]))
f2=F.G(range(10),F.Conj([p3,p4]))
r=directed_distance(F,f2,f1)
print(r)
