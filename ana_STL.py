#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:17:12 2018

@author: sadra
"""

"""
    The task here is to develop an automated framework for converting STL 
    formulae to MILP Constraints. 
    Ideally I want to do parsing.
"""

from copy import copy,deepcopy

from gurobipy import *
import random

bigM=20

class STL:
    def __init__(self,name="formula"):
        self.name=name
        self.composition=""
        self.type="formula" #it can be predicate
        self.children=[]
    
    def __repr__(self):
        return self.name
    
class STL_computation:
    def __init__(self,n,T=0):
        self.n=n   
        self.T=T+1
        self.model=Model("STL-spec")
        self.z={}
        self.y={}
        self.subformulas=[]
        self.Y=range(1,n+1)
        self.construct()
        self.n_predicates=0
    
    def __deepcopy__(self,memo):
        print("*"*80,"\n We are doing a deepcopy",str(memo))
        new=STL_computation(self.n,self.T)        
        new.model=self.model.copy()
        new.subformulas=self.subformulas
        return new

    
    def construct(self):
        self.rho=self.model.addVar(lb=0,ub=10,name="r")
        self.flag=1 # 1 for relaxation
                    # -1 for tightening
                    # 0 for holding constant
        for t in range(0,self.T):
            for y in self.Y:
                self.y[y,t]=self.model.addVar(lb=0,ub=1,name=str(y)+","+str(t))  
        self.model.update()

    def add_predicate(self,y,sign,c):
        self.n_predicates+=1
        f=STL("p"+str(self.n_predicates)+" over "+"y"+str(y))
        f.type="predicate"
        f.children=[f]
        
        self.subformulas.append(f)
        for t in range(self.T):
            self.z[f,t]=self.model.addVar(vtype=GRB.BINARY,name=str(f)+","+str(t))
        self.model.update()
        f.y=y
        f.sign=sign
        f.c=c
        return f
    
    def predicate_constraints(self,predicate,flag):
        f=predicate
        y_signal=predicate.y
        sign=predicate.sign
        c=predicate.c
        if sign=="<":
            for t in range(self.T):
                self.model.addConstr(self.y[y_signal,t]<=c+bigM-bigM*self.z[f,t]+flag*self.rho)
                self.model.addConstr(self.y[y_signal,t]>=c-bigM*self.z[f,t]+flag*self.rho)
        elif sign==">":
            for t in range(self.T):
                self.model.addConstr(self.y[y_signal,t]>=c-bigM+bigM*self.z[f,t]-flag*self.rho)
                self.model.addConstr(self.y[y_signal,t]<=c+bigM*self.z[f,t]-flag*self.rho)
        else:
            raise("Error: invalid form of predicate. Use only '<' or '>' ")
    
    def add_formulas(self,list_of_formulas):
        self.subformulas.extend(list_of_formulas)
        for f in list_of_formulas:
            for t in range(self.T):
                self.z[f,t]=self.model.addVar(lb=0,ub=1)
        self.model.update()
            
    def Conj(self,list_of_formulas):
        for f_in in list_of_formulas:
            if f_in not in self.subformulas:
                raise("Error:",f_in," one of the formulas not defined before")
        f_out=STL()
        self.subformulas.append(f_out)
        f_out.name=("phi_%d"%(len(self.subformulas)-self.n))
        for t in range(self.T):
            self.z[f_out,t]=self.model.addVar(lb=0,ub=1,name=str(f_out)+","+str(t))
        self.model.update()
        for t in range(self.T):
            s=LinExpr()
            for f_in in list_of_formulas:
                self.model.addConstr(self.z[f_out,t]<=self.z[f_in,t])
                s.add(self.z[f_in,t])
            self.model.addConstr(self.z[f_out,t]>=s-len(list_of_formulas)+1)
        f_out.composition="Conjunction of "+str(list_of_formulas)
        for formula in list_of_formulas:
            f_out.children.extend(formula.children)
        return f_out
            
    def Disj(self,list_of_formulas):
        for f_in in list_of_formulas:
            if f_in not in self.subformulas:
                raise("Error:",f_in," one of the formulas not defined before")
        f_out=STL()        
        self.subformulas.append(f_out)
        f_out.name=("phi_%d"%(len(self.subformulas)-self.n))
        for t in range(self.T):
            self.z[f_out,t]=self.model.addVar(lb=0,ub=1,name=str(f_out)+","+str(t))
        self.model.update()
        for t in range(self.T):
            s=LinExpr()
            for f_in in list_of_formulas:
                self.model.addConstr(self.z[f_out,t]>=self.z[f_in,t])
                s.add(self.z[f_in,t])
            self.model.addConstr(self.z[f_out,t]<=s)
        f_out.composition="Disjunction of "+str(list_of_formulas)
        for formula in list_of_formulas:
            f_out.children.extend(formula.children)
        return f_out            
        
    def G(self,I,f): 
        if f not in self.subformulas:
            raise("Error:",f," one of the formulas not defined before")
        f_out=STL()
        self.subformulas.append(f_out)
        f_out.name=("phi_%d"%(len(self.subformulas)-self.n))
        for t in range(self.T):
            self.z[f_out,t]=self.model.addVar(lb=0,ub=1,name=str(f_out)+","+str(t))
        self.model.update()
        for t in range(self.T-I[-1]):
            s=LinExpr()
            for tau in I:
                self.model.addConstr(self.z[f_out,t]<=self.z[f,t+tau])
                s.add(self.z[f,t+tau])
            self.model.addConstr(self.z[f_out,t]>=s-len(I)+1)
        f_out.composition="Always "+str(I)+" "+f.name
        f_out.children.extend(f.children)
        return f_out        
        
    def F(self,I,f):
        if f not in self.subformulas:
            raise("Error:",f," one of the formulas not defined before")
        f_out=STL()        
        self.subformulas.append(f_out)
        f_out.name=("phi_%d"%(len(self.subformulas)-self.n))
        for t in range(self.T):
            self.z[f_out,t]=self.model.addVar(lb=0,ub=1,name=str(f_out)+","+str(t))
        self.model.update()
        for t in range(self.T-I[-1]):
            s=LinExpr()
            for tau in I:
                self.model.addConstr(self.z[f_out,t]>=self.z[f,t+tau])
                s.add(self.z[f,t+tau])
            self.model.addConstr(self.z[f_out,t]<=s)
        f_out.composition="Eventually "+str(I)+" "+f.name
        f_out.children.extend(f.children)
        return f_out
    
    def best_signal(self,f):
        if f not in self.subformulas:
            raise("Error:",f," one of the formulas not defined before")
        for p in f.children:
            print(p)
            if p.type=="predicate":
                print(p)
                self.predicate_constraints(p,-1)
        self.model.addConstr(self.z[f,0]==1)
        self.model.setObjective(self.rho)
        self.model.optimize()
        for y,val in self.y.items():
            print(y,val.X)


    def worst_signal(self,f):
        if f not in self.subformulas:
            raise("Error:",f," one of the formulas not defined before")
        for p in f.children:
            print(p)
            if p.type=="predicate":
                print(p)
                self.predicate_constraints(p,-1)
        self.model.addConstr(self.z[f,0]==0)
        self.model.setObjective(-self.rho)
        self.model.optimize()
        for y,val in self.y.items():
            print(y,val.X)
        
    def directed_distance(self,f1,f2):
        if len(set(f1.children)&set(f2.children))>0:
            raise("Error: duplicated predicates:",set(f1.children)&set(f2.children))
        for p in f1.children:
            print(p)
            if p.type=="predicate":
                print(p)
                self.predicate_constraints(p,0)   
        for p in f2.children:
            print(p)
            if p.type=="predicate":
                print(p)
                self.predicate_constraints(p,1)
        self.model.setObjective(-self.rho)
        self.model.addConstr(self.z[f1,0]==1)
        self.model.addConstr(self.z[f2,0]==0)
        self.model.optimize()
        for y,val in self.y.items():
            print(y,val.X)
        if self.model.status==2:
            return self.rho.X
        else:
            return 0

def directed_distance(s,f1,f2):
    if len(set(f1.children)&set(f2.children))>0:
        raise("Error: duplicated predicates:",set(f1.children)&set(f2.children))
    for p in f1.children:
        print(p)
        if p.type=="predicate":
            print(p)
            s.predicate_constraints(p,0)   
    for p in f2.children:
        print(p)
        if p.type=="predicate":
            print(p)
            s.predicate_constraints(p,1)
    s.model.setObjective(-s.rho)
    s.model.addConstr(s.z[f1,0]==1)
    s.model.addConstr(s.z[f2,0]==0)
    s.model.optimize()
    if s.model.status==2:
        for y,val in s.y.items():
            print(y,val.X)
        return s.rho.X
    else:
        return 0