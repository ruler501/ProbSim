#!/usr/bin/env python3

import math
import operator
from collections import defaultdict
from fractions import Fraction
from functools import reduce 
from matplotlib import pyplot

zeroDice = False

if zeroDice:
    def sumProb(n, s):
        probList0 = [Fraction(1,s)]*s
        probListN = [Fraction(1,s)]*s
        for i in range(2, n + 1):
            tList = []
            for j in range(0, (s-1)*i + 1):
                tVal = 0
                for k in range(max(0,j-s+1), min(j,(s-1)*(i-1))+1):
                    tVal += probList0[j-k]*probListN[k]
                tList.append(tVal)
            probListN = tList
        return ((0,n*(s - 1)),probListN)
else:
    def sumProb(n, s):
        probList0 = [Fraction(1,s)]*s
        probListN = [Fraction(1,s)]*s
        for i in range(2, n + 1):
            tList = []
            for j in range(i, s*i + 1):
                tVal = 0
                for k in range(i-1, s*(i-1) + 1):
                    if 1 <= j-k and j - k <= s:
                        tVal += probList0[j-k-1]*probListN[k-i+1]
                tList.append(tVal)
            probListN = tList
        return ((n,n*s),probListN)
        
def prodProb(pN, pList):
    prodProb = []
    
    minN = pList[0][1][0][0]
    maxN = pList[0][1][0][1]
    
    for p in pList:
        minN = min(minN, p[1][0][0])
        maxN = max(maxN, p[1][0][1])
    
    for k in range(minN, maxN+1):
        tVal = Fraction()
        for p in pList:
            if p[1][0][0] <= k and k <= p[1][0][1]:
                tVal += pN[1][p[0] - pN[0][0]]*p[1][1][k-p[1][0][0]]
        prodProb.append(tVal)
    
    return ((minN,maxN),prodProb)

def aVal(p):
    av = 0
    for i in range(len(p[1])):
        av += p[1][i]*(p[0][0]+i)
    return av

def mVal(p):
#    print(p)
    tVal = 0
    for i in p[1]:
        tVal += i
#    print(tVal)
    med = p[1][0]
    i = 0
    while med < Fraction(1,2):
        i += 1
        med += p[1][i]
    return i + p[0][0]
    
def distVal(p,n):
    res = []
    i=0
    med = p[1][i]
    for j in range(n):
        while med < Fraction(j,n): 
            i += 1
            med += p[1][i]
        res.append(i + p[0][0])
    return res
    
def plotDists(p, *args):
    xPoints = []
    yPoints = []
    for i in range(p[0][1] - p[0][0]+1):
        xPoints.append(p[0][0]+i)
        yPoints.append(p[1][i])
    pargs = [list(xPoints), list(yPoints)]
    for d in args:
        xPoints = []
        yPoints = []
        for i in range(d[0][1] - d[0][0] + 1):
            xPoints.append(d[0][0]+i)
            yPoints.append(d[1][i])
        pargs.append(list(xPoints))
        pargs.append(list(yPoints))
    pyplot.plot(*pargs)
    pyplot.show()

plotDist = plotDists

def transformDists(tr,p,*args):
    v = defaultdict(Fraction)
    
    def recTransformDists(tr,targs,pargs,rargs):
        if rargs:
            for i in range(rargs[0][0][0], rargs[0][0][1] + 1):
                if rargs[0][1][i - rargs[0][0][0]] > 0:
                    recTransformDists(tr, [*targs, i], [*pargs, rargs[0][1][i - rargs[0][0][0]]], rargs[1:])
        else:
            v[tr(*targs)] += reduce(operator.mul, pargs)
    
    for i in range(p[0][0], p[0][1] + 1):
        recTransformDists(tr,[i,],[p[1][i - p[0][0]],],args)
    
    minVal = min(v.keys())
    maxVal = max(v.keys())
    r = ((minVal, maxVal), [])
    for i in range(minVal, maxVal + 1):
        r[1].append(v[i])
    return r

transformDist = transformDists

def distToBool(con, p):
    return transformDist(lambda x: int(con(x)), p)

def diffDists(p1,p2):
    minV = min(p1[0][0], p2[0][0])
    maxV = max(p1[0][1], p2[0][1])
    r = Fraction()
    for i in range(minV, maxV+1):
        if i < p1[0][0]:
            r += p2[1][i - p2[0][0]]**2
        elif i < p2[0][0]:
            r += p1[1][i - p1[0][0]]**2
        elif i > p1[0][1]:
            r += p2[1][i - p2[0][0]]**2
        elif i > p2[0][1]:
            r += p1[1][i - p1[0][0]]**2
        else:
            r += (p1[1][i - p1[0][0]] - p2[1][i - p2[0][0]])**2
    return r

def avDists(p1, p2):
    minV = min(p1[0][0], p2[0][0])
    maxV = max(p1[0][1], p2[0][1])
    res = ((minV, maxV), [])
    for i in range(minV, maxV+1):
        if i < p1[0][0]:
            r = p2[1][i - p2[0][0]]/2
        elif i < p2[0][0]:
            r = p1[1][i - p1[0][0]]/2
        elif i > p1[0][1]:
            r = p2[1][i - p2[0][0]]/2
        elif i > p2[0][1]:
            r = p1[1][i - p1[0][0]]/2
        else:
            r = (p1[1][i - p1[0][0]] + p2[1][i - p2[0][0]])/2
        res[1].append(r)
    return res

def distKLDist(p1, p2):
    minV = min(p1[0][0], p2[0][0])
    maxV = max(p1[0][1], p2[0][1])
    r = Fraction()
    for i in range(minV, maxV+1):
        if p2[0][0] <= i and i <= p2[0][1] and p1[0][0] <= i and i <= p1[0][1]:
            if p1[1][i - p1[0][0]] > 0 and p2[1][i - p2[0][0]] > 0:
                r += p1[1][i - p1[0][0]]*math.log(p1[1][i - p1[0][0]]/p2[1][i - p2[0][0]],2)
    return r

def distJSDist(p1,p2):
    av = avDists(p1,p2)
    return (distKLDist(p1,av) + distKLDist(p2,av))/2 

def cumDist(p):
    a = [sum(p[1][0:v+1]) for v in range(len(p[1]))]
    return (p[0],a)

def distToString(p):
    return "[" + "; ".join("({}, {:.2%})".format(i+p[0][0],float(v)) for i,v in enumerate(p[1])) + "]"

def printDist(p):
    print(distToString(p))
#def distToDists(tr, p, *args):
#     ind = [0]*(len(args)+1)
#     v = [None]*(len(args)+1)
#     v[0] = [None]*len(p[1]):
#     for i in range(len(args)):
#
#     while v[0] < len(p[1]):


def distToDists(tr, p, *args):
    v = defaultdict(Fraction)
    
    def recDistToDists(tr,targs,pargs,rargs):
        if rargs:
            for i in range(rargs[0][0][0], rargs[0][0][1] + 1):
                if rargs[0][1][i - rargs[0][0][0]] > 0:
                    recDistToDists(tr, [*targs, i], [*pargs, rargs[0][1][i - rargs[0][0][0]]], rargs[1:])
        else:
            dist = tr(*targs)
            d0 = reduce(operator.mul, pargs)
            for i,val in enumerate(dist[1]):
                v[dist[0][0] + i] += val*d0
    
    for i in range(p[0][0], p[0][1] + 1):
        recDistToDists(tr,[i,],[p[1][i - p[0][0]],],args)

    minVal = min(v.keys())
    maxVal = max(v.keys())
    r = ((minVal, maxVal), [])
    for i in range(minVal, maxVal + 1):
        r[1].append(v[i])
    return r

def obenDamage(x):
    if x == 20:
        return transformDist(lambda x,y: x+y+6, sumProb(2,8), sumProb(1,6))
    elif x + 10 > 14:
        return transformDist(lambda x,y: x+y+6, sumProb(1,8), sumProb(1,6))
    else:
        return sumProb(0,1)

def benDamage(x,y):
    d = sumProb(0,1)
    if x == 20:
        d = transformDist(lambda x,y: x+y+6, sumProb(2,8), sumProb(1,6))
    elif x + 10 > 14:
        d = transformDist(lambda x,y: x+y+6, sumProb(1,8), sumProb(1,6))
    if y == 20:
        d = transformDist(lambda x,y,z: x+y+z+6, sumProb(2,8), sumProb(1,6), d)
    elif y + 10 > 14:
        d = transformDist(lambda x,y,z: x+y+z+6, sumProb(1,8), sumProb(1,6), d)
    if x + 10 > 14 or y + 10 > 14:
        d = transformDist(lambda x,y: x+y, d, sumProb(1,8))
    return d

def newBenDamage(x):
    if x >= 19:
        return transformDist(lambda x: x+4, sumProb(2,8))
    elif x + 8 > 14:
        return transformDist(lambda x: x+4, sumProb(1,8))
    else:
        return sumProb(0,1)

def garrettDamage(x):
    if x == 20:
        return transformDist(lambda x: x+6, sumProb(4,6))
    elif x + 10 > 14:
        return transformDist(lambda x: x+6, sumProb(2,6))
    else:
        return sumProb(0,1)

def josephDamage(x,y):
    d = sumProb(0,1)
    if x == 20:
        d = transformDists(lambda x,y: x+y, sumProb(2,8), sumProb(4,6))
    elif x + 6 > 14:
        d = transformDists(lambda x,y: x+y, sumProb(1,8), sumProb(2,6))
    if y == 20 and x+6 < 15:
        d = transformDists(lambda x,y: x+y, sumProb(2,8), sumProb(4,6))
    elif y == 20:
        d = transformDists(lambda x,y: x+y, sumProb(2,8), d)
    elif y + 6 > 14 and x+6 < 15:
        d = transformDists(lambda x,y: x+y, sumProb(1,8), sumProb(2,6))
    elif y + 6 > 14:
        d = transformDists(lambda x,y: x+y, sumProb(1,8), d)
    return d

def ojosephDamage(x):
    d = sumProb(0,1)
    if x == 20:
        d = transformDists(lambda x,y: x+y, sumProb(3,8), sumProb(4,6))
    elif x + 6 > 14:
        d = transformDists(lambda x,y: x+y, sumProb(2,8), sumProb(2,6))
    return d

def dannyDamage(x,y,z):
    d = sumProb(0,1)
    h = 0
    c = False
    if x == 20:
        d = transformDist(lambda x: 2*x, sumProb(2,6))
        c = True
    elif x >= 8:
        h = 1
    if y == 20:
        c = True
        d = transformDists(lambda x,y: 2*x + y, sumProb(2,6), d)
    elif y >= 8:
        h += 1
    if z > 4 or c:
        for h in range(h):
            d = transformDists(lambda x,y: x + y, sumProb(2,6), d)
    return d

def devonDamage(x,z):
    d = sumProb(0,1)
    if x == 20:
        d = transformDist(lambda x,y: 2*(x + y), sumProb(3,4), sumProb(1,2))
    elif x > 1:
        if z > 4:
            d = transformDist(lambda x,y: x + y, sumProb(3,4), sumProb(1,2))
    return d

def weirdDist(*args):
    min1 = (7,-1)
    min2 = (7,-1)
    for i,val in enumerate(args):
        if val < min1[0]:
            min1 = (val,i)
        elif val < min2[0]:
            min2 = (val,i)
    s = 0
    for i,val in enumerate(args):
        if not (i == min1[1] or i == min2[1]):
            s += val
    return int(math.ceil(s/2))

def dropOne(*args):
    return reduce(operator.add, args) - min(args)



#d = distToDists(lambda x: transformDists(lambda y: int(math.ceil((y - x - 1)/4)), sumProb(x,2)), sumProb(1,167)) 
#print(d)
#print("Average: ", float(aVal(d)))
#print("Median: ", mVal(d))
#print("20 Piece distribution: ", distVal(d,20))
#print(distJSDist(sumProb(1,20),d))
#plotDist(d, sumProb(1,20))

#d6s = [sumProb(1,6)]*8
#d = transformDists(weirdDist, *d6s)
#d = transformDists(dropOne, *d6s)
#d = sumProb(50,4)
d6s = [sumProb(1,6)]*2
d = distToDists(lambda x,y: sumProb(x-1, y), *d6s)
#cd = cumDist(d)
printDist(d)
#printDist(cd)
#print("Average: ", float(aVal(d)))
#print("Median: ", mVal(d))
#print("10 Piece distribution: ", distVal(d,10))
#plotDist(d,cd)

#

#d = distToDists(devonDamage, sumProb(1,20), sumProb(1,10))
#print(d)
#print("Average: ", float(aVal(d)))
#print("Median: ", mVal(d))
#print("10 Piece distribution: ", distVal(d,10))
#plotDist(d)

#hd = sumProb(1,20)
#hd = transformDists(lambda x,y: min(x,y), sumProb(1,20), sumProb(1,20))
#d1 = distToDists(obenDamage, hd)
#
#hd = transformDists(lambda x,y: max(x,y), sumProb(1,20), sumProb(1,20))
#d2 = distToDists(ojosephDamage, hd)
#
#hd = transformDists(lambda x,y: max(x,y), sumProb(1,20), sumProb(1,20))
#d3 = distToDists(garrettDamage, hd)
##hd = sumProb(1,20)
##d = distToDists(garrettDamage, hd)
##d3 = transformDists(lambda x,y: x + y, d3, d)
#
#print(d1)
#print(d2)
#print(d3)
#print("Average: ", float(aVal(d1)))
#print("Average: ", float(aVal(d2)))
#print("Average: ", float(aVal(d3)))
#print("Median: ", mVal(d1))
#print("Median: ", mVal(d2))
#print("Median: ", mVal(d3))
#print("10 Piece distribution: ", distVal(d1,10))
#print("10 Piece distribution: ", distVal(d2,10))
#print("10 Piece distribution: ", distVal(d3,10))
#plotDist(d1, d2, d3)
""" 
tStr = input("Dice roll for determining pool ")
n,s1 = tStr.split('d')
n = int(n)
s1 = int(s1)
tStr = input("Size of dice in pool ")
s2 = int(tStr)

diceN = sumProb(n,s1)
listDiceM = []

if zeroDice:
    for i in range(0,n*(s1-1) + 1):
        listDiceM.append( (i,sumProb(i,s2)) )
else:
    for i in range(n,n*s1 + 1):
        listDiceM.append( (i,sumProb(i,s2)) )
        
#print(diceN)
#print(listDiceM)
#print()
jointP = prodProb(diceN,listDiceM)
#print(jointP)
print("Average: ", int(round(float(aVal(jointP)))))
print("Median: ", mVal(jointP))
print("10 Piece distribution: ", distVal(jointP,10))
plotDist(jointP)
"""
