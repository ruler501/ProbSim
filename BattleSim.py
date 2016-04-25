from ProbabilityDist import *

#d20 = sumProb(1,20)
#adv = transformDists(lambda x,y: max(x,y), d20, d20)
#dadv = transformDists(lambda x,y: min(x,y), d20, d20)
#
##jd = distToDists(josephDamage, adv, adv)
#gd = transformDists(lambda x,y: x+y, distToDists(garrettDamage, adv), distToDists(garrettDamage, d20))
#bd = distToDists(benDamage, dadv, dadv)
#def bossDamage(x):
#    if x == 20:
#        return transformDists(lambda x,y: x+y, sumProb(4,8), sumProb(6,6))
#    elif x >= 5:
#        return transformDists(lambda x,y: x+y, sumProb(2,8), sumProb(3,6))
#    return sumProb(0,1)

hd = transformDists(lambda x: x - 4, sumProb(6,4))
hl = 30

bd = sumProb(1,20)
bl = 30

#boss = distToDists(bossDamage, d20)

def lifeTransform(l, xd, yd):
    xl = l // (hl + 1)
    yl = l % (hl + 1)
    if xl > 0:
        yl -= xd
    if yl > 0:
        xl -= yd
    xl = max(xl, 0)
    yl = max(yl, 0)
    return xl*(hl + 1) + yl

life = transformDist(lambda x: x + bl*(hl + 1) + hl, sumProb(0,1))

lifes = [transformDists(lifeTransform, life, hd, bd)]

d = transformDist(lambda x: int(x//(hl+1) <= 0 and x%(hl+1) > 0), lifes[0])
endd = transformDist(lambda x: int(x//(hl+1) <= 0 or x%(hl+1) <= 0), lifes[0])
printDist(d)
printDist(endd)
print()

i = 0
while endd[1][0] > 0.001:
    lifes.append(transformDists(lifeTransform, lifes[i], hd, bd))
    i += 1
    d = transformDist(lambda x: int(x//(hl+1) <= 0 and x%(hl+1) > 0), lifes[i])
    endd = transformDist(lambda x: int(x//(hl+1) <= 0 or x%(hl+1) <= 0), lifes[i])
    printDist(d)
    printDist(endd)
    print()
