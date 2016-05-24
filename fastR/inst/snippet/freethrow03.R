pp <- dbinom(5,5,0.80); pp  # a) prob make 5 straight (= success)
dgeom(1,pp)   # b) succeed with 1 miss (correct answer)
0.20 * pp     # miss first then make 5 straight (INCORRECT)
1 - pgeom(1,pp)             # c) miss more than one shot before success
probs <- dgeom(0:15,pp)     # d) 
#
myplot <- xyplot(probs~0:15,main="Freddie",
                 xlab="misses",ylab="probability")
#################################################################
pp <- dbinom(5,5,0.70)    # a) prob make 5 straight (= success)
pp                             

dgeom(1,pp)     # b) succeed with 1 miss (correct answer)
0.20 * pp       # miss first then make 5 straight (_not_ the answer)

1 - pgeom(1,pp)           # c) miss more than one shot before success
probs <- dgeom(0:15,pp)   # d)
#
myplot <- xyplot(probs~0:15,main="Frank",
                 xlab="misses",ylab="probability")
