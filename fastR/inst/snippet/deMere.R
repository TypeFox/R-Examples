1 - dbinom(0,4,1/6)   # P(at least one 6 in 4 tries) 
pgeom(3,1/6)          # P(fail at most 3 times before getting a 6)
1 - dbinom(0,24,1/36) # P(at least one double 6 in 24 tries)
pgeom(23,1/36)        # P(fail at most 23 times before getting double 6)
