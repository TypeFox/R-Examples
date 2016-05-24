# ----------------------------------------------------------------------
# EXAMPLES FOR TESTING THE FUNCTIONALITY OF RANDOM SEQUENCE GENERATION
# ----------------------------------------------------------------------
# require(randomizeR)

# ------
# K = 2
# ------

# generate an object with the right class and parameters
# (this is not random)
pbr <- pbrPar(c(4, 4))
rar <- rarPar(10)
mp <- mpPar(10, 2)
bsd <- bsdPar(10, 2)
ebc <- ebcPar(10, 0.667) 
cr <- crPar(10)
had <- hadaPar(10)
ud <- udPar(10, 0,1)
tbd <- tbdPar(10)

# generate one sequence
# this is random
genSeq(pbr)
pbr

genSeq(rar)
genSeq(mp)
genSeq(bsd)
genSeq(ebc)
genSeq(cr)
genSeq(had)
genSeq(ud)
genSeq(tbd)

genSeq(pbr, 4)
genSeq(rar, 4)
genSeq(mp, 4)
genSeq(bsd, 4)
genSeq(ebc, 4)
genSeq(cr, 4)
genSeq(had, 4)
genSeq(ud, 4)
genSeq(tbd, 4)

# generate several sequences
# this is random

# generate all sequences
# this is not random
getAllSeq(pbr)
getAllSeq(rar)
getAllSeq(mp)
getAllSeq(bsd)
getAllSeq(ebc)
getAllSeq(cr)
getAllSeq(had)




# ---------------------
# K > 2 OR ratio!= 1:1
# ---------------------

# unequal allocation or K > 2 supported:

randPar(4, ratio = c(2, 1)) # expect error

cr <- crPar(10)
genSeq(cr)
genSeq(cr, 4)


randPar(10)
randPar(10, ratio = c(4, 3))
crPar(10, ratio = c(4, 3))


pbr<-pbrPar(c(4, 4), K = 3) # expect error
pbr<-pbrPar(c(4, 4),4) # each treament once per block
pbr<-pbrPar(c(4, 4),ratio = c(2, 1, 2),K = c(2, 3)) # expect error
pbr<-pbrPar(bc = c(3, 6),ratio = c(2, 1))
pbr<-pbrPar(bc = 3,ratio = c(2, 1))
genSeq(pbr)
genSeq(pbr, 4)

rar<-rarPar(10, K = 3) # expect error
rar<-rarPar(12, K = 3)
rar<-rarPar(12, ratio = c(1, 3)) # this should not work
rar<-rarPar(12,ratio = c(3, 4)) # this should not work
genSeq(rar)
genSeq(rar, 4)

mp<-mpPar(10, 2)     # works as usual
mp<-mpPar(10, 2, K = 3) # expect error: unused argument
mp<-new("mpPar", N = 10, mti = 2, K = 3,ratio = rep(1, 3)) # expect error
mp<-mpPar(10, 2,ratio = c(2, 3)) # works
genSeq(mp)
genSeq(mp,12)@M


# don't support K > 2 or ratio!=c(1,1)
# expect error: unused argument
# this means we have not implement K or ratio as possible input in the 
# constructor function.
# (we might want to change this behaviour)


bsd<-bsdPar(10, 2, K = 3) # expect error: unused argument
bsd<-new("bsdPar", N = 10, mti = 2, K = 3, ratio = rep(1, 3))# expect error
bsd<-new("bsdPar", N = 10, mti = 2, K = 2, ratio = c(1, 2))# expect error

ebc<-ebcPar(10, 0.667, K = 3) # expect error: unused argument
bsd<-new("ebcPar", N=10, K=2, ratio=c(1,2))# expect error

had<-hadaPar(10, K = 3) # expect error: unused argument
had<-hadaPar(10, ratio = c(1, 2)) # expect error: unused argument

had<-hadaPar(10, K = 3) # expect error: unused argument
ud<-udPar(10, 0,1, ratio=c(1, 2)) # expect error: unused argument

tbd<-tbdPar(10, K = 3) # expect error: unused argument
tbd<-tbdPar(10, ratio = c(1, 2)) # expect error: unused argument
