#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2013 the statnet development team
######################################################################

########################################################################
# This file contains the testing suite for the "query" methods, i.e.:
#     is.active
#     is.adjacent.active
#########################################################################

require(networkDynamic)
require(testthat)


#-------------------- IS.ACTIVE TESTS ---------------------------------
# Notes:
#  -these are meant to be exhaustive except where noted otherwise
#  -the list named "spells" created below is helpful for checking
#   that the tests are correctly written
#  -for understanding the comments, let
#      L=finite number lower than spell boundary,
#      H=finite number higher than spell boundary,
#      M=finite number within a spell boundary
#----------------------------------------------------------------------

# create a test network
anet <- network.initialize(100)
set.seed(10)
heads <- sample(100, 150, replace=TRUE)
set.seed(25)
tails <- sample(100, 150, replace=TRUE)
add.edges(anet, tails, heads)
anet.copy <- anet

cat("testing is.active ... ")

# proper behavior for bad inputs?  (not exhaustive)
#  -- bad network input
a1 = try(is.active(3, 3, 3), T)                  # properly handeled
#  -- bad onset
a2 = try(is.active(anet, "a", 3), T)             # properly handeled
a3 = try(is.active(anet, NULL, 3), T)            # properly handeled
# -- bad terminus
a4 = try(is.active(anet, 3, "b"), T)             # properly handeled
# -- bad edges
a5 = try(is.active(anet, 3, 10, e=224), T)       # properly handeled     
a6 = try(is.active(anet, 3, 10, e=-2), T)        # properly handeled    
a7 = try(is.active(anet, 3, 10, e=NULL), T)      # properly handeled    
# -- bad vertices
a8 = try(is.active(anet, 3, 10, v=224), T)       # properly handeled     
a9 = try(is.active(anet, 3, 10, v=-2), T)        # properly handeled    
# -- bad rule
a10 = try(is.active(anet, 3, 20, e=1, rule=10), T)   # properly handeled by R
a11 = try(is.active(anet, 15, 25, e=5, rule=c("all", "any")), T)  # properly handeled by R
# -- bad active.default
a12 = try(is.active(anet, 3, 10, e=1, active.default=4), T)   # returns 4
# -- bad onset & terminus combo
a13 = try(is.active(anet, 10, 3, e=5), T)        # properly handeled
# -- bad edges & vertices combo
a14 = try(is.active(anet, 3, 10, e=5, v=4), T)   # properly handeled


# now, for completely valid input ...
activate.edges(anet, c(-Inf, -Inf,  10, 10, 10),
                     c( Inf,   20, Inf, 20, 10), e=2:6)
activate.edges(anet, 10, 20, e=7)
activate.edges(anet, 20, 30, e=7)
activate.edges(anet, 30, 30, e=7)
activate.edges(anet, 40, 50, e=7)
deactivate.edges(anet, e=8)
spells = lapply(lapply(anet$mel[1:8],"[[","atl"),"[[","active")
  
# point queries that should be TRUE
b1  = is.active(anet, at=4, e= 1)      # for any &      NULL  
b2  = is.active(anet, at=3, e= 2)      # for M &     (-Inf, Inf)
b3  = is.active(anet, at=-Inf, e= 2)   # for -Inf &  (-Inf, Inf)
b4  = is.active(anet, at=Inf, e= 2)    # for Inf &   (-Inf, Inf)
b5  = is.active(anet, at=5, e=3)       # for M &     (-Inf, b)
b6  = is.active(anet, at=15, e=7)      # for M &     (-Inf, b)
b7  = is.active(anet, at=-Inf, e=3)    # for -Inf &  (-Inf, b)
b8  = is.active(anet, at=50, e=4)      # for M &     (a, Inf)
b9 = is.active(anet, at=45, e=7)       # for M &     (a, Inf)
b10 = is.active(anet, at=10, e=4)      # for a &     (a, Inf)
b11 = is.active(anet, at=20, e=7)      # for a &     (a, Inf)
b12 = is.active(anet, at=Inf, e=4)     # for Inf &   (a, Inf)
b13 = is.active(anet, at=15, e=4)      # for M &     (a, b)
b14 = is.active(anet, at=25, e=7)      # for M &     (a, b)
b15 = is.active(anet, at=10, e=4)      # for a &     (a, b)
b16 = is.active(anet, at=20, e=7)      # for a &     (a, b)
b17 = is.active(anet, at=10, e=6)      # for a &     (a, a)
b18 = is.active(anet, at=30, e=7)      # for a &     (a, a)

b.tests = paste("b", seq(1,18), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  stop(paste("is.active is returning improper FALSE values in tests",
             bad.tests))
}


# point queries that should be FALSE
c1  = is.active(anet, at=4, e=1,active.default=F)  # for any & NULL
c2  = is.active(anet, at=3, e=8)       # for 'never active' NULL spell
c3  = is.active(anet, at=25, e=3)      # for H &    (-Inf, b)
c4  = is.active(anet, at=55, e=7)      # for H &    (-Inf, b)
c5  = is.active(anet, at=20, e=3)      # for b &    (-Inf, b)
c6  = is.active(anet, at=50, e=7)      # for b &    (-Inf, b)
c7  = is.active(anet, at=Inf, e=3)     # for Inf &  (-Inf, b)
c8  = is.active(anet, at=5, e=4)       # for L &    (a, Inf)
c9  = is.active(anet, at=35, e=7)      # for L &    (a, Inf)
c10 = is.active(anet, at=-Inf, e=4)    # for -Inf & (a, Inf)
c11 = is.active(anet, at=5, e=5)       # for L &    (a, b)
c12 = is.active(anet, at=5, e=7)       # for L &    (a, b)
c13 = is.active(anet, at=25, e=5)      # for H &    (a, b)
c14 = is.active(anet, at=55, e=7)      # for H &    (a, b)
c15 = is.active(anet, at=20, e=5)      # for b &    (a, b)
c16 = is.active(anet, at=50, e=7)      # for b &    (a, b)
c17 = is.active(anet, at=0, e=6)       # for L &    (a, a)

c.tests = paste("c", seq(1,17), sep="")
c.results= sapply(c.tests, function(x){eval(parse(text=x))})
if(any(c.results)){
  bad.tests = paste("c", which(c.results), sep="", collapse=" ")
  stop(paste("is.active is returning improper TRUE values in tests",
             bad.tests))
}

# interval queries that should be TRUE using 'rule=all'
d1  = is.active(anet, 4, 5, e= 1, rule='all')       # for any &           NULL
d2  = is.active(anet, -Inf, Inf, e= 1, rule='all')  # for any &           NULL  
d3  = is.active(anet, -Inf, Inf, e= 2, rule='all')  # for (-Inf,Inf) & (-Inf, Inf)
d4  = is.active(anet, -Inf, -Inf, e= 2, rule='all') # for (-Inf,-Inf)& (-Inf, Inf)
d5  = is.active(anet, Inf, Inf, e= 2, rule='all')   # for (Inf,Inf) &  (-Inf, Inf)
d6  = is.active(anet, 0, Inf, e= 2, rule='all')     # for (M,Inf) &    (-Inf, Inf)
d7  = is.active(anet, -Inf, 0, e= 2, rule='all')    # for (-Inf,M) &   (-Inf, Inf)
d8  = is.active(anet, 0, 0, e= 2, rule='all')       # for (M,M) &      (-Inf, Inf)
d9  = is.active(anet, 0, 10, e= 2, rule='all')      # for (M,M) &      (-Inf, Inf)
d10 = is.active(anet, -Inf, -Inf, e=3, rule='all')  # for (-Inf,-Inf)& (-Inf, b)
d11 = is.active(anet, -Inf, 20, e=3, rule='all')    # for (-Inf,b) &   (-Inf, b)
d12  = is.active(anet, -Inf, 0, e=3, rule='all')    # for (-Inf,M) &   (-Inf, b)
d13 = is.active(anet, -5, 0, e=3, rule='all')       # for (M,M) &      (-Inf, b)
d14 = is.active(anet, 0, 20, e=3, rule='all')       # for (M,b) &      (-Inf, b)
d15 = is.active(anet, Inf, Inf, e=4, rule='all')    # for (Inf,Inf) &  (a, Inf)
d16 = is.active(anet, 10, Inf, e=4, rule='all')     # for (a,Inf) &    (a, Inf)
d17 = is.active(anet, 30, Inf, e=4, rule='all')     # for (M,Inf) &    (a, Inf)
d18 = is.active(anet, 10, 30, e=4, rule='all')      # for (a,M) &      (a, Inf)
d19 = is.active(anet, 30, 40, e=4, rule='all')      # for (M,M) &      (a, Inf)
d20 = is.active(anet, 10, 20, e=5, rule='all')      # for (a,b) &      (a, b)
d21 = is.active(anet, 10, 15, e=5, rule='all')      # for (a,M) &      (a, b)
d22 = is.active(anet, 15, 20, e=5, rule='all')      # for (M,b) &      (a, b)
d23 = is.active(anet, 10, 10, e=6, rule='all')      # for (a,a) &      (a, a)

d.tests = paste("d", seq(1,23), sep="")
d.results= sapply(d.tests, function(x){eval(parse(text=x))})
if(any(!d.results)){
  bad.tests = paste("d", which(!d.results), sep="", collapse=" ")
  stop(paste("is.active is returning improper FALSE values in tests",
             bad.tests))
}



# interval queries that should be TRUE using default 'rule=any'
e1  = is.active(anet, 4, 5, e= 1)       # for any &           NULL
e2  = is.active(anet, -Inf, Inf, e= 1)  # for any &           NULL  
e3  = is.active(anet, -Inf, Inf, e= 2)  # for (-Inf,Inf) & (-Inf, Inf)
e4  = is.active(anet, -Inf, -Inf, e=2)  # for (-Inf,-Inf)& (-Inf, Inf)
e5  = is.active(anet, Inf, Inf, e= 2)   # for (Inf,Inf) &  (-Inf, Inf)
e6  = is.active(anet, 0, Inf, e= 2)     # for (M,Inf) &    (-Inf, Inf)
e7  = is.active(anet, -Inf, 0, e= 2)    # for (-Inf,M) &   (-Inf, Inf)
e8  = is.active(anet, 0, 10, e= 2)      # for (M,M) &      (-Inf, Inf)
e9  = is.active(anet, -Inf, -Inf, e=3)  # for (-Inf,-Inf)& (-Inf, b)
e10 = is.active(anet, -Inf, 20, e=3)    # for (-Inf,b) &   (-Inf, b)
e11 = is.active(anet, -Inf, 10, e=3)    # for (-Inf,M) &   (-Inf, b)
e12 = is.active(anet, -Inf, 30, e=3)    # for (-Inf,H) &   (-Inf, b)
e13 = is.active(anet, 0, 20, e=3)       # for (M,b) &      (-Inf, b)
e14 = is.active(anet, 0, 10, e=3)       # for (M,M) &      (-Inf, b)
e15 = is.active(anet, 0, 30, e=3)       # for (M,H) &      (-Inf, b)
e16 = is.active(anet, Inf, Inf, e=4)    # for (Inf,Inf) &  (a, Inf)
e17 = is.active(anet, 10, Inf, e=4)     # for (a,Inf) &    (a, Inf)
e18 = is.active(anet, 20, Inf, e=4)     # for (M,Inf) &    (a, Inf)
e19 = is.active(anet, 0, Inf, e=4)      # for (L,Inf) &    (a, Inf)
e20 = is.active(anet, 10, 20, e=4)      # for (a,M) &      (a, Inf)
e21 = is.active(anet, 0, 20, e=4)       # for (L,M) &      (a, Inf)
e22 = is.active(anet, 20, 30, e=4)      # for (M,M) &      (a, Inf)
e23 = is.active(anet, -Inf, Inf, e=5)   # for (-Inf, Inf) & (a, b)
e24 = is.active(anet,  0, Inf, e=5)     # for (L, Inf) &    (a, b)
e25 = is.active(anet, 10, Inf, e=5)     # for (a, Inf) &    (a, b)
e26 = is.active(anet, 15, Inf, e=5)     # for (M, Inf) &    (a, b)
e27 = is.active(anet, -Inf, 30, e=5)    # for (-Inf, H) &   (a, b)
e28 = is.active(anet, -Inf, 20, e=5)    # for (-Inf, b) &   (a, b)
e29 = is.active(anet, -Inf, 15, e=5)    # for (-Inf, M) &   (a, b)
e30 = is.active(anet, 10, 20, e=5)      # for (a,b) &       (a, b)
e31 = is.active(anet, 10, 30, e=5)      # for (a,H) &       (a, b)
e32 = is.active(anet, 10, 15, e=5)      # for (a,M) &       (a, b)
e33 = is.active(anet, 0, 20, e=5)       # for (L,b) &       (a, b)
e34 = is.active(anet, 15, 20, e=5)      # for (M,b) &       (a, b)
e35 = is.active(anet, 0, 30, e=5)       # for (L,H) &       (a, b)
e36 = is.active(anet, 0, 15, e=5)       # for (L,M) &       (a, b)
e37 = is.active(anet, 15, 30, e=5)      # for (M,H) &       (a, b)
e38 = is.active(anet, 15, 18, e=5)      # for (M,M) &       (a, b)
e39 = is.active(anet, -Inf, Inf, e=6)   # for (-Inf,Inf) &  (a, a)
e40 = is.active(anet, 0, Inf, e=6)      # for (L,Inf) &     (a, a)
e41 = is.active(anet, -Inf, 20, e=6)    # for (-Inf,H) &    (a, a)
e42 = is.active(anet, 10, 20, e=6)      # for (a, H) &      (a, a)
e43 = is.active(anet, 30, 35, e=7)      # for (a, H) &      (a, a)
e44 = is.active(anet, 10, 10, e=6)      # for (a, a) &      (a, a)
e45 = is.active(anet, 0, 20, e=6)       # for (L,H) &       (a, a)

e.tests = paste("e", seq(1,45), sep="")
e.results= sapply(e.tests, function(x){eval(parse(text=x))})
if(any(!e.results)){
  bad.tests = paste("e", which(!e.results), sep="", collapse=" ")
  stop(paste("is.active is returning improper FALSE values in tests",
             bad.tests))
}


# interval queries that should be FALSE using 'rule=all'
f1  = is.active(anet, 4, 5, e= 1, rule='all', active.default=F)       # for any & NULL
f2  = is.active(anet, -Inf, Inf, e= 1, rule='all', active.default=F)  # for any & NULL  
f3  = is.active(anet, 3, 5, e= 8, rule= 'all')      # for 'never active' NULL spell
f4  = is.active(anet, -Inf, Inf, e= 8, rule= 'all') # for 'never active' NULL spell
f5  = is.active(anet, -Inf, Inf, e=3, rule='all')   # for (-Inf,Inf) & (-Inf, b)
f6  = is.active(anet, -Inf, 30, e=3, rule='all')    # for (-Inf,H) &   (-Inf, b)
f7 = is.active(anet, 0, 30, e=3, rule='all')        # for (M,H) &      (-Inf, b)
f8 = is.active(anet, 30, 40, e=3, rule='all')       # for (H,H) &      (-Inf, b)
f9 = is.active(anet, -Inf, Inf, e=4, rule='all')    # for (-Inf,Inf) & (a, Inf)
f10 = is.active(anet, 0, Inf, e=4, rule='all')      # for (L,Inf) &    (a, Inf)
f11 = is.active(anet, 0, 20, e=4, rule='all')       # for (L,M) &      (a, Inf)
f12 = is.active(anet, 0, 5, e=4, rule='all')        # for (L,L) &      (a, Inf)
f13 = is.active(anet, -Inf, Inf, e=5, rule='all')   # for (-Inf,Inf) & (a, b)
f14 = is.active(anet, -Inf, 0, e=5, rule='all')     # for (-Inf,L) &   (a, b)
f15 = is.active(anet, -Inf, 10, e=5, rule='all')    # for (-Inf,a) &   (a, b)
f16 = is.active(anet, -Inf, 15, e=5, rule='all')    # for (-Inf,M) &   (a, b)
f17 = is.active(anet, -Inf, 20, e=5, rule='all')    # for (-Inf,b) &   (a, b)
f18 = is.active(anet, -Inf, 30, e=5, rule='all')    # for (-Inf,H) &   (a, b)
f19 = is.active(anet, 0, Inf, e=5, rule='all')      # for (L,Inf) &    (a, b)
f20 = is.active(anet, 10, Inf, e=5, rule='all')     # for (a,Inf) &    (a, b)
f21 = is.active(anet, 15, Inf, e=5, rule='all')     # for (M,Inf) &    (a, b)
f22 = is.active(anet, 20, Inf, e=5, rule='all')     # for (b,Inf) &    (a, b)
f23 = is.active(anet, 30, Inf, e=5, rule='all')     # for (H,Inf) &    (a, b)
f24 = is.active(anet, 10, 30, e=5, rule='all')      # for (a,H) &      (a, b)
f25 = is.active(anet,  0, 10, e=5, rule='all')      # for (L,a) &      (a, b)
f26 = is.active(anet,  0, 20, e=5, rule='all')      # for (L,b) &      (a, b)
f27 = is.active(anet, 20, 30, e=5, rule='all')      # for (b,H) &      (a, b)
f28 = is.active(anet,  0, 15, e=5, rule='all')      # for (L,M) &      (a, b)
f29 = is.active(anet, 15, 30, e=5, rule='all')      # for (M,H) &      (a, b)
f30 = is.active(anet,  0, 30, e=5, rule='all')      # for (L,H) &      (a, b)
f31 = is.active(anet, -5,  5, e=5, rule='all')      # for (L,L) &      (a, b)
f32 = is.active(anet, 30, 40, e=5, rule='all')      # for (H,H) &      (a, b)
f33 = is.active(anet, -Inf,Inf, e=6, rule='all')    # for (-Inf,Inf) & (a, a)
f34 = is.active(anet, -Inf, 0, e=6, rule='all')     # for (-Inf,L) &   (a, a)
f35 = is.active(anet, -Inf, 10, e=6, rule='all')    # for (-Inf,a) &   (a, a)
f36 = is.active(anet, -Inf, 20, e=6, rule='all')    # for (-Inf,H) &   (a, a)
f37 = is.active(anet, 0, Inf, e=6, rule='all')      # for (L,Inf) &    (a, a)
f38 = is.active(anet, 30, Inf, e=6, rule='all')     # for (H,Inf) &    (a, a)
f39 = is.active(anet, 0, 10, e=6, rule='all')       # for (L,a) &      (a, a)
f40 = is.active(anet, 10, 20, e=6, rule='all')      # for (a,H) &      (a, a)
f41 = is.active(anet, 30, 35, e=7, rule='all')      # for (a,H) &      (a, a)
f42 = is.active(anet, -5, 0, e=6, rule='all')       # for (L,L) &      (a, a)
f43 = is.active(anet, 0, 40, e=6, rule='all')       # for (L,H) &      (a, a)
f44 = is.active(anet, 30, 40, e=6, rule='all')      # for (H,H) &      (a, a)

f.tests = paste("f", seq(1,43), sep="")
f.results= sapply(f.tests, function(x){eval(parse(text=x))})
if(any(f.results)){
  bad.tests = paste("f", which(f.results), sep="", collapse=" ")
  stop(paste("is.active is returning improper TRUE values in tests",
             bad.tests))
}


# interval queries that should be FALSE using default 'rule=any'
g1  = is.active(anet, 4, 5, e= 1, active.default=F)       # for any & NULL
g2  = is.active(anet, -Inf, Inf, e= 1, active.default=F)  # for any & NULL  
g3  = is.active(anet, 3, 5, e= 8)       # for 'never active' NULL spell
g4  = is.active(anet, -Inf, Inf, e= 8)  # for 'never active' NULL spell
g5  = is.active(anet, Inf, Inf, e=3)    # for (Inf, Inf) &  (-Inf, b)
g6  = is.active(anet, 30, Inf, e=3)     # for (H, Inf) &    (-Inf, b)
g7  = is.active(anet, 20, 30, e=3)      # for (b, H) &      (-Inf, b)
g8  = is.active(anet, 30, 40, e=3)      # for (H, H) &      (-Inf, b)
g9  = is.active(anet, -Inf, -Inf, e=4)  # for (-Inf,-Inf) & (a, Inf)
g10 = is.active(anet, -Inf, 0, e=4)     # for (-Inf,L) &    (a, Inf)
g11 = is.active(anet, 0, 10, e=4)       # for (L, a) &      (a, Inf)
g12 = is.active(anet, 0, 5, e=4)        # for (L, L) &      (a, Inf)
g13 = is.active(anet, -Inf, -Inf, e=5)  # for (-Inf,-Inf) & (a, b)
g14 = is.active(anet, Inf, Inf, e=5)    # for (Inf, Inf) &  (a, b)
g15 = is.active(anet, -Inf, 0, e=5)     # for (-Inf, L) &   (a, b)
g16 = is.active(anet, 30, Inf, e=5)     # for (H, Inf) &    (a, b)
g17 = is.active(anet, -5, 10, e=5)      # for (L,a) &       (a, b)
g18 = is.active(anet, 20, 30, e=5)      # for (b,H) &       (a, b)
g19 = is.active(anet, -5, 5, e=5)       # for (L,L) &       (a, b)
g20 = is.active(anet, 30, 40, e=5)      # for (H,H) &       (a, b)
g21 = is.active(anet, -Inf,-Inf, e=6)   # for (-Inf,-Inf) & (a, a)
g22 = is.active(anet, -Inf, 0, e=6)     # for (-Inf,L) &    (a, a)
g23 = is.active(anet, -Inf, 10, e=6)    # for (-Inf,a) &    (a, a)
g24 = is.active(anet, 30, Inf, e=6)     # for (H,Inf) &     (a, a)
g25 = is.active(anet, 0, 10, e=6)       # for (L,a) &       (a, a)
g26 = is.active(anet, -5, 0, e=6)       # for (L1,L2) &     (a, a)
g27 = is.active(anet, 30, 40, e=6)      # for (H1,H2) &     (a, a)

g.tests = paste("g", seq(1,27), sep="")
g.results= sapply(g.tests, function(x){eval(parse(text=x))})
if(any(g.results)){
  bad.tests = paste("g", which(g.results), sep="", collapse=" ")
  stop(paste("is.active is returning improper TRUE values in tests",
             bad.tests))
}

expect_equal(is.active(network.initialize(0),at=1),logical(0))

# check for warning when eid specified for deleted edge
test<-network.initialize(3)
add.edges(test,1:3,c(2,3,1))
delete.edges(test,eid=2)
expect_warning(is.active(test,e=1:3,onset=-Inf,terminus=Inf),'correspond to deleted edges')
expect_equal(suppressWarnings(is.active(test,e=1:3,onset=-Inf,terminus=Inf)),c(TRUE,TRUE))

# check that 'earliest' and 'latest' rules accepted as equivilent to 'any' #544
test<-network.initialize(1)
activate.vertices(test,at=1)
expect_true(is.active(test,onset=0,terminus=2,rule='earliest',v=1))
expect_true(is.active(test,onset=0,terminus=2,rule='latest',v=1))

cat("ok\n")


#----------------- IS.ADJACENT.ACTIVE TESTS-----------------------------
# Notes:
#  --this function is very brief, only really calling functions that
#    are well tested (get.edgeIDs by time, is.active by me), thus
#    the testing here is quite minimal
#-------------------------------------------------------------------

cat("testing is.adjacent.active ... ")
iov = as.matrix(anet, matrix.type="edgelist")[1:8,]

# tests that should return true for point queries
b1 = is.adjacent.active(anet, iov[2,1], iov[2,2], at=-Inf)  # at -Inf
b2 = is.adjacent.active(anet, iov[3,1], iov[3,2], at=-Inf)  # at -Inf
b3 = is.adjacent.active(anet, iov[2,1], iov[2,2], at=Inf)   # at Inf
b4 = is.adjacent.active(anet, iov[4,1], iov[4,2], at=Inf)   # at Inf
b5 = is.adjacent.active(anet, iov[4,1], iov[4,2], at=10)    # at a
b6 = is.adjacent.active(anet, iov[7,1], iov[7,2], at=10)    # at a
b7 = is.adjacent.active(anet, iov[7,1], iov[7,2], at=30)    # at a
b8 = is.adjacent.active(anet, iov[7,1], iov[7,2], at=30)    # at b
b9 = is.adjacent.active(anet, iov[7,1], iov[7,2], at=25)    # at M
b10 = is.adjacent.active(anet, iov[3,1], iov[3,2], at=10)   # at M
b11 = is.adjacent.active(anet, iov[4,1], iov[4,2], at=25)   # at M
b12 = is.adjacent.active(anet, iov[1,1], iov[1,2], at=25)   # at M

# tests that should return false for point queries b/c of non-adjacency
b13 = !is.adjacent.active(anet, iov[2,2], iov[2,1], at=-Inf)  # at -Inf
b14 = !is.adjacent.active(anet, iov[3,2], iov[3,1], at=-Inf)  # at -Inf
b15 = !is.adjacent.active(anet, iov[2,2], iov[2,1], at=Inf)   # at Inf
c16 = !is.adjacent.active(anet, iov[4,2], iov[4,1], at=-Inf)  # at Inf
b17 = !is.adjacent.active(anet, iov[4,2], iov[4,1], at=10)    # at a
b18 = !is.adjacent.active(anet, iov[7,2], iov[7,1], at=10)    # at a
b19 = !is.adjacent.active(anet, iov[7,2], iov[7,1], at=30)    # at a
b20 = !is.adjacent.active(anet, iov[7,2], iov[7,1], at=30)    # at b
b21 = !is.adjacent.active(anet, iov[7,2], iov[7,1], at=25)    # at M
b22 = !is.adjacent.active(anet, iov[3,2], iov[3,1], at=10)   # at M
b23 = !is.adjacent.active(anet, iov[4,2], iov[4,1], at=25)   # at M

# tests that should return false for point queries b/c of non-activity
b24 = !is.adjacent.active(anet, iov[4,1], iov[4,2], at=-Inf, active.default=F)  # at -Inf
b25 = !is.adjacent.active(anet, iov[6,1], iov[6,2], at=-Inf, active.default=F)  # at -Inf
b26 = !is.adjacent.active(anet, iov[3,1], iov[3,2], at=Inf, active.default=F)   # at Inf
b27 = !is.adjacent.active(anet, iov[5,1], iov[5,2], at=Inf, active.default=F)   # at Inf
b28 = !is.adjacent.active(anet, iov[3,1], iov[3,2], at=20, active.default=F)    # at b
b29 = !is.adjacent.active(anet, iov[5,1], iov[5,2], at=20, active.default=F)    # at b
b30 = !is.adjacent.active(anet, iov[3,1], iov[3,2], at=25, active.default=F)    # at H
b31 = !is.adjacent.active(anet, iov[7,1], iov[7,2], at=35, active.default=F)    # at H
b32 = !is.adjacent.active(anet, iov[1,1], iov[1,2], at=25, active.default=F)    # at H

b.tests = paste("b", seq(1,32), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  stop(paste("is.ajacent.active is incorrectly activating edges in tests",
             bad.tests))
}


# tests that should return true for interval queries
c1 = is.adjacent.active(anet, iov[1,1], iov[1,2], -Inf, Inf)  # over (-Inf,Inf)
c2 = is.adjacent.active(anet, iov[2,1], iov[2,2], -Inf, Inf)  # over (-Inf,Inf)
c3 = is.adjacent.active(anet, iov[3,1], iov[3,2], -Inf, 20)   # over (-Inf, b)
c4 = is.adjacent.active(anet, iov[3,1], iov[3,2], -Inf, 30)   # over (-Inf, H)
c5 = is.adjacent.active(anet, iov[4,1], iov[4,2], 0, Inf)     # over (L, Inf)
c6 = is.adjacent.active(anet, iov[4,1], iov[4,2], 10, Inf)    # over (a, Inf)
c7 = is.adjacent.active(anet, iov[5,1], iov[5,2], 10,20)      # over (a, b)
c8 = is.adjacent.active(anet, iov[7,1], iov[7,2], 10,30)      # over (a, b)
c9 = is.adjacent.active(anet, iov[6,1], iov[6,2], 10,10)      # over (a, a)
c10 = is.adjacent.active(anet, iov[7,1], iov[7,2], 0, 40)     # over (L, H)

# tests that should return false for interval queries b/c of non-adjacency
c11 = !is.adjacent.active(anet, iov[1,2], iov[1,1], -Inf, Inf)  # over (-Inf,Inf)
c12 = !is.adjacent.active(anet, iov[2,2], iov[2,1], -Inf, Inf)  # over (-Inf,Inf)
c13 = !is.adjacent.active(anet, iov[3,2], iov[3,1], -Inf, 20)   # over (-Inf, b)
c14 = !is.adjacent.active(anet, iov[3,2], iov[3,1], -Inf, 30)   # over (-Inf, H)
c15 = !is.adjacent.active(anet, iov[4,2], iov[4,1], 0, Inf)     # over (L, Inf)
c16 = !is.adjacent.active(anet, iov[4,2], iov[4,1], 10, Inf)    # over (a, Inf)
c17 = !is.adjacent.active(anet, iov[5,2], iov[5,1], 10,20)      # over (a, b)
c18 = !is.adjacent.active(anet, iov[7,2], iov[7,1], 10,30)      # over (a, b)
c19 = !is.adjacent.active(anet, iov[6,2], iov[6,1], 10,10)      # over (a, a)
c20 = !is.adjacent.active(anet, iov[7,2], iov[7,1], 0, 40)      # over (L, H)

# tests that should return false for interval queries b/c of non-activity
c21 = !is.adjacent.active(anet, iov[1,1], iov[1,2], -Inf, Inf, active.default=F)  # over (-Inf,Inf)
c22 = !is.adjacent.active(anet, iov[8,1], iov[8,2], -Inf, Inf, active.default=F)  # over (-Inf,Inf)
c23 = !is.adjacent.active(anet, iov[3,1], iov[3,2], -Inf, 30, rule="all", active.default=F)  # over (-Inf, H)
c24 = !is.adjacent.active(anet, iov[4,1], iov[4,2], 0, Inf, rule="all", active.default=F)    # over (L, Inf)
c25 = !is.adjacent.active(anet, iov[5,1], iov[5,2], 10, Inf, rule="all", active.default=F)   # over (a, Inf)
c26 = !is.adjacent.active(anet, iov[5,1], iov[5,2], 20,20, active.default=F)                 # over (b, b)
c27 = !is.adjacent.active(anet, iov[5,1], iov[5,2], 0, 5, active.default=F)                  # over (L, L)
c28 = !is.adjacent.active(anet, iov[7,1], iov[7,2], 0,35, rule="all", active.default=F)      # over (L, H)

c.tests = paste("c", seq(1,28), sep="")
c.results= sapply(c.tests, function(x){eval(parse(text=x))})
if(any(!c.results)){
  bad.tests = paste("c", which(!b.results), sep="", collapse=" ")
  stop(paste("activate.edges is incorrectly activating edges in tests",
             bad.tests))
}

expect_false(is.adjacent.active(network.initialize(0),1))


# ----- tests for spells.overlap ------

  
# |----|
# |----|
expect_true(spells.overlap(c(1,2),c(1,2)))

#   |----|
# |----|
expect_true(spells.overlap(c(1,3),c(0,2)))

# |----|
#   |----|
expect_true(spells.overlap(c(0,2),c(1,3)))

# |----|
#  |--|
expect_true(spells.overlap(c(0,3),c(1,2)))

#  |--|
# |----|
expect_true(spells.overlap(c(1,2),c(0,3)))

#  |--|
# |---|
expect_true(spells.overlap(c(1,3),c(0,3)))

# |---|
#  |--|
expect_true(spells.overlap(c(0,3),c(1,3)))

# |--|
# |---|
expect_true(spells.overlap(c(0,2),c(0,3)))

# |---|
# |--|
expect_true(spells.overlap(c(0,3),c(0,2)))


# null spell
# valid spell
expect_false(spells.overlap(c(-Inf,-Inf),c(1,3)))
expect_false(spells.overlap(c(Inf,Inf),c(1,3)))

# valid spell
# null spell
expect_false(spells.overlap(c(1,3),c(-Inf,-Inf)))
expect_false(spells.overlap(c(1,3),c(Inf,Inf)))

# |---|
#       |--|
expect_false(spells.overlap(c(0,2),c(3,5)))

#      |---|
# |--|
expect_false(spells.overlap(c(3,5),c(0,2)))

#    |---|
# |--|
expect_false(spells.overlap(c(3,5),c(0,3)))

# |---|
#     |--|
expect_false(spells.overlap(c(0,3),c(3,5)))

# point queries

# |
#  |--|
expect_false(spells.overlap(c(0,0),c(1,3)))

#  |
#  |--|
expect_true(spells.overlap(c(1,1),c(1,3)))

#   |
#  |--|
expect_true(spells.overlap(c(2,2),c(1,3)))

#     |
#  |--|
expect_false(spells.overlap(c(3,3),c(1,3)))

#       |
#  |--|
expect_false(spells.overlap(c(5,5),c(1,3)))

# |
#  |
expect_false(spells.overlap(c(1,1),c(2,2)))

#  |
#  |
expect_true(spells.overlap(c(1,1),c(1,1)))

#   |
#  |
expect_false(spells.overlap(c(2,2),c(1,1)))

#  |--|
# |
expect_false(spells.overlap(c(1,3),c(0,0)))

#  |--|
#  |
expect_true(spells.overlap(c(1,3),c(1,1)))

#  |--|
#   |
expect_true(spells.overlap(c(1,3),c(2,2)))

#  |--|
#     |
expect_false(spells.overlap(c(1,3),c(3,3)))

#  |--|
#       |
expect_false(spells.overlap(c(1,3),c(4,4)))

# ------ tests for spells.hit ----

expect_error(spells.hit(c(1,2), matrix(ncol=2,nrow=0)),regexp = 'must have at least one row')
expect_error(spells.hit(c(1,2,3), c(1,2)),regexp = 'must have exactly two elements')
expect_error(spells.hit(c(1,2), matrix(ncol=3,nrow=1)),regexp = 'must have exactly two columns')

expect_equal(spells.hit(c(1,2), matrix(c(3,5,
                                         6,7),ncol=2,byrow=TRUE)),-1)
expect_equal(spells.hit(c(1,4), matrix(c(3,5,
                                         6,7),ncol=2,byrow=TRUE)),1)
expect_equal(spells.hit(c(5,5), matrix(c(3,5,
                                         6,7),ncol=2,byrow=TRUE)),-1)
expect_equal(spells.hit(c(1,7), matrix(c(3,5,
                                         6,7),ncol=2,byrow=TRUE)),1)
expect_equal(spells.hit(c(5.5,7), matrix(c(3,5,
                                         6,7),ncol=2,byrow=TRUE)),2)


cat("ok\n")








