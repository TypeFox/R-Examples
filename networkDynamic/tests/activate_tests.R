#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2013 the statnet development team
######################################################################


########################################################################
# This file contains the testing suite for the "activate" methods, i.e.:
#      activate.edges        deactivate.edges
#      activate.vertices     deactivate.vertices     
#########################################################################

require(networkDynamic)
require(testthat)


#------------------- ACTIVATE.EDGES TESTS ------------------------
# Notes:
#  --except for the error-checking code, the rest of the tests
#    are meant to be exhaustive
#-----------------------------------------------------------------

cat("testing activate.edges ... ")
anet <- network.initialize(100)
set.seed(10)
heads <- sample(100, 150, replace=TRUE)
set.seed(25)
tails <- sample(100, 150, replace=TRUE)
add.edges(anet, tails, heads)
anet.copy <- anet

# proper behavior for bad inputs?
#  -- bad network input
a1 = try(activate.edges(3, e=1), T)                    # properly handeled
#  -- bad at times
a2 = try(activate.edges(anet, at="a", e=2), T)         # properly handeled
#  -- bad onset
a3 = try(activate.edges(anet, "a", e=3), T)            # properly handeled
a4 = try(activate.edges(anet, NULL, 3, e=4), T)        # properly handeled
# -- bad terminus
a5 = try(activate.edges(anet, 3, "b", e=5), T)         # properly handeled
# -- bad length
a6 = try(activate.edges(anet, 3, length=-9, e=6), T)   # properly handeled
a7 = try(activate.edges(anet, 3, length="r", e=7), T)  # properly handeled
# -- bad edges
a8 = try(activate.edges(anet, 3, 10, e=174), T)        # properly handeled
a9 = try(activate.edges(anet, 3, 10, e=-2), T)         # properly handeled
a10 = try(activate.edges(anet, 3, 10, e=NULL), T)      # properly handeled
a11 = try(activate.edges(anet, 3, 10, e="hello"), T)   # properly handeled
# -- bad onset & terminus combo
a12 = try(activate.edges(anet, 10, 3, e=8), T)         # properly handeled
# -- not fully specified intervals
a13 = try(activate.edges(anet, 10, e=9), T)            # properly handeled
a14 = try(activate.edges(anet, terminus=10, e=10), T)  # properly handeled
a15 = try(activate.edges(anet, length=10, e=11), T)    # properly handeled


# good input
anet <- anet.copy
# for edges lacking an active attribute
activate.edges(anet, Inf, Inf, e=1)         # add (Inf, Inf) spell
b1 = is.null(anet$mel[[1]]$atl$active)
activate.edges(anet, -Inf, -Inf, e=2)       # add (-Inf, -Inf) spell
b2 = is.null(anet$mel[[2]]$atl$active)
activate.edges(anet, -Inf, Inf, e=3)        # add (-Inf, Inf) spell
b3 = identical(anet$mel[[3]]$atl$active,
              matrix(c(-Inf, Inf), 1,2))
activate.edges(anet, -Inf, Inf, e=4)        # add (-Inf, Inf) spell
b4 = identical(anet$mel[[4]]$atl$active,
              matrix(c(-Inf, Inf), 1,2))
activate.edges(anet, -Inf, 10,  e=5)        # add (-Inf, b) spell
b5 = identical(anet$mel[[5]]$atl$active,
              matrix(c(-Inf, 10), 1,2))
activate.edges(anet, 0, Inf, e=6)           # add (a, Inf) spell
b6 = identical(anet$mel[[6]]$atl$active,
              matrix(c(0, Inf), 1,2))
activate.edges(anet, -10, Inf, e=7)         # add (a, Inf) spell
b7 = identical(anet$mel[[7]]$atl$active,
              matrix(c(-10, Inf), 1,2))
activate.edges(anet, 10, 20, e=8)           # add (a, b) spell
b8 = identical(anet$mel[[8]]$atl$active,
              matrix(c(10, 20), 1,2))
activate.edges(anet, 10, 10, e=9)           # add (a, a) spell
b8 = identical(anet$mel[[9]]$atl$active,
              matrix(c(10, 10), 1,2))

b.tests = paste("b", seq(1,8), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  stop(paste("activate.edges is incorrectly activating edges in tests",
             bad.tests))
}


# for edges already having an active attribute
anet <- anet.copy

# Notes:
#    - a, b are upper, lower boundary points
#    - L is a finite number lower than a
#    - M is a finite number in between a and b
#    - H is a finite number higher than b
#    -Gij is a finite number that is in the gap b/t intervals i and j

# tests for activation at a point
for (i in 1:3)
  anet$mel[[i]]$atl$active <- matrix(c(Inf,Inf), 1,2)
for (i in 4:6)
  anet$mel[[i]]$atl$active <- matrix(c(-Inf,-Inf), 1,2)
activate.edges(anet, -Inf,  Inf, e=7:9)
activate.edges(anet, -Inf,   20, e=10:13)
activate.edges(anet,   10,  Inf, e=14:17)
activate.edges(anet,   10,   20, e=18:24)
activate.edges(anet,   10,   10, e=25:26)

activate.edges(anet, at=c(-Inf, Inf, 0, -Inf, Inf,  0, -Inf, Inf,  0, -Inf, 20,   0, 30, Inf, 10, 30,  0, -Inf, Inf, 10, 20, 15, 30,  0, 10, 20),
                    e = c(   1,   2, 3,    4,   5,  6,    7,   8,  9,   10, 11,  12, 13,  14, 15, 16, 17,   18,  19, 20, 21, 22, 23, 24, 25, 26))  
c1  = identical(anet$mel[[1]]$atl$active,      # add -Inf to  (Inf, Inf)
      matrix(c(Inf, Inf),1,2))    
c2  = identical(anet$mel[[2]]$atl$active,      # add Inf  to  (Inf, Inf)
      matrix(c(Inf, Inf),1,2))    
c3  = identical(anet$mel[[3]]$atl$active,      # add M    to  (Inf, Inf)
      matrix(c(0, 0),1,2))  
c4  = identical(anet$mel[[4]]$atl$active,      # add -Inf to  (-Inf, -Inf)
      matrix(c(-Inf, -Inf),1,2))    
c5  = identical(anet$mel[[5]]$atl$active,      # add Inf  to  (-Inf, -Inf)
      matrix(c(-Inf, -Inf),1,2))    
c6  = identical(anet$mel[[6]]$atl$active,      # add M    to  (-Inf, -Inf)
      matrix(c(0, 0),1,2))    
c7  = identical(anet$mel[[7]]$atl$active,      # add -Inf to  (-Inf, Inf)
      matrix(c(-Inf, Inf),1,2))    
c8  = identical(anet$mel[[8]]$atl$active,      # add Inf  to  (-Inf, Inf)
      matrix(c(-Inf, Inf),1,2))    
c9  = identical(anet$mel[[9]]$atl$active,      # add M    to  (-Inf, Inf)
      matrix(c(-Inf, Inf),1,2))    
c10 = identical(anet$mel[[10]]$atl$active,     # add -Inf to  (-Inf, b)
      matrix(c(-Inf, 20),1,2))      
c11 = identical(anet$mel[[11]]$atl$active,     # add b    to  (-Inf, b)
      matrix(c(-Inf,20,20,20),2,2))    
c12 = identical(anet$mel[[12]]$atl$active,     # add M    to  (-Inf, b)
      matrix(c(-Inf, 20),1,2))    
c13 = identical(anet$mel[[13]]$atl$active,     # add H    to  (-Inf, b)
      matrix(c(-Inf,30,20,30),2,2))      
c14 = identical(anet$mel[[14]]$atl$active,     # add Inf  to  (a, Inf)
      matrix(c(10, Inf),1,2))    
c15 = identical(anet$mel[[15]]$atl$active,     # add a    to  (a, Inf)
      matrix(c(10, Inf),1,2))    
c16 = identical(anet$mel[[16]]$atl$active,     # add M    to  (a, Inf)
      matrix(c(10, Inf),1,2))    
c17 = identical(anet$mel[[17]]$atl$active,     # add L    to  (a, Inf)
      matrix(c(0,10,0,Inf),2,2))    
c18 = identical(anet$mel[[18]]$atl$active,     # add -Inf to  (a,b)
      matrix(c(10,20),1,2))    
c19 = identical(anet$mel[[19]]$atl$active,     # add Inf  to  (a,b)
      matrix(c(10, 20),1,2))    
c20 = identical(anet$mel[[20]]$atl$active,     # add a    to  (a,b)
      matrix(c(10,20),1,2))    
c21 = identical(anet$mel[[21]]$atl$active,     # add b    to  (a,b)
      matrix(c(10,20,20,20),2,2))    
c22 = identical(anet$mel[[22]]$atl$active,     # add M    to  (a,b)
      matrix(c(10, 20),1,2))    
c23 = identical(anet$mel[[23]]$atl$active,     # add H    to  (a,b)
      matrix(c(10,30,20,30),2,2))    
c24 = identical(anet$mel[[24]]$atl$active,     # add L    to  (a,b)    
      matrix(c(0,10,0,20),2,2))  
c25 = identical(anet$mel[[25]]$atl$active,     # add a    to  (a,a)    
      matrix(c(10,10),1,2))  
c26 = identical(anet$mel[[26]]$atl$active,     # add H    to  (a,a)    
      matrix(c(10,20,10,20),2,2))  
  
c.tests = paste("c", seq(1,26), sep="")
c.results= sapply(c.tests, function(x){eval(parse(text=x))})
if(any(!c.results)){
  bad.tests = paste("c", which(!c.results), sep="", collapse=" ")
  stop(paste("activate.edges is incorrectly activating edges in tests",
             bad.tests))
}


anet <- anet.copy
for (i in 1:7)
  anet$mel[[i]]$atl$active <- matrix(c(Inf,Inf), 1,2)
activate.edges(anet, c(Inf, -Inf, -Inf, -Inf,  10, 10, 10),
                     c(Inf, -Inf,  Inf,   10, Inf, 20, 10),
                 e = c(  1,    2,    3,    4,   5,  6,  7))
d1  = identical(anet$mel[[1]]$atl$active,   # add (Inf, Inf)   to (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d2  = identical(anet$mel[[2]]$atl$active,   # add (-Inf, -Inf) to (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d3  = identical(anet$mel[[3]]$atl$active,   # add (-Inf, Inf)  to (Inf, Inf) spell
      matrix(c(-Inf, Inf), 1,2))
d4  = identical(anet$mel[[4]]$atl$active,   # add (-Inf, b)    to (Inf, Inf) spell
      matrix(c(-Inf, 10), 1,2))
d5  = identical(anet$mel[[5]]$atl$active,   # add (a, Inf)     to (Inf, Inf) spell
      matrix(c(10, Inf), 1,2))
d6  = identical(anet$mel[[6]]$atl$active,   # add (a, b)       to (Inf, Inf) spell
      matrix(c(10, 20), 1,2))
d7  = identical(anet$mel[[7]]$atl$active,   # add (a, a)       to (Inf, Inf) spell
      matrix(c(10, 10), 1,2))


for (i in 8:14)
  anet$mel[[i]]$atl$active <- matrix(c(-Inf, -Inf), 1,2)
activate.edges(anet, c(Inf, -Inf, -Inf, -Inf,  10, 10, 20),
                     c(Inf, -Inf,  Inf,   10, Inf, 20, 20),
                 e = c(  8,    9,   10,  11,   12, 13, 14))
d8  = identical(anet$mel[[8]]$atl$active,   # add (Inf, Inf)   to (-Inf, -Inf) spell
      matrix(c(-Inf, -Inf), 1,2))
d9  = identical(anet$mel[[9]]$atl$active,   # add (-Inf, -Inf) to (-Inf, -Inf) spell
      matrix(c(-Inf, -Inf), 1,2))
d10  = identical(anet$mel[[10]]$atl$active, # add (-Inf, Inf)  to (-Inf, -Inf) spell
      matrix(c(-Inf, Inf), 1,2))
d11 = identical(anet$mel[[11]]$atl$active,  # add (-Inf, H)    to (-Inf, -Inf) spell
      matrix(c(-Inf, 10), 1,2))
d12 = identical(anet$mel[[12]]$atl$active,  # add (H, Inf)     to (-Inf, -Inf) spell
      matrix(c(10,  Inf), 1,2, byrow=T))
d13 = identical(anet$mel[[13]]$atl$active,  # add (H1, H2)       to -(Inf, -Inf) spell
      matrix(c(10,  20), 1,2, byrow=T))
d14 = identical(anet$mel[[14]]$atl$active,  # add (H1, H1)       to -(Inf, -Inf) spell
      matrix(c(20,  20), 1,2, byrow=T))


activate.edges(anet, -Inf, 10, e=seq(15,32))
activate.edges(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  20,  0, 10, 10,  0, -5,  0, 20, 20, 10),
                     c(Inf, -Inf,  Inf,    0,   10,   20, Inf, Inf, Inf, 10, 10, 20,  0,  0, 20, 30, 20, 10),
                 e = c( 15,   16,   17,   18,   19,   20,  21,  22,  23, 24, 25, 26, 27, 28, 29, 30, 31, 32))
d15 = identical(anet$mel[[15]]$atl$active,   # add (Inf, Inf)   to (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d16 = identical(anet$mel[[16]]$atl$active,   # add (-Inf, -Inf) to (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d17 = identical(anet$mel[[17]]$atl$active,   # add (-Inf, Inf)  to (-Inf, b) spell
      matrix(c(-Inf, Inf), 1,2))
d18 = identical(anet$mel[[18]]$atl$active,   # add (-Inf, M)    to (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d19 = identical(anet$mel[[19]]$atl$active,   # add (-Inf, b)    to (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2)) 
d20 = identical(anet$mel[[20]]$atl$active,   # add (-Inf, H)    to (-Inf, b) spell
      matrix(c(-Inf, 20), 1,2))
d21 = identical(anet$mel[[21]]$atl$active,   # add (M, Inf)     to (-Inf, b) spell
      matrix(c(-Inf, Inf), 1,2))
d22 = identical(anet$mel[[22]]$atl$active,   # add (b, Inf)     to (-Inf, b) spell
      matrix(c(-Inf, Inf), 1,2))
d23 = identical(anet$mel[[23]]$atl$active,   # add (H, Inf)     to (-Inf, b) spell
      matrix(c(-Inf,20,10,Inf),2,2))
d24 = identical(anet$mel[[24]]$atl$active,   # add (M, b)       to (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d25 = identical(anet$mel[[25]]$atl$active,   # add (b, b)       to (-Inf, b) spell
      matrix(c(-Inf, 10,10,10), 2,2))
d26 = identical(anet$mel[[26]]$atl$active,   # add (b, H)       to (-Inf, b) spell
      matrix(c(-Inf, 20), 1,2))
d27 = identical(anet$mel[[27]]$atl$active,   # add (M, M)       to -(Inf, b) spell
      matrix(c(-Inf, 10), 1,2))  
d28 = identical(anet$mel[[28]]$atl$active,   # add (M1, M2)     to -(Inf, b) spell
      matrix(c(-Inf, 10), 1,2))    
d29 = identical(anet$mel[[29]]$atl$active,   # add (M, H)       to -(Inf, b) spell
      matrix(c(-Inf, 20), 1,2))  
d30 = identical(anet$mel[[30]]$atl$active,   # add (H1, H2)     to -(Inf, b) spell
      matrix(c(-Inf,20,10,30),2,2))
d31 = identical(anet$mel[[31]]$atl$active,   # add (H, H)       to -(Inf, b) spell
      matrix(c(-Inf,20,10,20),2,2))
d32 = identical(anet$mel[[32]]$atl$active,   # add (b, b)       to -(Inf, b) spell
      matrix(c(-Inf,10,10,10),2,2))


activate.edges(anet, 10, Inf, e=seq(33,49))
activate.edges(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  20,  0, 10, 10,  0, -5, 20,  0, 20),
                     c(Inf, -Inf,  Inf,    0,   10,   20, Inf, Inf, Inf, 10, 10, 20,  5, 20, 30,  0, 20),
                 e = c( 33,   34,   35,   36,   37,   38,  39,  40,  41, 42, 43, 44, 45, 46, 47, 48, 49))
d33 = identical(anet$mel[[33]]$atl$active,   # add (Inf, Inf)   to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d34 = identical(anet$mel[[34]]$atl$active,   # add (-Inf, -Inf) to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d35 = identical(anet$mel[[35]]$atl$active,   # add (-Inf, Inf)  to (a, Inf) spell
      matrix(c(-Inf, Inf), 1,2))
d36 = identical(anet$mel[[36]]$atl$active,   # add (-Inf, L)    to (a, Inf) spell
      matrix(c(-Inf,10,0,Inf),2,2))
d37 = identical(anet$mel[[37]]$atl$active,   # add (-Inf, a)    to (a, Inf) spell
      matrix(c(-Inf, Inf), 1,2)) 
d38 = identical(anet$mel[[38]]$atl$active,   # add (-Inf, M)    to (a, Inf) spell
      matrix(c(-Inf, Inf), 1,2))
d39 = identical(anet$mel[[39]]$atl$active,   # add (L, Inf)     to (a, Inf) spell
      matrix(c(0, Inf), 1,2))
d40 = identical(anet$mel[[40]]$atl$active,   # add (a, Inf)     to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d41 = identical(anet$mel[[41]]$atl$active,   # add (M, Inf)     to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d42 = identical(anet$mel[[42]]$atl$active,   # add (L, a)       to (a, Inf) spell
      matrix(c(0, Inf), 1,2))
d43 = identical(anet$mel[[43]]$atl$active,   # add (a, a)       to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d44 = identical(anet$mel[[44]]$atl$active,   # add (a, M)       to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d45 = identical(anet$mel[[45]]$atl$active,   # add (L1, L2)     to (a, Inf) spell
      matrix(c(0,10,5,Inf),2,2))
d46 = identical(anet$mel[[46]]$atl$active,   # add (L, M)       to (a, Inf) spell
      matrix(c(-5, Inf), 1,2))    
d47 = identical(anet$mel[[47]]$atl$active,   # add (M1, M2)     to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d48 = identical(anet$mel[[48]]$atl$active,   # add (L, L)       to (a, Inf) spell
      matrix(c(0,10,0,Inf), 2,2))
d49 = identical(anet$mel[[49]]$atl$active,   # add (H, H)       to (a, Inf) spell
      matrix(c(10, Inf), 1,2))


activate.edges(anet, 10, 20, e=seq(50,76))
activate.edges(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  15,  20,  30,  0, 10, 10, 10, 10,  0, 20, 20, -5,  0,  0, 12, 15, 30),
                     c(Inf, -Inf,  Inf,    0,   10,   15,   20,   30, Inf, Inf, Inf, Inf, Inf, 10, 10, 15, 20, 30, 20, 20, 30,  0, 15, 30, 15, 30, 40),
                 e = c( 50,   51,   52,   53,   54,   55,   56,   57,  58,  59,  60,  61,  62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76))
d50 = identical(anet$mel[[50]]$atl$active,   # add (Inf, Inf)   to (a, b) spell
      matrix(c(10,20), 1,2))
d51 = identical(anet$mel[[51]]$atl$active,   # add (-Inf, -Inf) to (a, b) spell
      matrix(c(10,20), 1,2))  
d52 = identical(anet$mel[[52]]$atl$active,   # add (-Inf, Inf)  to (a, b) spell
      matrix(c(-Inf, Inf), 1,2))
d53 = identical(anet$mel[[53]]$atl$active,   # add (-Inf, L)    to (a, b) spell
      matrix(c(-Inf,10,0,20),2,2))
d54 = identical(anet$mel[[54]]$atl$active,   # add (-Inf, a)    to (a, b) spell
      matrix(c(-Inf, 20), 1,2)) 
d55 = identical(anet$mel[[55]]$atl$active,   # add (-Inf, M)    to (a, b) spell
      matrix(c(-Inf, 20), 1,2))
d56 = identical(anet$mel[[56]]$atl$active,   # add (-Inf, b)    to (a, b) spell
      matrix(c(-Inf, 20), 1,2)) 
d57 = identical(anet$mel[[57]]$atl$active,   # add (-Inf, H)    to (a, b) spell
      matrix(c(-Inf, 30), 1,2))
d58 = identical(anet$mel[[58]]$atl$active,   # add (L, Inf)     to (a, b) spell
      matrix(c(0, Inf), 1,2))
d59 = identical(anet$mel[[59]]$atl$active,   # add (a, Inf)     to (a, b) spell
      matrix(c(10, Inf), 1,2))
d60 = identical(anet$mel[[60]]$atl$active,   # add (M, Inf)     to (a, b) spell
      matrix(c(10, Inf), 1,2))
d61 = identical(anet$mel[[61]]$atl$active,   # add (b, Inf)     to (a, b) spell
      matrix(c(10, Inf), 1,2))
d62 = identical(anet$mel[[62]]$atl$active,   # add (H, Inf)     to (a, b) spell
      matrix(c(10,30,20,Inf),2,2))
d63 = identical(anet$mel[[63]]$atl$active,   # add (L, a)       to (a, b) spell
      matrix(c(0, 20), 1,2))
d64 = identical(anet$mel[[64]]$atl$active,   # add (a, a)       to (a, b) spell
      matrix(c(10, 20), 1,2))
d65 = identical(anet$mel[[65]]$atl$active,   # add (a, M)       to (a, b) spell
      matrix(c(10, 20), 1,2))
d66 = identical(anet$mel[[66]]$atl$active,   # add (a, b)       to (a, b) spell
      matrix(c(10, 20), 1,2))
d67 = identical(anet$mel[[67]]$atl$active,   # add (a, H)       to (a, b) spell
      matrix(c(10, 30), 1,2))
d68 = identical(anet$mel[[68]]$atl$active,   # add (L, b)       to (a, b) spell
      matrix(c(0, 20), 1,2))
d69 = identical(anet$mel[[69]]$atl$active,   # add (b, b)       to (a, b) spell
      matrix(c(10, 20,20,20), 2,2))
d70 = identical(anet$mel[[70]]$atl$active,   # add (b, H)       to (a, b) spell
      matrix(c(10, 30), 1,2))
d71 = identical(anet$mel[[71]]$atl$active,   # add (L1, L2)     to (a, b) spell
      matrix(c(-5,10,0,20),2,2))
d72 = identical(anet$mel[[72]]$atl$active,   # add (L, M)       to (a, b) spell
      matrix(c(0, 20), 1,2))    
d73 = identical(anet$mel[[73]]$atl$active,   # add (L, H)       to (a, b) spell
      matrix(c(0, 30), 1,2))    
d74 = identical(anet$mel[[74]]$atl$active,   # add (M1, M2)     to (a, b) spell
      matrix(c(10, 20), 1,2))
d75 = identical(anet$mel[[75]]$atl$active,   # add (M, H)       to (a, b) spell
      matrix(c(10, 30), 1,2))    
d76 = identical(anet$mel[[76]]$atl$active,   # add (H1, H2)     to (a, b) spell
      matrix(c(10,30,20,40),2,2))

activate.edges(anet, -Inf,  10, e=seq(77, 100))
activate.edges(anet,   20,  30, e=seq(77, 100))
activate.edges(anet,   30,  30, e=seq(77, 100))
activate.edges(anet,   40,  50, e=seq(77, 100))
activate.edges(anet,   60, Inf, e=seq(77, 100))
activate.edges(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,  20,  25,  30,  35, 20, 20, 20, 30, 20, 30, 30, 30, 15, 15, 15,  0,   0),
                     c(Inf,  Inf,   20,   25,   30,   35,   55, Inf, Inf, Inf, Inf, 40, 45, 50, 35, 55, 40, 45, 55, 35, 45, 55, 25,  55),
                 e = c( 77,   78,   79,   80,   81,   82,   83,  84,  85,  86,  87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100))
dmat0 = matrix(c(-Inf,20,30,40,60,10,30,30,50,Inf),5,2)
dmat1 = matrix(c(-Inf,30,40,60,30,30,50,Inf),4,2)
dmat2 = matrix(c(-Inf,40,60,35,50,Inf),3,2)
dmat3 = matrix(c(-Inf,20,60,10,55,Inf),3,2)
dmat4 = matrix(c(-Inf,20,10,Inf),2,2)
dmat5 = matrix(c(-Inf,20,30,35,10,30,30,Inf),4,2)
dmat6 = matrix(c(-Inf,20,60,10,50,Inf),3,2)
dmat7 = matrix(c(-Inf,20,60,10,55,Inf),3,2)
dmat8 = matrix(c(-Inf,15,40,60,10,35,50,Inf),4,2)
dmat9 = matrix(c(-Inf,15,60,10,50,Inf),3,2)
dmat10 = matrix(c(-Inf,15,60,10,55,Inf),3,2)
dmat11 = matrix(c(-Inf,20,40,60,10,35,50,Inf),4,2)
dmat12 = matrix(c(-Inf,60,55,Inf),2,2)

d77 = identical(anet$mel[[77]]$atl$active, dmat0)    # add (Inf, Inf)   to set of spells
d78 = identical(anet$mel[[78]]$atl$active,           # add (-Inf, Inf)  to set of spells
      matrix(c(-Inf, Inf), 1,2))
d79 = identical(anet$mel[[79]]$atl$active, dmat1)    # add (-Inf, a2)   to set of spells
d80 = identical(anet$mel[[80]]$atl$active, dmat1)    # add (-Inf, M2)   to set of spells
d81 = identical(anet$mel[[81]]$atl$active, dmat1)    # add (-Inf, b2)   to set of spells
d82 = identical(anet$mel[[82]]$atl$active, dmat2)    # add (-Inf, G23)  to set of spells
d83 = identical(anet$mel[[83]]$atl$active, dmat12)   # add (-Inf, G34)  to set of spells
d84 = identical(anet$mel[[84]]$atl$active, dmat4)    # add (a2, Inf)    to set of spells
d85 = identical(anet$mel[[85]]$atl$active, dmat4)    # add (M2, Inf)    to set of spells
d86 = identical(anet$mel[[86]]$atl$active, dmat4)    # add (b2, Inf)    to set of spells
d87 = identical(anet$mel[[87]]$atl$active, dmat5)    # add (G23, Inf)   to set of spells
d88 = identical(anet$mel[[88]]$atl$active, dmat6)    # add (a2, a3)     to set of spells
d89 = identical(anet$mel[[89]]$atl$active, dmat6)    # add (a2, M3)     to set of spells
d90 = identical(anet$mel[[90]]$atl$active, dmat6)    # add (a2, b3)     to set of spells
d91 = identical(anet$mel[[91]]$atl$active, dmat11)   # add (a2, G23)    to set of spells
d92 = identical(anet$mel[[92]]$atl$active, dmat3)    # add (a2, G34)    to set of spells
d93 = identical(anet$mel[[93]]$atl$active, dmat6)    # add (b2, a3)     to set of spells
d94 = identical(anet$mel[[94]]$atl$active, dmat6)    # add (b2, M3)     to set of spells
d95 = identical(anet$mel[[95]]$atl$active, dmat7)    # add (b2, G34)    to set of spells
d96 = identical(anet$mel[[96]]$atl$active, dmat8)    # add (G12, G23)   to set of spells
d97 = identical(anet$mel[[97]]$atl$active, dmat9)    # add (G12, M3)    to set of spells
d98 = identical(anet$mel[[98]]$atl$active, dmat10)   # add (G12, G34)   to set of spells
d99 = identical(anet$mel[[99]]$atl$active, dmat1)    # add (M1, M2)     to set of spells
d100 = identical(anet$mel[[100]]$atl$active, dmat12) # add (M1, G34)    to set of spells

               
d.tests = paste("d", seq(1,94), sep="")
d.results= sapply(d.tests, function(x){eval(parse(text=x))})
if(any(!d.results)){
  bad.tests = paste("d", which(!d.results), sep="", collapse=" ")
  stop(paste("activate.edges is incorrectly activating edges in tests",
             bad.tests))
}


# other tests for 'activate.edges'
anet <- anet.copy
mat0=matrix(c(-Inf, Inf),1,2)
# default behavior
activate.edges(anet, e=1)                     # no inputs, no activity matrix
e1 = identical(anet$mel[[1]]$atl$active,mat0)
activate.edges(anet, 3, 10, e=2)               # no inputs, activity matrix
activate.edges(anet, e=2)                     
e2 = identical(anet$mel[[2]]$atl$active,mat0)
anet$mel[[3]]$atl$active <- matrix(c(Inf,Inf), 1,2) # no inputs, matrix of null spell
activate.edges(anet, e=3) 
e3 = identical(anet$mel[[3]]$atl$active,mat0)
#activate.edges(anet, e=2, onset=3)            # only onset
#e2 = identical(anet$mel[[2]]$atl$active,
#  matrix(c(3, Inf),1,2))
#activate.edges(anet, e=3, terminus=3)         # only terminus
#e3 = identical(anet$mel[[3]]$atl$active,
#  matrix(c(-Inf, 3),1,2))
# ignoring 'null' spells
activate.edges(anet, Inf, Inf, e=4)          # single call (Inf, Inf)
e4 = is.null(anet$mel[[4]]$atl$active)
activate.edges(anet, -Inf, -Inf, e=5)        # single call (-Inf, -Inf)
e5 = is.null(anet$mel[[5]]$atl$active)
activate.edges(anet, c(-Inf, -Inf, -Inf,  3, Inf), # mixed call, even lengths
                     c( Inf, -Inf,  10, Inf, Inf),
                     e=seq(6,10))
e6 = is.null(anet$mel[[7]]$atl$active) &
     is.null(anet$mel[[10]]$atl$active)
activate.edges(anet, c(-Inf, -Inf, Inf),     # mixed call, uneven lengths
                     c( Inf, -Inf),
                     e=seq(11,14))
e7 = identical(anet$mel[[11]]$atl$active,
        matrix(c(-Inf, Inf),1,2)) &&
     is.null(anet$mel[[12]]$atl$active) &&
     is.null(anet$mel[[13]]$atl$active) &&
     is.null(anet$mel[[14]]$atl$active)

e.tests = paste("e", seq(1,7), sep="")
e.results= sapply(e.tests, function(x){eval(parse(text=x))})
if(any(!e.results)){
  bad.tests = paste("e", which(!e.results), sep="", collapse=" ")
  stop(paste("activate.edges is incorrectly activating edges in tests",
             bad.tests))
}
# using 'at' rather than 'onset' and 'terminus'
activate.edges(anet, -Inf,  10, e=seq(15, 24))
activate.edges(anet,   20,  30, e=seq(15, 24))
activate.edges(anet,   40, Inf, e=seq(15, 24))
activate.edges(anet, c(Inf, -Inf, 20, 30, 25), c(Inf, -Inf, 20, 30, 25), e=15:19)
activate.edges(anet, at=c(Inf, -Inf, 20, 30, 25), e=20:24)
f1 = identical(anet$mel[[15]]$atl$active, anet$mel[[20]]$atl$active)
f2 = identical(anet$mel[[16]]$atl$active, anet$mel[[21]]$atl$active)
f3 = identical(anet$mel[[17]]$atl$active, anet$mel[[22]]$atl$active)
f4 = identical(anet$mel[[18]]$atl$active, anet$mel[[23]]$atl$active)
f5 = identical(anet$mel[[19]]$atl$active, anet$mel[[24]]$atl$active)
# or 'length' rather than 'terminus'
activate.edges(anet, -Inf,  10, e=seq(25, 43))
activate.edges(anet,   20,  30, e=seq(25, 43))
activate.edges(anet,   40, Inf, e=seq(25, 43))
activate.edges(anet, c(15, 15, 15, 15, 15, 20, 20, 20, 20),
                     c(18, 20, 25, 30, 35, 25, 30, 35, 40),
                   e=c(25, 26, 27, 28, 29, 30, 31, 32, 33))
activate.edges(anet, c(15, 15, 15, 15, 15, 20, 20, 20, 20),
             length =c( 3,  5, 10, 15, 20,  5, 10, 15, 20),
                  e =c(35, 36, 37, 38, 39, 40, 41, 42, 43))
f6  = identical(anet$mel[[25]]$atl$active, anet$mel[[35]]$atl$active)
f7  = identical(anet$mel[[26]]$atl$active, anet$mel[[36]]$atl$active)
f8  = identical(anet$mel[[27]]$atl$active, anet$mel[[37]]$atl$active)
f9  = identical(anet$mel[[28]]$atl$active, anet$mel[[38]]$atl$active)
f10 = identical(anet$mel[[29]]$atl$active, anet$mel[[39]]$atl$active)
f11 = identical(anet$mel[[30]]$atl$active, anet$mel[[40]]$atl$active)
f12 = identical(anet$mel[[31]]$atl$active, anet$mel[[41]]$atl$active)
f13 = identical(anet$mel[[32]]$atl$active, anet$mel[[42]]$atl$active)
f14 = identical(anet$mel[[33]]$atl$active, anet$mel[[43]]$atl$active)
# possible oddities around inclusion of instantaneous points
activate.edges(anet, 10, 20, e=44:48)
activate.edges(anet, at=20, e=44:48)
activate.edges(anet, at=30, e=46:48)
activate.edges(anet, at=10, e=49:50)
activate.edges(anet, at=20, e=49:50)
activate.edges(anet, at=30, e=50)
               
activate.edges(anet, 20, 30, e=44)
f15 = identical(anet$mel[[44]]$atl$active,
        matrix(c(10, 30),1,2))
activate.edges(anet, 15, 30, e=45)
f16 = identical(anet$mel[[45]]$atl$active,
        matrix(c(10, 30),1,2))
activate.edges(anet, 20, 30, e=46)
f17 = identical(anet$mel[[46]]$atl$active,
        matrix(c(10,30,30,30),2,2))
activate.edges(anet, 20, 40, e=47)
f18 = identical(anet$mel[[47]]$atl$active,
        matrix(c(10, 40),1,2))
activate.edges(anet, 30, 40, e=48)
f19 = identical(anet$mel[[48]]$atl$active,
        matrix(c(10,20,30,20,20,40),3,2))
activate.edges(anet, 20,20, e=49)
f20 = identical(anet$mel[[49]]$atl$active,
        matrix(c(10,20,10,20),2,2))
activate.edges(anet, 20,20, e=50)               
f21 = identical(anet$mel[[50]]$atl$active,
        matrix(c(10,20,30,10,20,30),3,2))
               
f.tests = paste("f", seq(1,21), sep="")
f.results= sapply(f.tests, function(x){eval(parse(text=x))})
if(any(!f.results)){
  bad.tests = paste("f", which(!f.results), sep="", collapse=" ")
  stop(paste("activate.edges is incorrectly activating edges in tests",
             bad.tests))
}

# no edges, so should do nothing (but not crash)
expect_equal(network.edgecount(activate.edges(network.initialize(0))),0)

# test activate multiple spells for same edge
test<-network.initialize(3)
test[1,2]<-1
activate.edges(test,onset=c(1,2),terminus=c(2,3),e=c(1,1))
expect_equal(as.numeric(as.data.frame(test)[,1:4]),c(1,3,1,2),info='test that activate edges merges spells if e includes repeated ids')

cat("ok\n")


#------------------- ACTIVATE.VERTICES TESTS ------------------------
# Notes:
#  --except for the error-checking code, the rest of the tests
#    are meant to be exhaustive
#-----------------------------------------------------------------

cat("testing activate.vertices ... ")
anet <- anet.copy

# proper behavior for bad inputs?
#  -- bad network input
a1 = try(activate.vertices(3, v=1), T)                    # properly handeled
#  -- bad at times
a2 = try(activate.vertices(anet, at="a", v=2), T)         # properly handeled
#  -- bad onset
a3 = try(activate.vertices(anet, "a", v=3), T)            # properly handeled
a4 = try(activate.vertices(anet, NULL, 3, v=4), T)        # properly handeled
# -- bad terminus
a5 = try(activate.vertices(anet, 3, "b", v=5), T)         # properly handeled
# -- bad length
a6 = try(activate.vertices(anet, 3, length=-9, v=6), T)   # properly handeled
a7 = try(activate.vertices(anet, 3, length="r", v=7), T)  # properly handeled
# -- bad vertices
a8 = try(activate.vertices(anet, 3, 10, v=174), T)        # properly handeled
a9 = try(activate.vertices(anet, 3, 10, v=-2), T)         # properly handeled
a10 = try(activate.vertices(anet, 3, 10, v=NULL), T)      # properly handeled
a11 = try(activate.vertices(anet, 3, 10, v="hello"), T)   # properly handeled
# -- bad onset & terminus combo
a12 = try(activate.vertices(anet, 10, 3, v=8), T)         # properly handeled
# -- not fully specified intervals
a13 = try(activate.vertices(anet, 10, v=9), T)            # properly handeled
a14 = try(activate.vertices(anet, terminus=10, v=10), T)  # properly handeled
a15 = try(activate.vertices(anet, length=10, v=11), T)    # properly handeled


# good input
anet <- anet.copy
# for vertices lacking an active attribute
activate.vertices(anet, Inf, Inf, v=1)         # add (Inf, Inf) spell
b1 = is.null(anet$val[[1]]$active)
activate.vertices(anet, -Inf, -Inf, v=2)       # add (-Inf, -Inf) spell
b2 = is.null(anet$val[[2]]$active)
activate.vertices(anet, -Inf, Inf, v=3)        # add (-Inf, Inf) spell
b3 = identical(anet$val[[3]]$active,
              matrix(c(-Inf, Inf), 1,2))
activate.vertices(anet, -Inf, Inf, v=4)        # add (-Inf, Inf) spell
b4 = identical(anet$val[[4]]$active,
              matrix(c(-Inf, Inf), 1,2))
activate.vertices(anet, -Inf, 10,  v=5)        # add (-Inf, b) spell
b5 = identical(anet$val[[5]]$active,
              matrix(c(-Inf, 10), 1,2))
activate.vertices(anet, 0, Inf, v=6)           # add (a, Inf) spell
b6 = identical(anet$val[[6]]$active,
              matrix(c(0, Inf), 1,2))
activate.vertices(anet, -10, Inf, v=7)         # add (a, Inf) spell
b7 = identical(anet$val[[7]]$active,
              matrix(c(-10, Inf), 1,2))
activate.vertices(anet, 10, 20, v=8)           # add (a, b) spell
b8 = identical(anet$val[[8]]$active,
              matrix(c(10, 20), 1,2))
activate.vertices(anet, 10, 10, v=9)           # add (a, a) spell
b8 = identical(anet$val[[9]]$active,
              matrix(c(10, 10), 1,2))

b.tests = paste("b", seq(1,8), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  stop(paste("activate.vertices is incorrectly activating vertices in tests",
             bad.tests))
}


# for vertices already having an active attribute
anet <- anet.copy

# Notes:
#    - a, b are upper, lower boundary points
#    - L is a finite number lower than a
#    - M is a finite number in between a and b
#    - H is a finite number higher than b
#    -Gij is a finite number that is in the gap b/t intervals i and j

# tests for activation at a point
for (i in 1:3)
  anet$val[[i]]$active <- matrix(c(Inf,Inf), 1,2)
for (i in 4:6)
  anet$val[[i]]$active <- matrix(c(-Inf,-Inf), 1,2)
activate.vertices(anet, -Inf,  Inf, v=7:9)
activate.vertices(anet, -Inf,   20, v=10:13)
activate.vertices(anet,   10,  Inf, v=14:17)
activate.vertices(anet,   10,   20, v=18:24)
activate.vertices(anet,   10,   10, v=25:26)

activate.vertices(anet, at=c(-Inf, Inf, 0, -Inf, Inf,  0, -Inf, Inf,  0, -Inf, 20,   0, 30, Inf, 10, 30,  0, -Inf, Inf, 10, 20, 15, 30,  0, 10, 20),
                    v = c(   1,   2, 3,    4,   5,  6,    7,   8,  9,   10, 11,  12, 13,  14, 15, 16, 17,   18,  19, 20, 21, 22, 23, 24, 25, 26))  
c1  = identical(anet$val[[1]]$active,      # add -Inf to  (Inf, Inf)
      matrix(c(Inf, Inf),1,2))    
c2  = identical(anet$val[[2]]$active,      # add Inf  to  (Inf, Inf)
      matrix(c(Inf, Inf),1,2))    
c3  = identical(anet$val[[3]]$active,      # add M    to  (Inf, Inf)
      matrix(c(0, 0),1,2))  
c4  = identical(anet$val[[4]]$active,      # add -Inf to  (-Inf, -Inf)
      matrix(c(-Inf, -Inf),1,2))    
c5  = identical(anet$val[[5]]$active,      # add Inf  to  (-Inf, -Inf)
      matrix(c(-Inf, -Inf),1,2))    
c6  = identical(anet$val[[6]]$active,      # add M    to  (-Inf, -Inf)
      matrix(c(0, 0),1,2))    
c7  = identical(anet$val[[7]]$active,      # add -Inf to  (-Inf, Inf)
      matrix(c(-Inf, Inf),1,2))    
c8  = identical(anet$val[[8]]$active,      # add Inf  to  (-Inf, Inf)
      matrix(c(-Inf, Inf),1,2))    
c9  = identical(anet$val[[9]]$active,      # add M    to  (-Inf, Inf)
      matrix(c(-Inf, Inf),1,2))    
c10 = identical(anet$val[[10]]$active,     # add -Inf to  (-Inf, b)
      matrix(c(-Inf, 20),1,2))      
c11 = identical(anet$val[[11]]$active,     # add b    to  (-Inf, b)
      matrix(c(-Inf,20,20,20),2,2))    
c12 = identical(anet$val[[12]]$active,     # add M    to  (-Inf, b)
      matrix(c(-Inf, 20),1,2))    
c13 = identical(anet$val[[13]]$active,     # add H    to  (-Inf, b)
      matrix(c(-Inf,30,20,30),2,2))      
c14 = identical(anet$val[[14]]$active,     # add Inf  to  (a, Inf)
      matrix(c(10, Inf),1,2))    
c15 = identical(anet$val[[15]]$active,     # add a    to  (a, Inf)
      matrix(c(10, Inf),1,2))    
c16 = identical(anet$val[[16]]$active,     # add M    to  (a, Inf)
      matrix(c(10, Inf),1,2))    
c17 = identical(anet$val[[17]]$active,     # add L    to  (a, Inf)
      matrix(c(0,10,0,Inf),2,2))    
c18 = identical(anet$val[[18]]$active,     # add -Inf to  (a,b)
      matrix(c(10,20),1,2))    
c19 = identical(anet$val[[19]]$active,     # add Inf  to  (a,b)
      matrix(c(10, 20),1,2))    
c20 = identical(anet$val[[20]]$active,     # add a    to  (a,b)
      matrix(c(10,20),1,2))    
c21 = identical(anet$val[[21]]$active,     # add b    to  (a,b)
      matrix(c(10,20,20,20),2,2))    
c22 = identical(anet$val[[22]]$active,     # add M    to  (a,b)
      matrix(c(10, 20),1,2))    
c23 = identical(anet$val[[23]]$active,     # add H    to  (a,b)
      matrix(c(10,30,20,30),2,2))    
c24 = identical(anet$val[[24]]$active,     # add L    to  (a,b)    
      matrix(c(0,10,0,20),2,2))  
c25 = identical(anet$val[[25]]$active,     # add a    to  (a,a)    
      matrix(c(10,10),1,2))  
c26 = identical(anet$val[[26]]$active,     # add H    to  (a,a)    
      matrix(c(10,20,10,20),2,2))  
  
c.tests = paste("c", seq(1,26), sep="")
c.results= sapply(c.tests, function(x){eval(parse(text=x))})
if(any(!c.results)){
  bad.tests = paste("c", which(!c.results), sep="", collapse=" ")
  stop(paste("activate.vertices is incorrectly activating vertices in tests",
             bad.tests))
}


anet <- anet.copy
for (i in 1:7)
  anet$val[[i]]$active <- matrix(c(Inf,Inf), 1,2)
activate.vertices(anet, c(Inf, -Inf, -Inf, -Inf,  10, 10, 10),
                     c(Inf, -Inf,  Inf,   10, Inf, 20, 10),
                 v = c(  1,    2,    3,    4,   5,  6,  7))
d1  = identical(anet$val[[1]]$active,   # add (Inf, Inf)   to (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d2  = identical(anet$val[[2]]$active,   # add (-Inf, -Inf) to (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d3  = identical(anet$val[[3]]$active,   # add (-Inf, Inf)  to (Inf, Inf) spell
      matrix(c(-Inf, Inf), 1,2))
d4  = identical(anet$val[[4]]$active,   # add (-Inf, b)    to (Inf, Inf) spell
      matrix(c(-Inf, 10), 1,2))
d5  = identical(anet$val[[5]]$active,   # add (a, Inf)     to (Inf, Inf) spell
      matrix(c(10, Inf), 1,2))
d6  = identical(anet$val[[6]]$active,   # add (a, b)       to (Inf, Inf) spell
      matrix(c(10, 20), 1,2))
d7  = identical(anet$val[[7]]$active,   # add (a, a)       to (Inf, Inf) spell
      matrix(c(10, 10), 1,2))


for (i in 8:14)
  anet$val[[i]]$active <- matrix(c(-Inf, -Inf), 1,2)
activate.vertices(anet, c(Inf, -Inf, -Inf, -Inf,  10, 10, 20),
                     c(Inf, -Inf,  Inf,   10, Inf, 20, 20),
                 v = c(  8,    9,   10,  11,   12, 13, 14))
d8  = identical(anet$val[[8]]$active,   # add (Inf, Inf)   to (-Inf, -Inf) spell
      matrix(c(-Inf, -Inf), 1,2))
d9  = identical(anet$val[[9]]$active,   # add (-Inf, -Inf) to (-Inf, -Inf) spell
      matrix(c(-Inf, -Inf), 1,2))
d10  = identical(anet$val[[10]]$active, # add (-Inf, Inf)  to (-Inf, -Inf) spell
      matrix(c(-Inf, Inf), 1,2))
d11 = identical(anet$val[[11]]$active,  # add (-Inf, H)    to (-Inf, -Inf) spell
      matrix(c(-Inf, 10), 1,2))
d12 = identical(anet$val[[12]]$active,  # add (H, Inf)     to (-Inf, -Inf) spell
      matrix(c(10,  Inf), 1,2, byrow=T))
d13 = identical(anet$val[[13]]$active,  # add (H1, H2)       to -(Inf, -Inf) spell
      matrix(c(10,  20), 1,2, byrow=T))
d14 = identical(anet$val[[14]]$active,  # add (H1, H1)       to -(Inf, -Inf) spell
      matrix(c(20,  20), 1,2, byrow=T))


activate.vertices(anet, -Inf, 10, v=seq(15,32))
activate.vertices(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  20,  0, 10, 10,  0, -5,  0, 20, 20, 10),
                     c(Inf, -Inf,  Inf,    0,   10,   20, Inf, Inf, Inf, 10, 10, 20,  0,  0, 20, 30, 20, 10),
                 v = c( 15,   16,   17,   18,   19,   20,  21,  22,  23, 24, 25, 26, 27, 28, 29, 30, 31, 32))
d15 = identical(anet$val[[15]]$active,   # add (Inf, Inf)   to (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d16 = identical(anet$val[[16]]$active,   # add (-Inf, -Inf) to (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d17 = identical(anet$val[[17]]$active,   # add (-Inf, Inf)  to (-Inf, b) spell
      matrix(c(-Inf, Inf), 1,2))
d18 = identical(anet$val[[18]]$active,   # add (-Inf, M)    to (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d19 = identical(anet$val[[19]]$active,   # add (-Inf, b)    to (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2)) 
d20 = identical(anet$val[[20]]$active,   # add (-Inf, H)    to (-Inf, b) spell
      matrix(c(-Inf, 20), 1,2))
d21 = identical(anet$val[[21]]$active,   # add (M, Inf)     to (-Inf, b) spell
      matrix(c(-Inf, Inf), 1,2))
d22 = identical(anet$val[[22]]$active,   # add (b, Inf)     to (-Inf, b) spell
      matrix(c(-Inf, Inf), 1,2))
d23 = identical(anet$val[[23]]$active,   # add (H, Inf)     to (-Inf, b) spell
      matrix(c(-Inf,20,10,Inf),2,2))
d24 = identical(anet$val[[24]]$active,   # add (M, b)       to (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d25 = identical(anet$val[[25]]$active,   # add (b, b)       to (-Inf, b) spell
      matrix(c(-Inf, 10,10,10), 2,2))
d26 = identical(anet$val[[26]]$active,   # add (b, H)       to (-Inf, b) spell
      matrix(c(-Inf, 20), 1,2))
d27 = identical(anet$val[[27]]$active,   # add (M, M)       to -(Inf, b) spell
      matrix(c(-Inf, 10), 1,2))  
d28 = identical(anet$val[[28]]$active,   # add (M1, M2)     to -(Inf, b) spell
      matrix(c(-Inf, 10), 1,2))    
d29 = identical(anet$val[[29]]$active,   # add (M, H)       to -(Inf, b) spell
      matrix(c(-Inf, 20), 1,2))  
d30 = identical(anet$val[[30]]$active,   # add (H1, H2)     to -(Inf, b) spell
      matrix(c(-Inf,20,10,30),2,2))
d31 = identical(anet$val[[31]]$active,   # add (H, H)       to -(Inf, b) spell
      matrix(c(-Inf,20,10,20),2,2))
d32 = identical(anet$val[[32]]$active,   # add (b, b)       to -(Inf, b) spell
      matrix(c(-Inf,10,10,10),2,2))


activate.vertices(anet, 10, Inf, v=seq(33,49))
activate.vertices(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  20,  0, 10, 10,  0, -5, 20,  0, 20),
                     c(Inf, -Inf,  Inf,    0,   10,   20, Inf, Inf, Inf, 10, 10, 20,  5, 20, 30,  0, 20),
                 v = c( 33,   34,   35,   36,   37,   38,  39,  40,  41, 42, 43, 44, 45, 46, 47, 48, 49))
d33 = identical(anet$val[[33]]$active,   # add (Inf, Inf)   to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d34 = identical(anet$val[[34]]$active,   # add (-Inf, -Inf) to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d35 = identical(anet$val[[35]]$active,   # add (-Inf, Inf)  to (a, Inf) spell
      matrix(c(-Inf, Inf), 1,2))
d36 = identical(anet$val[[36]]$active,   # add (-Inf, L)    to (a, Inf) spell
      matrix(c(-Inf,10,0,Inf),2,2))
d37 = identical(anet$val[[37]]$active,   # add (-Inf, a)    to (a, Inf) spell
      matrix(c(-Inf, Inf), 1,2)) 
d38 = identical(anet$val[[38]]$active,   # add (-Inf, M)    to (a, Inf) spell
      matrix(c(-Inf, Inf), 1,2))
d39 = identical(anet$val[[39]]$active,   # add (L, Inf)     to (a, Inf) spell
      matrix(c(0, Inf), 1,2))
d40 = identical(anet$val[[40]]$active,   # add (a, Inf)     to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d41 = identical(anet$val[[41]]$active,   # add (M, Inf)     to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d42 = identical(anet$val[[42]]$active,   # add (L, a)       to (a, Inf) spell
      matrix(c(0, Inf), 1,2))
d43 = identical(anet$val[[43]]$active,   # add (a, a)       to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d44 = identical(anet$val[[44]]$active,   # add (a, M)       to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d45 = identical(anet$val[[45]]$active,   # add (L1, L2)     to (a, Inf) spell
      matrix(c(0,10,5,Inf),2,2))
d46 = identical(anet$val[[46]]$active,   # add (L, M)       to (a, Inf) spell
      matrix(c(-5, Inf), 1,2))    
d47 = identical(anet$val[[47]]$active,   # add (M1, M2)     to (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d48 = identical(anet$val[[48]]$active,   # add (L, L)       to (a, Inf) spell
      matrix(c(0,10,0,Inf), 2,2))
d49 = identical(anet$val[[49]]$active,   # add (H, H)       to (a, Inf) spell
      matrix(c(10, Inf), 1,2))


activate.vertices(anet, 10, 20, v=seq(50,76))
activate.vertices(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  15,  20,  30,  0, 10, 10, 10, 10,  0, 20, 20, -5,  0,  0, 12, 15, 30),
                     c(Inf, -Inf,  Inf,    0,   10,   15,   20,   30, Inf, Inf, Inf, Inf, Inf, 10, 10, 15, 20, 30, 20, 20, 30,  0, 15, 30, 15, 30, 40),
                 v = c( 50,   51,   52,   53,   54,   55,   56,   57,  58,  59,  60,  61,  62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76))
d50 = identical(anet$val[[50]]$active,   # add (Inf, Inf)   to (a, b) spell
      matrix(c(10,20), 1,2))
d51 = identical(anet$val[[51]]$active,   # add (-Inf, -Inf) to (a, b) spell
      matrix(c(10,20), 1,2))  
d52 = identical(anet$val[[52]]$active,   # add (-Inf, Inf)  to (a, b) spell
      matrix(c(-Inf, Inf), 1,2))
d53 = identical(anet$val[[53]]$active,   # add (-Inf, L)    to (a, b) spell
      matrix(c(-Inf,10,0,20),2,2))
d54 = identical(anet$val[[54]]$active,   # add (-Inf, a)    to (a, b) spell
      matrix(c(-Inf, 20), 1,2)) 
d55 = identical(anet$val[[55]]$active,   # add (-Inf, M)    to (a, b) spell
      matrix(c(-Inf, 20), 1,2))
d56 = identical(anet$val[[56]]$active,   # add (-Inf, b)    to (a, b) spell
      matrix(c(-Inf, 20), 1,2)) 
d57 = identical(anet$val[[57]]$active,   # add (-Inf, H)    to (a, b) spell
      matrix(c(-Inf, 30), 1,2))
d58 = identical(anet$val[[58]]$active,   # add (L, Inf)     to (a, b) spell
      matrix(c(0, Inf), 1,2))
d59 = identical(anet$val[[59]]$active,   # add (a, Inf)     to (a, b) spell
      matrix(c(10, Inf), 1,2))
d60 = identical(anet$val[[60]]$active,   # add (M, Inf)     to (a, b) spell
      matrix(c(10, Inf), 1,2))
d61 = identical(anet$val[[61]]$active,   # add (b, Inf)     to (a, b) spell
      matrix(c(10, Inf), 1,2))
d62 = identical(anet$val[[62]]$active,   # add (H, Inf)     to (a, b) spell
      matrix(c(10,30,20,Inf),2,2))
d63 = identical(anet$val[[63]]$active,   # add (L, a)       to (a, b) spell
      matrix(c(0, 20), 1,2))
d64 = identical(anet$val[[64]]$active,   # add (a, a)       to (a, b) spell
      matrix(c(10, 20), 1,2))
d65 = identical(anet$val[[65]]$active,   # add (a, M)       to (a, b) spell
      matrix(c(10, 20), 1,2))
d66 = identical(anet$val[[66]]$active,   # add (a, b)       to (a, b) spell
      matrix(c(10, 20), 1,2))
d67 = identical(anet$val[[67]]$active,   # add (a, H)       to (a, b) spell
      matrix(c(10, 30), 1,2))
d68 = identical(anet$val[[68]]$active,   # add (L, b)       to (a, b) spell
      matrix(c(0, 20), 1,2))
d69 = identical(anet$val[[69]]$active,   # add (b, b)       to (a, b) spell
      matrix(c(10, 20,20,20), 2,2))
d70 = identical(anet$val[[70]]$active,   # add (b, H)       to (a, b) spell
      matrix(c(10, 30), 1,2))
d71 = identical(anet$val[[71]]$active,   # add (L1, L2)     to (a, b) spell
      matrix(c(-5,10,0,20),2,2))
d72 = identical(anet$val[[72]]$active,   # add (L, M)       to (a, b) spell
      matrix(c(0, 20), 1,2))    
d73 = identical(anet$val[[73]]$active,   # add (L, H)       to (a, b) spell
      matrix(c(0, 30), 1,2))    
d74 = identical(anet$val[[74]]$active,   # add (M1, M2)     to (a, b) spell
      matrix(c(10, 20), 1,2))
d75 = identical(anet$val[[75]]$active,   # add (M, H)       to (a, b) spell
      matrix(c(10, 30), 1,2))    
d76 = identical(anet$val[[76]]$active,   # add (H1, H2)     to (a, b) spell
      matrix(c(10,30,20,40),2,2))

activate.vertices(anet, -Inf,  10, v=seq(77, 100))
activate.vertices(anet,   20,  30, v=seq(77, 100))
activate.vertices(anet,   30,  30, v=seq(77, 100))
activate.vertices(anet,   40,  50, v=seq(77, 100))
activate.vertices(anet,   60, Inf, v=seq(77, 100))
activate.vertices(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,  20,  25,  30,  35, 20, 20, 20, 30, 20, 30, 30, 30, 15, 15, 15,  0,   0),
                     c(Inf,  Inf,   20,   25,   30,   35,   55, Inf, Inf, Inf, Inf, 40, 45, 50, 35, 55, 40, 45, 55, 35, 45, 55, 25,  55),
                 v = c( 77,   78,   79,   80,   81,   82,   83,  84,  85,  86,  87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100))
dmat0 = matrix(c(-Inf,20,30,40,60,10,30,30,50,Inf),5,2)
dmat1 = matrix(c(-Inf,30,40,60,30,30,50,Inf),4,2)
dmat2 = matrix(c(-Inf,40,60,35,50,Inf),3,2)
dmat3 = matrix(c(-Inf,20,60,10,55,Inf),3,2)
dmat4 = matrix(c(-Inf,20,10,Inf),2,2)
dmat5 = matrix(c(-Inf,20,30,35,10,30,30,Inf),4,2)
dmat6 = matrix(c(-Inf,20,60,10,50,Inf),3,2)
dmat7 = matrix(c(-Inf,20,60,10,55,Inf),3,2)
dmat8 = matrix(c(-Inf,15,40,60,10,35,50,Inf),4,2)
dmat9 = matrix(c(-Inf,15,60,10,50,Inf),3,2)
dmat10 = matrix(c(-Inf,15,60,10,55,Inf),3,2)
dmat11 = matrix(c(-Inf,20,40,60,10,35,50,Inf),4,2)
dmat12 = matrix(c(-Inf,60,55,Inf),2,2)

d77 = identical(anet$val[[77]]$active, dmat0)    # add (Inf, Inf)   to set of spells
d78 = identical(anet$val[[78]]$active,           # add (-Inf, Inf)  to set of spells
      matrix(c(-Inf, Inf), 1,2))
d79 = identical(anet$val[[79]]$active, dmat1)    # add (-Inf, a2)   to set of spells
d80 = identical(anet$val[[80]]$active, dmat1)    # add (-Inf, M2)   to set of spells
d81 = identical(anet$val[[81]]$active, dmat1)    # add (-Inf, b2)   to set of spells
d82 = identical(anet$val[[82]]$active, dmat2)    # add (-Inf, G23)  to set of spells
d83 = identical(anet$val[[83]]$active, dmat12)   # add (-Inf, G34)  to set of spells
d84 = identical(anet$val[[84]]$active, dmat4)    # add (a2, Inf)    to set of spells
d85 = identical(anet$val[[85]]$active, dmat4)    # add (M2, Inf)    to set of spells
d86 = identical(anet$val[[86]]$active, dmat4)    # add (b2, Inf)    to set of spells
d87 = identical(anet$val[[87]]$active, dmat5)    # add (G23, Inf)   to set of spells
d88 = identical(anet$val[[88]]$active, dmat6)    # add (a2, a3)     to set of spells
d89 = identical(anet$val[[89]]$active, dmat6)    # add (a2, M3)     to set of spells
d90 = identical(anet$val[[90]]$active, dmat6)    # add (a2, b3)     to set of spells
d91 = identical(anet$val[[91]]$active, dmat11)   # add (a2, G23)    to set of spells
d92 = identical(anet$val[[92]]$active, dmat3)    # add (a2, G34)    to set of spells
d93 = identical(anet$val[[93]]$active, dmat6)    # add (b2, a3)     to set of spells
d94 = identical(anet$val[[94]]$active, dmat6)    # add (b2, M3)     to set of spells
d95 = identical(anet$val[[95]]$active, dmat7)    # add (b2, G34)    to set of spells
d96 = identical(anet$val[[96]]$active, dmat8)    # add (G12, G23)   to set of spells
d97 = identical(anet$val[[97]]$active, dmat9)    # add (G12, M3)    to set of spells
d98 = identical(anet$val[[98]]$active, dmat10)   # add (G12, G34)   to set of spells
d99 = identical(anet$val[[99]]$active, dmat1)    # add (M1, M2)     to set of spells
d100 = identical(anet$val[[100]]$active, dmat12) # add (M1, G34)    to set of spells

               
d.tests = paste("d", seq(1,94), sep="")
d.results= sapply(d.tests, function(x){eval(parse(text=x))})
if(any(!d.results)){
  bad.tests = paste("d", which(!d.results), sep="", collapse=" ")
  stop(paste("activate.vertices is incorrectly activating vertices in tests",
             bad.tests))
}


# other tests for 'activate.vertices'
anet <- anet.copy
mat0=matrix(c(-Inf, Inf),1,2)
# default behavior
activate.vertices(anet, v=1)                     # no inputs, no activity matrix
e1 = identical(anet$val[[1]]$active,mat0)
activate.vertices(anet, 3, 10, v=2)               # no inputs, activity matrix
activate.vertices(anet, v=2)                     
e2 = identical(anet$val[[2]]$active,mat0)
anet$val[[3]]$active <- matrix(c(Inf,Inf), 1,2) # no inputs, matrix of null spell
activate.vertices(anet, v=3) 
e3 = identical(anet$val[[3]]$active,mat0)
#activate.vertices(anet, v=2, onset=3)            # only onset
#e2 = identical(anet$val[[2]]$active,
#  matrix(c(3, Inf),1,2))
#activate.vertices(anet, v=3, terminus=3)         # only terminus
#e3 = identical(anet$val[[3]]$active,
#  matrix(c(-Inf, 3),1,2))
# ignoring 'null' spells
activate.vertices(anet, Inf, Inf, v=4)          # single call (Inf, Inf)
e4 = is.null(anet$val[[4]]$active)
activate.vertices(anet, -Inf, -Inf, v=5)        # single call (-Inf, -Inf)
e5 = is.null(anet$val[[5]]$active)
activate.vertices(anet, c(-Inf, -Inf, -Inf,  3, Inf), # mixed call, even lengths
                     c( Inf, -Inf,  10, Inf, Inf),
                     v=seq(6,10))
e6 = is.null(anet$val[[7]]$active) &
     is.null(anet$val[[10]]$active)
activate.vertices(anet, c(-Inf, -Inf, Inf),     # mixed call, uneven lengths
                     c( Inf, -Inf),
                     v=seq(11,14))
e7 = identical(anet$val[[11]]$active,
        matrix(c(-Inf, Inf),1,2)) &&
     is.null(anet$val[[12]]$active) &&
     is.null(anet$val[[13]]$active) &&
     is.null(anet$val[[14]]$active)

e.tests = paste("e", seq(1,7), sep="")
e.results= sapply(e.tests, function(x){eval(parse(text=x))})
if(any(!e.results)){
  bad.tests = paste("e", which(!e.results), sep="", collapse=" ")
  stop(paste("activate.vertices is incorrectly activating vertices in tests",
             bad.tests))
}
# using 'at' rather than 'onset' and 'terminus'
activate.vertices(anet, -Inf,  10, v=seq(15, 24))
activate.vertices(anet,   20,  30, v=seq(15, 24))
activate.vertices(anet,   40, Inf, v=seq(15, 24))
activate.vertices(anet, c(Inf, -Inf, 20, 30, 25), c(Inf, -Inf, 20, 30, 25), v=15:19)
activate.vertices(anet, at=c(Inf, -Inf, 20, 30, 25), v=20:24)
f1 = identical(anet$val[[15]]$active, anet$val[[20]]$active)
f2 = identical(anet$val[[16]]$active, anet$val[[21]]$active)
f3 = identical(anet$val[[17]]$active, anet$val[[22]]$active)
f4 = identical(anet$val[[18]]$active, anet$val[[23]]$active)
f5 = identical(anet$val[[19]]$active, anet$val[[24]]$active)
# or 'length' rather than 'terminus'
activate.vertices(anet, -Inf,  10, v=seq(25, 43))
activate.vertices(anet,   20,  30, v=seq(25, 43))
activate.vertices(anet,   40, Inf, v=seq(25, 43))
activate.vertices(anet, c(15, 15, 15, 15, 15, 20, 20, 20, 20),
                     c(18, 20, 25, 30, 35, 25, 30, 35, 40),
                   v=c(25, 26, 27, 28, 29, 30, 31, 32, 33))
activate.vertices(anet, c(15, 15, 15, 15, 15, 20, 20, 20, 20),
             length =c( 3,  5, 10, 15, 20,  5, 10, 15, 20),
                  v =c(35, 36, 37, 38, 39, 40, 41, 42, 43))
f6  = identical(anet$val[[25]]$active, anet$val[[35]]$active)
f7  = identical(anet$val[[26]]$active, anet$val[[36]]$active)
f8  = identical(anet$val[[27]]$active, anet$val[[37]]$active)
f9  = identical(anet$val[[28]]$active, anet$val[[38]]$active)
f10 = identical(anet$val[[29]]$active, anet$val[[39]]$active)
f11 = identical(anet$val[[30]]$active, anet$val[[40]]$active)
f12 = identical(anet$val[[31]]$active, anet$val[[41]]$active)
f13 = identical(anet$val[[32]]$active, anet$val[[42]]$active)
f14 = identical(anet$val[[33]]$active, anet$val[[43]]$active)
# possible oddities around inclusion of instantaneous points
activate.vertices(anet, 10, 20, v=44:48)
activate.vertices(anet, at=20, v=44:48)
activate.vertices(anet, at=30, v=46:48)
activate.vertices(anet, at=10, v=49:50)
activate.vertices(anet, at=20, v=49:50)
activate.vertices(anet, at=30, v=50)
               
activate.vertices(anet, 20, 30, v=44)
f15 = identical(anet$val[[44]]$active,
        matrix(c(10, 30),1,2))
activate.vertices(anet, 15, 30, v=45)
f16 = identical(anet$val[[45]]$active,
        matrix(c(10, 30),1,2))
activate.vertices(anet, 20, 30, v=46)
f17 = identical(anet$val[[46]]$active,
        matrix(c(10,30,30,30),2,2))
activate.vertices(anet, 20, 40, v=47)
f18 = identical(anet$val[[47]]$active,
        matrix(c(10, 40),1,2))
activate.vertices(anet, 30, 40, v=48)
f19 = identical(anet$val[[48]]$active,
        matrix(c(10,20,30,20,20,40),3,2))
activate.vertices(anet, 20,20, v=49)
f20 = identical(anet$val[[49]]$active,
        matrix(c(10,20,10,20),2,2))
activate.vertices(anet, 20,20, v=50)               
f21 = identical(anet$val[[50]]$active,
        matrix(c(10,20,30,10,20,30),3,2))
               
f.tests = paste("f", seq(1,21), sep="")
f.results= sapply(f.tests, function(x){eval(parse(text=x))})
if(any(!f.results)){
  bad.tests = paste("f", which(!f.results), sep="", collapse=" ")
  stop(paste("activate.vertices is incorrectly activating vertices in tests",
             bad.tests))
}

# should do nothing, since no verts to activate
expect_equal(is.active(activate.vertices(network.initialize(0),at=1),at=1),logical(0))

# test activating multiple spells per vertex
test <- network.initialize(3)
activate.vertices(test,onset=0:3,terminus=1:4,v=c(1,2,3,1))
expect_equal(get.vertex.activity(test, as.spellList=TRUE)[,1:3],as.data.frame(matrix(c(0,1,1, 3,4,1, 1,2,2, 2,3,3),byrow=TRUE,ncol=3)),check.attributes=FALSE, info="check activating multiple spells on a single vertex")


cat("ok\n")



#------------------- DEACTIVATE.EDGES TESTS ------------------------
# Notes:
#  --except for the error-checking code, the rest of the tests
#    are meant to be exhaustive
#-----------------------------------------------------------------

cat("testing deactivate.edges ... ")
anet <- anet.copy

# proper behavior for bad inputs?
#  -- bad network input
a1 = try(deactivate.edges(3, e=1), T)                    # properly handeled
#  -- bad at times
a2 = try(deactivate.edges(anet, at="a", e=2), T)         # properly handeled
#  -- bad onset
a3 = try(deactivate.edges(anet, "a", e=3), T)            # properly handeled
a4 = try(deactivate.edges(anet, NULL, 3, e=4), T)        # properly handeled
# -- bad terminus
a5 = try(deactivate.edges(anet, 3, "b", e=5), T)         # properly handeled
# -- bad length
a6 = try(deactivate.edges(anet, 3, length=-9, e=6), T)   # properly handeled
a7 = try(deactivate.edges(anet, 3, length="r", e=7), T)  # properly handeled
# -- bad edges
a8 = try(deactivate.edges(anet, 3, 10, e=174), T)        # properly handeled
a9 = try(deactivate.edges(anet, 3, 10, e=-2), T)         # properly handeled
a10 = try(deactivate.edges(anet, 3, 10, e=NULL), T)      # properly handeled
a11 = try(deactivate.edges(anet, 3, 10, e="hello"), T)   # properly handeled
# -- bad onset & terminus combo
a12 = try(deactivate.edges(anet, 10, 3, e=8), T)         # properly handeled
# -- not fully specified intervals
a13 = try(deactivate.edges(anet, 10, e=9), T)            # properly handeled
a14 = try(deactivate.edges(anet, terminus=10, e=10), T)  # properly handeled
a15 = try(deactivate.edges(anet, length=10, e=11), T)    # properly handeled


# good input
anet <- anet.copy
mat0 <- matrix(c(Inf, Inf), 1,2)
# for edges lacking an active attribute
deactivate.edges(anet, Inf, Inf, e=1)         # deact (Inf, Inf) spell
b1 = identical(anet$mel[[1]]$atl$active,mat0)
deactivate.edges(anet, -Inf, -Inf, e=2)       # deact (-Inf, -Inf) spell
b2 = identical(anet$mel[[2]]$atl$active,mat0)
deactivate.edges(anet, -Inf, Inf, e=3)        # deact (-Inf, Inf) spell
b3 = identical(anet$mel[[3]]$atl$active,mat0)
deactivate.edges(anet, -Inf, Inf, e=4)        # deact (-Inf, Inf) spell
b4 = identical(anet$mel[[4]]$atl$active,mat0)
deactivate.edges(anet, -Inf, 10,  e=5)        # deact (-Inf, b) spell
b5 = identical(anet$mel[[5]]$atl$active,
              matrix(c(10, Inf), 1,2))
deactivate.edges(anet, 0, Inf, e=6)           # deact (a, Inf) spell
b6 = identical(anet$mel[[6]]$atl$active,
              matrix(c(-Inf, 0), 1,2))
deactivate.edges(anet, -10, Inf, e=7)         # deact (a, Inf) spell
b7 = identical(anet$mel[[7]]$atl$active,
              matrix(c(-Inf, -10), 1,2))
deactivate.edges(anet, 10, 20, e=8)           # deact (a, b) spell
b8 = identical(anet$mel[[8]]$atl$active,
              matrix(c(-Inf,20,10,Inf), 2,2))

b.tests = paste("b", seq(1,8), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  stop(paste("deactivate.edges is incorrectly activating edges in tests",
             bad.tests))
}


# for edges already having an active attribute
anet <- anet.copy

# Notes:
#    - a, b are upper, lower boundary points
#    - L is a finite number lower than a
#    - M is a finite number in between a and b
#    - H is a finite number higher than b
#    -Gij is a finite number that is in the gap b/t intervals i and j

# tests for deactivation at a point
for (i in 1:3)
  anet$mel[[i]]$atl$active <- matrix(c(Inf,Inf), 1,2)
for (i in 4:6)
  anet$mel[[i]]$atl$active <- matrix(c(-Inf,-Inf), 1,2)
activate.edges(anet, -Inf,  Inf, e=7:9)
activate.edges(anet, -Inf,   20, e=10:13)
activate.edges(anet,   20,   20, e=10:13)
activate.edges(anet,   10,  Inf, e=14:17)
activate.edges(anet,   10,   20, e=18:24)
activate.edges(anet,   20,   20, e=18:24)
activate.edges(anet,   10,   10, e=25:26)

deactivate.edges(anet, at=c(-Inf, Inf, 0, -Inf, Inf,  0, -Inf, Inf,  0, -Inf, 20,   0, 30, Inf, 10, 30,  0, -Inf, Inf, 10, 20, 15, 30, 0,  0,  10),
                      e = c(   1,   2, 3,    4,   5,  6,    7,   8,  9,   10, 11,  12, 13,  14, 15, 16, 17,   18,  19, 20, 21, 22, 23, 24, 25, 26))  
c1  = identical(anet$mel[[1]]$atl$active,    # deact -Inf in  (Inf, Inf)
      matrix(c(Inf, Inf),1,2))
c2  = identical(anet$mel[[2]]$atl$active,    # deact Inf  in  (Inf, Inf)
      matrix(c(Inf, Inf),1,2))
c3  = identical(anet$mel[[3]]$atl$active,    # deact M    in  (Inf, Inf)
      matrix(c(Inf, Inf),1,2))  
c4  = identical(anet$mel[[4]]$atl$active,    # deact -Inf in  (-Inf, -Inf)
      matrix(c(Inf, Inf),1,2))
c5  = identical(anet$mel[[5]]$atl$active,    # deact Inf  in  (-Inf, -Inf)
      matrix(c(Inf, Inf),1,2))
c6  = identical(anet$mel[[6]]$atl$active,    # deact M    in  (-Inf, -Inf)
      matrix(c(-Inf, -Inf),1,2))    
c7  = identical(anet$mel[[7]]$atl$active,    # deact -Inf in  (-Inf, Inf)
      matrix(c(Inf, Inf),1,2))    
c8  = identical(anet$mel[[8]]$atl$active,    # deact Inf  in  (-Inf, Inf)
      matrix(c(Inf, Inf),1,2))    
c9  = identical(anet$mel[[9]]$atl$active,    # deact M    in  (-Inf, Inf)
      matrix(c(-Inf, Inf),1,2))    
c10 = identical(anet$mel[[10]]$atl$active,   # deact -Inf in  (-Inf, b)
      matrix(c(Inf,Inf),1,2))    
c11 = identical(anet$mel[[11]]$atl$active,   # deact b    in  (-Inf, b)
      matrix(c(-Inf, 20),1,2))    
c12 = identical(anet$mel[[12]]$atl$active,   # deact M    in  (-Inf, b)
      matrix(c(-Inf,20,20,20),2,2))      
c13 = identical(anet$mel[[13]]$atl$active,   # deact H    in  (-Inf, b)
      matrix(c(-Inf,20,20,20),2,2))    
c14 = identical(anet$mel[[14]]$atl$active,   # deact Inf  in  (a, Inf)
      matrix(c(Inf, Inf),1,2))    
c15 = identical(anet$mel[[15]]$atl$active,   # deact a    in  (a, Inf)
      matrix(c(10, Inf),1,2))    
c16 = identical(anet$mel[[16]]$atl$active,   # deact M    in  (a, Inf)
      matrix(c(10, Inf),1,2))    
c17 = identical(anet$mel[[17]]$atl$active,   # deact L    in  (a, Inf)
      matrix(c(10, Inf),1,2))    
c18 = identical(anet$mel[[18]]$atl$active,   # deact -Inf in  (a,b)
      matrix(c(Inf,Inf),1,2))    
c19 = identical(anet$mel[[19]]$atl$active,   # deact Inf  in  (a,b)
      matrix(c(Inf,Inf),1,2))    
c20 = identical(anet$mel[[20]]$atl$active,   # deact a    in  (a,b)
      matrix(c(10,20,20,20),2,2))    
c21 = identical(anet$mel[[21]]$atl$active,   # deact b    in  (a,b)
      matrix(c(10, 20),1,2))    
c22 = identical(anet$mel[[22]]$atl$active,   # deact M    in  (a,b)
      matrix(c(10,20,20,20),2,2))    
c23 = identical(anet$mel[[23]]$atl$active,   # deact H    in  (a,b)
      matrix(c(10,20,20,20),2,2))    
c24 = identical(anet$mel[[24]]$atl$active,   # deact L    in  (a,b)
      matrix(c(10,20,20,20),2,2))    
c25 = identical(anet$mel[[25]]$atl$active,   # deact L    in  (a,a)
      matrix(c(10,10),1,2))    
c26 = identical(anet$mel[[26]]$atl$active,   # deact a    in  (a,a)
      matrix(c(Inf,Inf),1,2))    
  
c.tests = paste("c", seq(1,26), sep="")
c.results= sapply(c.tests, function(x){eval(parse(text=x))})
if(any(!c.results)){
  bad.tests = paste("c", which(!c.results), sep="", collapse=" ")
  stop(paste("deactivate.edges is incorrectly activating edges in tests",
             bad.tests))
}


anet <- anet.copy
for (i in 1:6)
  anet$mel[[i]]$atl$active <- matrix(c(Inf,Inf), 1,2)
deactivate.edges(anet, c(Inf, -Inf, -Inf, -Inf,  10, 10),
                       c(Inf, -Inf,  Inf,   10, Inf, 20),
                   e = c(  1,    2,    3,    4,   5,  6))
d1  = identical(anet$mel[[1]]$atl$active,   # deact (Inf, Inf)   in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d2  = identical(anet$mel[[2]]$atl$active,   # deact (-Inf, -Inf) in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d3  = identical(anet$mel[[3]]$atl$active,   # deact (-Inf, Inf)  in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d4  = identical(anet$mel[[4]]$atl$active,   # deact (-Inf, b)    in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d5  = identical(anet$mel[[5]]$atl$active,   # deact (a, Inf)     in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d6  = identical(anet$mel[[6]]$atl$active,   # deact (a, b)       in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))


for (i in 7:12)
  anet$mel[[i]]$atl$active <- matrix(c(-Inf, -Inf), 1,2)
deactivate.edges(anet, c(Inf, -Inf, -Inf, -Inf,  10, 10),
                       c(Inf, -Inf,  Inf,   10, Inf, 20),
                   e = c(  7,    8,    9,   10,  11,  12))
d7  = identical(anet$mel[[7]]$atl$active,   # deact (Inf, Inf)   in (-Inf, -Inf) spell
      matrix(c(Inf, Inf), 1,2))
d8  = identical(anet$mel[[8]]$atl$active,   # deact (-Inf, -Inf) in (-Inf, -Inf) spell
      matrix(c(Inf, Inf), 1,2))
d9  = identical(anet$mel[[9]]$atl$active,   # deact (-Inf, Inf)  in (-Inf, -Inf) spell
      matrix(c(Inf, Inf), 1,2))
d10 = identical(anet$mel[[10]]$atl$active,  # deact (-Inf, H)    in (-Inf, -Inf) spell
      matrix(c(-Inf, -Inf), 1,2))
d11 = identical(anet$mel[[11]]$atl$active,  # deact (H, Inf)     in (-Inf, -Inf) spell
      matrix(c(-Inf, -Inf), 1,2))
d12 = identical(anet$mel[[12]]$atl$active,  # deact (H, H)       in (-Inf, -Inf) spell
      matrix(c(-Inf, -Inf), 1,2))

for (i in 13:18)
  anet$mel[[i]]$atl$active <- matrix(c(-Inf, Inf), 1,2)
deactivate.edges(anet, c(Inf, -Inf, -Inf, -Inf,  10, 10),
                       c(Inf, -Inf,  Inf,   10, Inf, 20),
                   e = c( 13,   14,   15,   16,  17, 18))
d13 = identical(anet$mel[[13]]$atl$active,   # deact (Inf, Inf)   in (-Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d14 = identical(anet$mel[[14]]$atl$active,   # deact (-Inf, -Inf) in (-Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d15 = identical(anet$mel[[15]]$atl$active,   # deact (-Inf, Inf)  in (-Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d16 = identical(anet$mel[[16]]$atl$active,   # deact (-Inf, H)    in (-Inf, Inf) spell
      matrix(c(10, Inf), 1,2))
d17 = identical(anet$mel[[17]]$atl$active,   # deact (H, Inf)     in (-Inf, Inf) spell
      matrix(c(-Inf, 10), 1,2))
d18 = identical(anet$mel[[18]]$atl$active,   # deact (H1, H2)     in (-Inf, Inf) spell
      matrix(c(-Inf,20,10,Inf),2,2))


activate.edges(anet, -Inf, 10, e=seq(19,34))
deactivate.edges(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  20,  0, 10, 10,  0, -5,  0, 20),
                       c(Inf, -Inf,  Inf,    0,   10,   20, Inf, Inf, Inf, 10, 10, 20,  0,  0, 20, 30),
                   e = c( 19,   20,   21,   22,   23,   24,  25,  26,  27, 28, 29, 31, 31, 32, 33, 34))
d19 = identical(anet$mel[[19]]$atl$active,   # deact (Inf, Inf)   in (-Inf, b) spell
      matrix(c(Inf, Inf), 1,2))
d20 = identical(anet$mel[[20]]$atl$active,   # deact (-Inf, -Inf) in (-Inf, b) spell
      matrix(c(Inf, Inf), 1,2))
d21 = identical(anet$mel[[21]]$atl$active,   # deact (-Inf, Inf)  in (-Inf, b) spell
      matrix(c(Inf, Inf), 1,2))
d22 = identical(anet$mel[[22]]$atl$active,   # deact (-Inf, M)    in (-Inf, b) spell
      matrix(c(0, 10), 1,2))
d23 = identical(anet$mel[[23]]$atl$active,   # deact (-Inf, b)    in (-Inf, b) spell
      matrix(c(Inf, Inf), 1,2)) 
d24 = identical(anet$mel[[24]]$atl$active,   # deact (-Inf, H)    in (-Inf, b) spell
      matrix(c(Inf, Inf), 1,2))
d25 = identical(anet$mel[[25]]$atl$active,   # deact (M, Inf)     in (-Inf, b) spell
      matrix(c(-Inf, 0), 1,2))
d26 = identical(anet$mel[[26]]$atl$active,   # deact (b, Inf)     in (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d27 = identical(anet$mel[[27]]$atl$active,   # deact (H, Inf)     in (-Inf, b) spell
      matrix(c(-Inf,10),1,2))
d28 = identical(anet$mel[[28]]$atl$active,   # deact (M, b)       in (-Inf, b) spell
      matrix(c(-Inf, 0), 1,2))
d29 = identical(anet$mel[[29]]$atl$active,   # deact (b, b)       in (-Inf, b) spell
      matrix(c(-Inf, 10),1,2))
d30 = identical(anet$mel[[30]]$atl$active,   # deact (b, H)       in (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d31 = identical(anet$mel[[31]]$atl$active,   # deact (M, M)       in (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))  
d32 = identical(anet$mel[[32]]$atl$active,   # deact (M1, M2)     in (-Inf, b) spell
      matrix(c(-Inf,0,-5,10), 2,2))    
d33 = identical(anet$mel[[33]]$atl$active,   # deact (M, H)       in (-Inf, b) spell
      matrix(c(-Inf, 0), 1,2))  
d34 = identical(anet$mel[[34]]$atl$active,   # deact (H, H)       in (-Inf, b) spell
      matrix(c(-Inf,10),1,2))


activate.edges(anet, 10, Inf, e=seq(35,49))
deactivate.edges(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  20,  0, 10, 10,  0, -5, 20),
                       c(Inf, -Inf,  Inf,    0,   10,   20, Inf, Inf, Inf, 10, 10, 20,  5, 20, 30),
                   e = c( 35,   36,   37,   38,   39,   40,  41,  42,  43, 44, 45, 46, 47, 48, 49))
d35 = identical(anet$mel[[35]]$atl$active,   # deact (Inf, Inf)   in (a, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d36 = identical(anet$mel[[36]]$atl$active,   # deact (-Inf, -Inf) in (a, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d37 = identical(anet$mel[[37]]$atl$active,   # deact (-Inf, Inf)  in (a, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d38 = identical(anet$mel[[38]]$atl$active,   # deact (-Inf, L)    in (a, Inf) spell
      matrix(c(10,Inf),1,2))
d39 = identical(anet$mel[[39]]$atl$active,   # deact (-Inf, a)    in (a, Inf) spell
      matrix(c(10, Inf), 1,2)) 
d40 = identical(anet$mel[[40]]$atl$active,   # deact (-Inf, M)    in (a, Inf) spell
      matrix(c(20, Inf), 1,2))
d41 = identical(anet$mel[[41]]$atl$active,   # deact (L, Inf)     in (a, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d42 = identical(anet$mel[[42]]$atl$active,   # deact (a, Inf)     in (a, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d43 = identical(anet$mel[[43]]$atl$active,   # deact (M, Inf)     in (a, Inf) spell
      matrix(c(10, 20), 1,2))
d44 = identical(anet$mel[[44]]$atl$active,   # deact (L, a)       in (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d45 = identical(anet$mel[[45]]$atl$active,   # deact (a, a)       in (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d46 = identical(anet$mel[[46]]$atl$active,   # deact (a, M)       in (a, Inf) spell
      matrix(c(20, Inf), 1,2))
d47 = identical(anet$mel[[47]]$atl$active,   # deact (L1, L2)     in (a, Inf) spell
      matrix(c(10,Inf),1,2))
d48 = identical(anet$mel[[48]]$atl$active,   # deact (L, M)       in (a, Inf) spell
      matrix(c(20, Inf), 1,2))    
d49 = identical(anet$mel[[49]]$atl$active,   # deact (M1, M2)     in (a, Inf) spell
      matrix(c(10,30,20,Inf), 2,2))


activate.edges(anet, 10, 20, e=seq(50,76))
deactivate.edges(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  15,  20,  30,  0, 10, 10, 10, 10,  0, 20, 20, -5,  0,  0, 12, 15, 30),
                       c(Inf, -Inf,  Inf,    0,   10,   15,   20,   30, Inf, Inf, Inf, Inf, Inf, 10, 10, 15, 20, 30, 20, 20, 30,  0, 15, 30, 15, 30, 40),
                   e = c( 50,   51,   52,   53,   54,   55,   56,   57,  58,  59,  60,  61,  62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76))
d50 = identical(anet$mel[[50]]$atl$active,   # deact (Inf, Inf)   in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d51 = identical(anet$mel[[51]]$atl$active,   # deact (-Inf, -Inf) in (a, b) spell
      matrix(c(Inf, Inf), 1,2))  
d52 = identical(anet$mel[[52]]$atl$active,   # deact (-Inf, Inf)  in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d53 = identical(anet$mel[[53]]$atl$active,   # deact (-Inf, L)    in (a, b) spell
      matrix(c(10,20),1,2))
d54 = identical(anet$mel[[54]]$atl$active,   # deact (-Inf, a)    in (a, b) spell
      matrix(c(10, 20), 1,2)) 
d55 = identical(anet$mel[[55]]$atl$active,   # deact (-Inf, M)    in (a, b) spell
      matrix(c(15, 20), 1,2))
d56 = identical(anet$mel[[56]]$atl$active,   # deact (-Inf, b)    in (a, b) spell
      matrix(c(Inf, Inf), 1,2)) 
d57 = identical(anet$mel[[57]]$atl$active,   # deact (-Inf, H)    in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d58 = identical(anet$mel[[58]]$atl$active,   # deact (L, Inf)     in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d59 = identical(anet$mel[[59]]$atl$active,   # deact (a, Inf)     in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d60 = identical(anet$mel[[60]]$atl$active,   # deact (M, Inf)     in (a, b) spell
      matrix(c(10, 15), 1,2))
d61 = identical(anet$mel[[61]]$atl$active,   # deact (b, Inf)     in (a, b) spell
      matrix(c(10, 20), 1,2))
d62 = identical(anet$mel[[62]]$atl$active,   # deact (H, Inf)     in (a, b) spell
      matrix(c(10,20),1,2))
d63 = identical(anet$mel[[63]]$atl$active,   # deact (L, a)       in (a, b) spell
      matrix(c(10, 20), 1,2))
d64 = identical(anet$mel[[64]]$atl$active,   # deact (a, a)       in (a, b) spell
      matrix(c(10, 20), 1,2))
d65 = identical(anet$mel[[65]]$atl$active,   # deact (a, M)       in (a, b) spell
      matrix(c(15, 20), 1,2))
d66 = identical(anet$mel[[66]]$atl$active,   # deact (a, b)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d67 = identical(anet$mel[[67]]$atl$active,   # deact (a, H)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d68 = identical(anet$mel[[68]]$atl$active,   # deact (L, b)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d69 = identical(anet$mel[[69]]$atl$active,   # deact (b, b)       in (a, b) spell
      matrix(c(10, 20),1,2))
d70 = identical(anet$mel[[70]]$atl$active,   # deact (b, H)       in (a, b) spell
      matrix(c(10, 20), 1,2))
d71 = identical(anet$mel[[71]]$atl$active,   # deact (L1, L2)     in (a, b) spell
      matrix(c(10,20),1,2))
d72 = identical(anet$mel[[72]]$atl$active,   # deact (L, M)       in (a, b) spell
      matrix(c(15, 20), 1,2))    
d73 = identical(anet$mel[[73]]$atl$active,   # deact (L, H)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))    
d74 = identical(anet$mel[[74]]$atl$active,   # deact (M1, M2)     in (a, b) spell
      matrix(c(10,15,12,20), 2,2))
d75 = identical(anet$mel[[75]]$atl$active,   # deact (M, H)       in (a, b) spell
      matrix(c(10, 15), 1,2))    
d76 = identical(anet$mel[[76]]$atl$active,   # deact (H1, H2)     in (a, b) spell
      matrix(c(10,20),1,2))

activate.edges(anet, at=10, e=seq(77,91))
deactivate.edges(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  15,  0, 10, 10, -5,  0, 30),
                       c(Inf, -Inf,  Inf,    0,   10,   15, Inf, Inf, Inf, 10, 10, 15,  0, 15, 40),
                   e = c( 77,   78,   79,   80,   81,   82,  83,  84,  85, 86, 87, 88, 89, 90, 91))
d77 = identical(anet$mel[[77]]$atl$active,   # deact (Inf, Inf)   in (a, a) spell
      matrix(c(Inf, Inf), 1,2))
d78 = identical(anet$mel[[78]]$atl$active,   # deact (-Inf, -Inf) in (a, a) spell
      matrix(c(Inf, Inf), 1,2))  
d79 = identical(anet$mel[[79]]$atl$active,   # deact (-Inf, Inf)  in (a, a) spell
      matrix(c(Inf, Inf), 1,2))
d80 = identical(anet$mel[[80]]$atl$active,   # deact (-Inf, L)    in (a, a) spell
      matrix(c(10,10),1,2))
d81 = identical(anet$mel[[81]]$atl$active,   # deact (-Inf, a)    in (a, a) spell
      matrix(c(10,10), 1,2)) 
d82 = identical(anet$mel[[82]]$atl$active,   # deact (-Inf, H)    in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d83 = identical(anet$mel[[83]]$atl$active,   # deact (L, Inf)     in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d84 = identical(anet$mel[[84]]$atl$active,   # deact (a, Inf)     in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d85 = identical(anet$mel[[85]]$atl$active,   # deact (H, Inf)     in (a, b) spell
      matrix(c(10,10),1,2))
d86 = identical(anet$mel[[86]]$atl$active,   # deact (L, a)       in (a, b) spell
      matrix(c(10, 10), 1,2))
d87 = identical(anet$mel[[87]]$atl$active,   # deact (a, a)       in (a, b) spell
      matrix(c(Inf,Inf), 1,2))
d88 = identical(anet$mel[[88]]$atl$active,   # deact (a, H)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d89 = identical(anet$mel[[89]]$atl$active,   # deact (L1, L2)     in (a, b) spell
      matrix(c(10,10),1,2))
d90 = identical(anet$mel[[90]]$atl$active,   # deact (L, H)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))    
d91 = identical(anet$mel[[91]]$atl$active,   # deact (H1, H2)     in (a, b) spell
      matrix(c(10,10),1,2))

               
d.tests = paste("d", seq(1,91), sep="")
d.results= sapply(d.tests, function(x){eval(parse(text=x))})
if(any(!d.results)){
  bad.tests = paste("d", which(!d.results), sep="", collapse=" ")
  stop(paste("deactivate.edges is incorrectly activating edges in tests",
             bad.tests))
}

                 
anet<-anet.copy                 
activate.edges(anet, -Inf,  10, e=seq(1,28))
activate.edges(anet,   20,  30, e=seq(1,28))
activate.edges(anet,   30,  30, e=seq(1,28))
activate.edges(anet,   40,  50, e=seq(1,28))
activate.edges(anet,   60, Inf, e=seq(1,28))
deactivate.edges(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,  20,  25,  30,  35, 20, 20, 20, 20, 30, 30, 30, 15, 15, 15,  0,  0, 22,  0, 70, -Inf, 60),
                       c(Inf,  Inf,   20,   25,   30,   35,   55, Inf, Inf, Inf, Inf, 40, 45, 50, 55, 40, 45, 55, 35, 45, 55, 25, 55, 25,  5, 80,  10, Inf),
                   e = c(  1,    2,    3,    4,    5,    6,    7,   8,   9,  10,  11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,  27,  28))
emat0 = matrix(c(Inf,Inf),1,2)
emat1 = matrix(c(20,30,40,60,30,30,50,Inf),4,2)
emat2 = matrix(c(25,30,40,60,30,30,50,Inf),4,2)
emat3 = matrix(c(40,60,50,Inf),2,2)
emat4 = matrix(c(60,Inf),1,2)
emat5 = matrix(c(-Inf,10),1,2)
emat6 = matrix(c(-Inf,20,10,25),2,2)
emat7 = matrix(c(-Inf,20,10,30),2,2)
emat8 = matrix(c(-Inf,40,60,10,50,Inf),3,2)
emat9 = matrix(c(-Inf,45,60,10,50,Inf),3,2)
emat10 = matrix(c(-Inf,60,10,Inf),2,2)
emat11 = matrix(c(-Inf,20,40,60,10,30,50,Inf),4,2)
emat12 = matrix(c(-Inf,20,45,60,10,30,50,Inf),4,2)
emat13 = matrix(c(-Inf,20,60,10,30,Inf),3,2)
emat14 = matrix(c(-Inf,60,55,Inf),2,2)
emat15 = matrix(c(-Inf,20,30,10,30,30),3,2)
emat16 = matrix(c(30,40,60,30,50,Inf),3,2)
emat17 = matrix(c(-Inf,25,30,40,60,0,30,30,50,Inf),5,2)
emat18 = matrix(c(-Inf,60,0,Inf),2,2)
emat19 = matrix(c(-Inf,20,25,30,40,60,10,22,30,30,50,Inf),6,2)
emat20 = matrix(c(-Inf,5,20,30,40,60,0,10,30,30,50,Inf),6,2)
emat21 = matrix(c(-Inf,20,30,40,60,80,10,30,30,50,70,Inf),6,2)
emat22 = matrix(c(-Inf,20,30,40,10,30,30,50),4,2)

e1 = identical(anet$mel[[1]]$atl$active, emat0)    # deact (Inf, Inf)   in set of spells
e2 = identical(anet$mel[[2]]$atl$active, emat0)    # deact (-Inf, Inf)  in set of spells
e3 = identical(anet$mel[[3]]$atl$active, emat1)    # deact (-Inf, a2)   in set of spells
e4 = identical(anet$mel[[4]]$atl$active, emat2)    # deact (-Inf, M2)   in set of spells
e5 = identical(anet$mel[[5]]$atl$active, emat16)   # deact (-Inf, b2)   in set of spells
e6 = identical(anet$mel[[6]]$atl$active, emat3)    # deact (-Inf, G23)  in set of spells
e7 = identical(anet$mel[[7]]$atl$active, emat4)    # deact (-Inf, G34)  in set of spells
e8 = identical(anet$mel[[8]]$atl$active, emat5)    # deact (a2, Inf)    in set of spells
e9 = identical(anet$mel[[9]]$atl$active, emat6)    # deact (M2, Inf)    in set of spells
e10 = identical(anet$mel[[10]]$atl$active, emat7)  # deact (b2, Inf)    in set of spells
e11 = identical(anet$mel[[11]]$atl$active, emat15) # deact (G23, Inf)   in set of spells
e12 = identical(anet$mel[[12]]$atl$active, emat8)  # deact (a2, a3)     in set of spells
e13 = identical(anet$mel[[13]]$atl$active, emat9)  # deact (a2, M3)     in set of spells
e14 = identical(anet$mel[[14]]$atl$active, emat10) # deact (a2, b3)     in set of spells
e15= identical(anet$mel[[15]]$atl$active, emat10)  # deact (a2, G34)    in set of spells
e16= identical(anet$mel[[16]]$atl$active, emat11)  # deact (b2, a3)     in set of spells
e17= identical(anet$mel[[17]]$atl$active, emat12)  # deact (b2, M3)     in set of spells
e18= identical(anet$mel[[18]]$atl$active, emat13)  # deact (b2, G34)    in set of spells
e19 = identical(anet$mel[[19]]$atl$active, emat8)  # deact (G12, G23)   in set of spells
e20 = identical(anet$mel[[20]]$atl$active, emat9)  # deact (G12, M3)    in set of spells
e21 = identical(anet$mel[[21]]$atl$active, emat10) # deact (G12, G34)   in set of spells
e22 = identical(anet$mel[[22]]$atl$active, emat17) # deact (M1, M2)     in set of spells
e23 = identical(anet$mel[[23]]$atl$active, emat18) # deact (M1, G34)    in set of spells
e24 = identical(anet$mel[[24]]$atl$active, emat19) # deact (M2, M2)     in set of spells
e25 = identical(anet$mel[[25]]$atl$active, emat20) # deact (M1, M1)     in set of spells
e26 = identical(anet$mel[[26]]$atl$active, emat21) # deact (M6, M6)     in set of spells
e27 = identical(anet$mel[[27]]$atl$active, emat1)  # deact 1st spell    in set of spells
e28 = identical(anet$mel[[28]]$atl$active, emat22) # deact last spell   in set of spells
               
e.tests = paste("e", seq(1,28), sep="")
e.results= sapply(e.tests, function(x){eval(parse(text=x))})
if(any(!e.results)){
  bad.tests = paste("e", which(!e.results), sep="", collapse=" ")
  stop(paste("deactivate.edges is incorrectly activating edges in tests",
             bad.tests))
}


# other tests for 'deactivate.edges'
anet <- anet.copy
# default behavior
deactivate.edges(anet, e=1)                         # no inputs, no activity matrix
f1 = identical(anet$mel[[1]]$atl$active,emat0)
activate.edges(anet, 3, 10, e=2)                    # no inputs, activity matrix
deactivate.edges(anet, e=2)                     
f2 = identical(anet$mel[[2]]$atl$active,emat0)
anet$mel[[3]]$atl$active <- matrix(c(-Inf,-Inf), 1,2) # no inputs, matrix of null spell
deactivate.edges(anet, e=3) 
f3 = identical(anet$mel[[3]]$atl$active,emat0)
# using 'at' rather than 'onset' and 'terminus'
# ignored unless previously activated
activate.edges(anet, at=10, e=seq(4, 9))
deactivate.edges(anet, c(0, 10, 20), c(0, 10, 20), e=4:6)
deactivate.edges(anet, at=c(0, 10, 20), e=7:9)
f4 = identical(anet$mel[[4]]$atl$active, anet$mel[[7]]$atl$active)
f5 = identical(anet$mel[[5]]$atl$active, anet$mel[[8]]$atl$active)
f6 = identical(anet$mel[[6]]$atl$active, anet$mel[[9]]$atl$active)
# or 'length' rather than 'terminus'
activate.edges(anet, -Inf,  10, e=seq(25, 43))
activate.edges(anet,   20,  30, e=seq(25, 43))
activate.edges(anet,   40, Inf, e=seq(25, 43))
deactivate.edges(anet, c(15, 15, 15, 15, 15, 20, 20, 20, 20),
                       c(18, 20, 25, 30, 35, 25, 30, 35, 40),
                     e=c(25, 26, 27, 28, 29, 30, 31, 32, 33))
deactivate.edges(anet, c(15, 15, 15, 15, 15, 20, 20, 20, 20),
               length =c( 3,  5, 10, 15, 20,  5, 10, 15, 20),
                    e =c(35, 36, 37, 38, 39, 40, 41, 42, 43))
f7  = identical(anet$mel[[25]]$atl$active, anet$mel[[35]]$atl$active)
f8  = identical(anet$mel[[26]]$atl$active, anet$mel[[36]]$atl$active)
f9  = identical(anet$mel[[27]]$atl$active, anet$mel[[37]]$atl$active)
f10 = identical(anet$mel[[28]]$atl$active, anet$mel[[38]]$atl$active)
f11 = identical(anet$mel[[29]]$atl$active, anet$mel[[39]]$atl$active)
f12 = identical(anet$mel[[30]]$atl$active, anet$mel[[40]]$atl$active)
f13 = identical(anet$mel[[31]]$atl$active, anet$mel[[41]]$atl$active)
f14 = identical(anet$mel[[32]]$atl$active, anet$mel[[42]]$atl$active)
f15 = identical(anet$mel[[33]]$atl$active, anet$mel[[43]]$atl$active)


f.tests = paste("f", seq(1,15), sep="")
f.results= sapply(f.tests, function(x){eval(parse(text=x))})
if(any(!f.results)){
  bad.tests = paste("f", which(!f.results), sep="", collapse=" ")
  stop(paste("deactivate.edges is incorrectly activating edges in tests",
             bad.tests))
}

expect_equal(network.edgecount.active(deactivate.edges(network.initialize(0))),0)

cat("ok\n")

#------------------- DEACTIVATE.VERTICES TESTS ------------------------
# Notes:
#  --except for the error-checking code, the rest of the tests
#    are meant to be exhaustive
#-----------------------------------------------------------------

cat("testing deactivate.vertices ... ")
anet <- anet.copy

# proper behavior for bad inputs?
#  -- bad network input
a1 = try(deactivate.vertices(3, v=1), T)                    # properly handeled
#  -- bad at times
a2 = try(deactivate.vertices(anet, at="a", v=2), T)         # properly handeled
#  -- bad onset
a3 = try(deactivate.vertices(anet, "a", v=3), T)            # properly handeled
a4 = try(deactivate.vertices(anet, NULL, 3, v=4), T)        # properly handeled
# -- bad terminus
a5 = try(deactivate.vertices(anet, 3, "b", v=5), T)         # properly handeled
# -- bad length
a6 = try(deactivate.vertices(anet, 3, length=-9, v=6), T)   # properly handeled
a7 = try(deactivate.vertices(anet, 3, length="r", v=7), T)  # properly handeled
# -- bad vertices
a8 = try(deactivate.vertices(anet, 3, 10, v=174), T)        # properly handeled
a9 = try(deactivate.vertices(anet, 3, 10, v=-2), T)         # properly handeled
a10 = try(deactivate.vertices(anet, 3, 10, v=NULL), T)      # properly handeled
a11 = try(deactivate.vertices(anet, 3, 10, v="hello"), T)   # properly handeled
# -- bad onset & terminus combo
a12 = try(deactivate.vertices(anet, 10, 3, v=8), T)         # properly handeled
# -- not fully specified intervals
a13 = try(deactivate.vertices(anet, 10, v=9), T)            # properly handeled
a14 = try(deactivate.vertices(anet, terminus=10, v=10), T)  # properly handeled
a15 = try(deactivate.vertices(anet, length=10, v=11), T)    # properly handeled


# good input
anet <- anet.copy
mat0 <- matrix(c(Inf, Inf), 1,2)
# for vertices lacking an active attribute
deactivate.vertices(anet, Inf, Inf, v=1)         # deact (Inf, Inf) spell
b1 = identical(anet$val[[1]]$active,mat0)
deactivate.vertices(anet, -Inf, -Inf, v=2)       # deact (-Inf, -Inf) spell
b2 = identical(anet$val[[2]]$active,mat0)
deactivate.vertices(anet, -Inf, Inf, v=3)        # deact (-Inf, Inf) spell
b3 = identical(anet$val[[3]]$active,mat0)
deactivate.vertices(anet, -Inf, Inf, v=4)        # deact (-Inf, Inf) spell
b4 = identical(anet$val[[4]]$active,mat0)
deactivate.vertices(anet, -Inf, 10,  v=5)        # deact (-Inf, b) spell
b5 = identical(anet$val[[5]]$active,
              matrix(c(10, Inf), 1,2))
deactivate.vertices(anet, 0, Inf, v=6)           # deact (a, Inf) spell
b6 = identical(anet$val[[6]]$active,
              matrix(c(-Inf, 0), 1,2))
deactivate.vertices(anet, -10, Inf, v=7)         # deact (a, Inf) spell
b7 = identical(anet$val[[7]]$active,
              matrix(c(-Inf, -10), 1,2))
deactivate.vertices(anet, 10, 20, v=8)           # deact (a, b) spell
b8 = identical(anet$val[[8]]$active,
              matrix(c(-Inf,20,10,Inf), 2,2))

b.tests = paste("b", seq(1,8), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  stop(paste("deactivate.vertices is incorrectly activating vertices in tests",
             bad.tests))
}


# for vertices already having an active attribute
anet <- anet.copy

# Notes:
#    - a, b are upper, lower boundary points
#    - L is a finite number lower than a
#    - M is a finite number in between a and b
#    - H is a finite number higher than b
#    -Gij is a finite number that is in the gap b/t intervals i and j

# tests for deactivation at a point
for (i in 1:3)
  anet$val[[i]]$active <- matrix(c(Inf,Inf), 1,2)
for (i in 4:6)
  anet$val[[i]]$active <- matrix(c(-Inf,-Inf), 1,2)
activate.vertices(anet, -Inf,  Inf, v=7:9)
activate.vertices(anet, -Inf,   20, v=10:13)
activate.vertices(anet,   20,   20, v=10:13)
activate.vertices(anet,   10,  Inf, v=14:17)
activate.vertices(anet,   10,   20, v=18:24)
activate.vertices(anet,   20,   20, v=18:24)
activate.vertices(anet,   10,   10, v=25:26)

deactivate.vertices(anet, at=c(-Inf, Inf, 0, -Inf, Inf,  0, -Inf, Inf,  0, -Inf, 20,   0, 30, Inf, 10, 30,  0, -Inf, Inf, 10, 20, 15, 30, 0,  0,  10),
                      v = c(   1,   2, 3,    4,   5,  6,    7,   8,  9,   10, 11,  12, 13,  14, 15, 16, 17,   18,  19, 20, 21, 22, 23, 24, 25, 26))  
c1  = identical(anet$val[[1]]$active,    # deact -Inf in  (Inf, Inf)
      matrix(c(Inf, Inf),1,2))
c2  = identical(anet$val[[2]]$active,    # deact Inf  in  (Inf, Inf)
      matrix(c(Inf, Inf),1,2))
c3  = identical(anet$val[[3]]$active,    # deact M    in  (Inf, Inf)
      matrix(c(Inf, Inf),1,2))  
c4  = identical(anet$val[[4]]$active,    # deact -Inf in  (-Inf, -Inf)
      matrix(c(Inf, Inf),1,2))
c5  = identical(anet$val[[5]]$active,    # deact Inf  in  (-Inf, -Inf)
      matrix(c(Inf, Inf),1,2))
c6  = identical(anet$val[[6]]$active,    # deact M    in  (-Inf, -Inf)
      matrix(c(-Inf, -Inf),1,2))    
c7  = identical(anet$val[[7]]$active,    # deact -Inf in  (-Inf, Inf)
      matrix(c(Inf, Inf),1,2))    
c8  = identical(anet$val[[8]]$active,    # deact Inf  in  (-Inf, Inf)
      matrix(c(Inf, Inf),1,2))    
c9  = identical(anet$val[[9]]$active,    # deact M    in  (-Inf, Inf)
      matrix(c(-Inf, Inf),1,2))    
c10 = identical(anet$val[[10]]$active,   # deact -Inf in  (-Inf, b)
      matrix(c(Inf, Inf),1,2))    
c11 = identical(anet$val[[11]]$active,   # deact b    in  (-Inf, b)
      matrix(c(-Inf, 20),1,2))    
c12 = identical(anet$val[[12]]$active,   # deact M    in  (-Inf, b)
      matrix(c(-Inf,20,20,20),2,2))      
c13 = identical(anet$val[[13]]$active,   # deact H    in  (-Inf, b)
      matrix(c(-Inf,20,20,20),2,2))    
c14 = identical(anet$val[[14]]$active,   # deact Inf  in  (a, Inf)
      matrix(c(Inf, Inf),1,2))    
c15 = identical(anet$val[[15]]$active,   # deact a    in  (a, Inf)
      matrix(c(10, Inf),1,2))    
c16 = identical(anet$val[[16]]$active,   # deact M    in  (a, Inf)
      matrix(c(10, Inf),1,2))    
c17 = identical(anet$val[[17]]$active,   # deact L    in  (a, Inf)
      matrix(c(10, Inf),1,2))    
c18 = identical(anet$val[[18]]$active,   # deact -Inf in  (a,b)
      matrix(c(Inf, Inf),1,2))    
c19 = identical(anet$val[[19]]$active,   # deact Inf  in  (a,b)
      matrix(c(Inf, Inf),1,2))    
c20 = identical(anet$val[[20]]$active,   # deact a    in  (a,b)
      matrix(c(10,20,20,20),2,2))    
c21 = identical(anet$val[[21]]$active,   # deact b    in  (a,b)
      matrix(c(10, 20),1,2))    
c22 = identical(anet$val[[22]]$active,   # deact M    in  (a,b)
      matrix(c(10,20,20,20),2,2))    
c23 = identical(anet$val[[23]]$active,   # deact H    in  (a,b)
      matrix(c(10,20,20,20),2,2))    
c24 = identical(anet$val[[24]]$active,   # deact L    in  (a,b)
      matrix(c(10,20,20,20),2,2))    
c25 = identical(anet$val[[25]]$active,   # deact L    in  (a,a)
      matrix(c(10,10),1,2))    
c26 = identical(anet$val[[26]]$active,   # deact a    in  (a,a)
      matrix(c(Inf,Inf),1,2))    
  
c.tests = paste("c", seq(1,26), sep="")
c.results= sapply(c.tests, function(x){eval(parse(text=x))})
if(any(!c.results)){
  bad.tests = paste("c", which(!c.results), sep="", collapse=" ")
  stop(paste("deactivate.vertices is incorrectly activating vertices in tests",
             bad.tests))
}


anet <- anet.copy
for (i in 1:6)
  anet$val[[i]]$active <- matrix(c(Inf,Inf), 1,2)
deactivate.vertices(anet, c(Inf, -Inf, -Inf, -Inf,  10, 10),
                       c(Inf, -Inf,  Inf,   10, Inf, 20),
                   v = c(  1,    2,    3,    4,   5,  6))
d1  = identical(anet$val[[1]]$active,   # deact (Inf, Inf)   in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d2  = identical(anet$val[[2]]$active,   # deact (-Inf, -Inf) in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d3  = identical(anet$val[[3]]$active,   # deact (-Inf, Inf)  in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d4  = identical(anet$val[[4]]$active,   # deact (-Inf, b)    in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d5  = identical(anet$val[[5]]$active,   # deact (a, Inf)     in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d6  = identical(anet$val[[6]]$active,   # deact (a, b)       in (Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))


for (i in 7:12)
  anet$val[[i]]$active <- matrix(c(-Inf, -Inf), 1,2)
deactivate.vertices(anet, c(Inf, -Inf, -Inf, -Inf,  10, 10),
                       c(Inf, -Inf,  Inf,   10, Inf, 20),
                   v = c(  7,    8,    9,   10,  11,  12))
d7  = identical(anet$val[[7]]$active,   # deact (Inf, Inf)   in (-Inf, -Inf) spell
      matrix(c(Inf, Inf), 1,2))
d8  = identical(anet$val[[8]]$active,   # deact (-Inf, -Inf) in (-Inf, -Inf) spell
      matrix(c(Inf, Inf), 1,2))
d9  = identical(anet$val[[9]]$active,   # deact (-Inf, Inf)  in (-Inf, -Inf) spell
      matrix(c(Inf, Inf), 1,2))
d10 = identical(anet$val[[10]]$active,  # deact (-Inf, H)    in (-Inf, -Inf) spell
      matrix(c(-Inf, -Inf), 1,2))
d11 = identical(anet$val[[11]]$active,  # deact (H, Inf)     in (-Inf, -Inf) spell
      matrix(c(-Inf, -Inf), 1,2))
d12 = identical(anet$val[[12]]$active,  # deact (H, H)       in (-Inf, -Inf) spell
      matrix(c(-Inf, -Inf), 1,2))

for (i in 13:18)
  anet$val[[i]]$active <- matrix(c(-Inf, Inf), 1,2)
deactivate.vertices(anet, c(Inf, -Inf, -Inf, -Inf,  10, 10),
                       c(Inf, -Inf,  Inf,   10, Inf, 20),
                   v = c( 13,   14,   15,   16,  17, 18))
d13 = identical(anet$val[[13]]$active,   # deact (Inf, Inf)   in (-Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d14 = identical(anet$val[[14]]$active,   # deact (-Inf, -Inf) in (-Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d15 = identical(anet$val[[15]]$active,   # deact (-Inf, Inf)  in (-Inf, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d16 = identical(anet$val[[16]]$active,   # deact (-Inf, H)    in (-Inf, Inf) spell
      matrix(c(10, Inf), 1,2))
d17 = identical(anet$val[[17]]$active,   # deact (H, Inf)     in (-Inf, Inf) spell
      matrix(c(-Inf, 10), 1,2))
d18 = identical(anet$val[[18]]$active,   # deact (H1, H2)     in (-Inf, Inf) spell
      matrix(c(-Inf,20,10,Inf),2,2))


activate.vertices(anet, -Inf, 10, v=seq(19,34))
deactivate.vertices(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  20,  0, 10, 10,  0, -5,  0, 20),
                       c(Inf, -Inf,  Inf,    0,   10,   20, Inf, Inf, Inf, 10, 10, 20,  0,  0, 20, 30),
                   v = c( 19,   20,   21,   22,   23,   24,  25,  26,  27, 28, 29, 31, 31, 32, 33, 34))
d19 = identical(anet$val[[19]]$active,   # deact (Inf, Inf)   in (-Inf, b) spell
      matrix(c(Inf, Inf), 1,2))
d20 = identical(anet$val[[20]]$active,   # deact (-Inf, -Inf) in (-Inf, b) spell
      matrix(c(Inf, Inf), 1,2))
d21 = identical(anet$val[[21]]$active,   # deact (-Inf, Inf)  in (-Inf, b) spell
      matrix(c(Inf, Inf), 1,2))
d22 = identical(anet$val[[22]]$active,   # deact (-Inf, M)    in (-Inf, b) spell
      matrix(c(0, 10), 1,2))
d23 = identical(anet$val[[23]]$active,   # deact (-Inf, b)    in (-Inf, b) spell
      matrix(c(Inf, Inf), 1,2)) 
d24 = identical(anet$val[[24]]$active,   # deact (-Inf, H)    in (-Inf, b) spell
      matrix(c(Inf, Inf), 1,2))
d25 = identical(anet$val[[25]]$active,   # deact (M, Inf)     in (-Inf, b) spell
      matrix(c(-Inf, 0), 1,2))
d26 = identical(anet$val[[26]]$active,   # deact (b, Inf)     in (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d27 = identical(anet$val[[27]]$active,   # deact (H, Inf)     in (-Inf, b) spell
      matrix(c(-Inf,10),1,2))
d28 = identical(anet$val[[28]]$active,   # deact (M, b)       in (-Inf, b) spell
      matrix(c(-Inf, 0), 1,2))
d29 = identical(anet$val[[29]]$active,   # deact (b, b)       in (-Inf, b) spell
      matrix(c(-Inf, 10),1,2))
d30 = identical(anet$val[[30]]$active,   # deact (b, H)       in (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))
d31 = identical(anet$val[[31]]$active,   # deact (M, M)       in (-Inf, b) spell
      matrix(c(-Inf, 10), 1,2))  
d32 = identical(anet$val[[32]]$active,   # deact (M1, M2)     in (-Inf, b) spell
      matrix(c(-Inf,0,-5,10), 2,2))    
d33 = identical(anet$val[[33]]$active,   # deact (M, H)       in (-Inf, b) spell
      matrix(c(-Inf, 0), 1,2))  
d34 = identical(anet$val[[34]]$active,   # deact (H, H)       in (-Inf, b) spell
      matrix(c(-Inf,10),1,2))


activate.vertices(anet, 10, Inf, v=seq(35,49))
deactivate.vertices(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  20,  0, 10, 10,  0, -5, 20),
                       c(Inf, -Inf,  Inf,    0,   10,   20, Inf, Inf, Inf, 10, 10, 20,  5, 20, 30),
                   v = c( 35,   36,   37,   38,   39,   40,  41,  42,  43, 44, 45, 46, 47, 48, 49))
d35 = identical(anet$val[[35]]$active,   # deact (Inf, Inf)   in (a, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d36 = identical(anet$val[[36]]$active,   # deact (-Inf, -Inf) in (a, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d37 = identical(anet$val[[37]]$active,   # deact (-Inf, Inf)  in (a, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d38 = identical(anet$val[[38]]$active,   # deact (-Inf, L)    in (a, Inf) spell
      matrix(c(10,Inf),1,2))
d39 = identical(anet$val[[39]]$active,   # deact (-Inf, a)    in (a, Inf) spell
      matrix(c(10, Inf), 1,2)) 
d40 = identical(anet$val[[40]]$active,   # deact (-Inf, M)    in (a, Inf) spell
      matrix(c(20, Inf), 1,2))
d41 = identical(anet$val[[41]]$active,   # deact (L, Inf)     in (a, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d42 = identical(anet$val[[42]]$active,   # deact (a, Inf)     in (a, Inf) spell
      matrix(c(Inf, Inf), 1,2))
d43 = identical(anet$val[[43]]$active,   # deact (M, Inf)     in (a, Inf) spell
      matrix(c(10, 20), 1,2))
d44 = identical(anet$val[[44]]$active,   # deact (L, a)       in (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d45 = identical(anet$val[[45]]$active,   # deact (a, a)       in (a, Inf) spell
      matrix(c(10, Inf), 1,2))
d46 = identical(anet$val[[46]]$active,   # deact (a, M)       in (a, Inf) spell
      matrix(c(20, Inf), 1,2))
d47 = identical(anet$val[[47]]$active,   # deact (L1, L2)     in (a, Inf) spell
      matrix(c(10,Inf),1,2))
d48 = identical(anet$val[[48]]$active,   # deact (L, M)       in (a, Inf) spell
      matrix(c(20, Inf), 1,2))    
d49 = identical(anet$val[[49]]$active,   # deact (M1, M2)     in (a, Inf) spell
      matrix(c(10,30,20,Inf), 2,2))


activate.vertices(anet, 10, 20, v=seq(50,76))
deactivate.vertices(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  15,  20,  30,  0, 10, 10, 10, 10,  0, 20, 20, -5,  0,  0, 12, 15, 30),
                       c(Inf, -Inf,  Inf,    0,   10,   15,   20,   30, Inf, Inf, Inf, Inf, Inf, 10, 10, 15, 20, 30, 20, 20, 30,  0, 15, 30, 15, 30, 40),
                   v = c( 50,   51,   52,   53,   54,   55,   56,   57,  58,  59,  60,  61,  62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76))
d50 = identical(anet$val[[50]]$active,   # deact (Inf, Inf)   in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d51 = identical(anet$val[[51]]$active,   # deact (-Inf, -Inf) in (a, b) spell
      matrix(c(Inf, Inf), 1,2))  
d52 = identical(anet$val[[52]]$active,   # deact (-Inf, Inf)  in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d53 = identical(anet$val[[53]]$active,   # deact (-Inf, L)    in (a, b) spell
      matrix(c(10,20),1,2))
d54 = identical(anet$val[[54]]$active,   # deact (-Inf, a)    in (a, b) spell
      matrix(c(10, 20), 1,2)) 
d55 = identical(anet$val[[55]]$active,   # deact (-Inf, M)    in (a, b) spell
      matrix(c(15, 20), 1,2))
d56 = identical(anet$val[[56]]$active,   # deact (-Inf, b)    in (a, b) spell
      matrix(c(Inf, Inf), 1,2)) 
d57 = identical(anet$val[[57]]$active,   # deact (-Inf, H)    in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d58 = identical(anet$val[[58]]$active,   # deact (L, Inf)     in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d59 = identical(anet$val[[59]]$active,   # deact (a, Inf)     in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d60 = identical(anet$val[[60]]$active,   # deact (M, Inf)     in (a, b) spell
      matrix(c(10, 15), 1,2))
d61 = identical(anet$val[[61]]$active,   # deact (b, Inf)     in (a, b) spell
      matrix(c(10, 20), 1,2))
d62 = identical(anet$val[[62]]$active,   # deact (H, Inf)     in (a, b) spell
      matrix(c(10,20),1,2))
d63 = identical(anet$val[[63]]$active,   # deact (L, a)       in (a, b) spell
      matrix(c(10, 20), 1,2))
d64 = identical(anet$val[[64]]$active,   # deact (a, a)       in (a, b) spell
      matrix(c(10, 20), 1,2))
d65 = identical(anet$val[[65]]$active,   # deact (a, M)       in (a, b) spell
      matrix(c(15, 20), 1,2))
d66 = identical(anet$val[[66]]$active,   # deact (a, b)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d67 = identical(anet$val[[67]]$active,   # deact (a, H)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d68 = identical(anet$val[[68]]$active,   # deact (L, b)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d69 = identical(anet$val[[69]]$active,   # deact (b, b)       in (a, b) spell
      matrix(c(10, 20),1,2))
d70 = identical(anet$val[[70]]$active,   # deact (b, H)       in (a, b) spell
      matrix(c(10, 20), 1,2))
d71 = identical(anet$val[[71]]$active,   # deact (L1, L2)     in (a, b) spell
      matrix(c(10,20),1,2))
d72 = identical(anet$val[[72]]$active,   # deact (L, M)       in (a, b) spell
      matrix(c(15, 20), 1,2))    
d73 = identical(anet$val[[73]]$active,   # deact (L, H)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))    
d74 = identical(anet$val[[74]]$active,   # deact (M1, M2)     in (a, b) spell
      matrix(c(10,15,12,20), 2,2))
d75 = identical(anet$val[[75]]$active,   # deact (M, H)       in (a, b) spell
      matrix(c(10, 15), 1,2))    
d76 = identical(anet$val[[76]]$active,   # deact (H1, H2)     in (a, b) spell
      matrix(c(10,20),1,2))

activate.vertices(anet, at=10, v=seq(77,91))
deactivate.vertices(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf,   0,  10,  15,  0, 10, 10, -5,  0, 30),
                       c(Inf, -Inf,  Inf,    0,   10,   15, Inf, Inf, Inf, 10, 10, 15,  0, 15, 40),
                   v = c( 77,   78,   79,   80,   81,   82,  83,  84,  85, 86, 87, 88, 89, 90, 91))
d77 = identical(anet$val[[77]]$active,   # deact (Inf, Inf)   in (a, a) spell
      matrix(c(Inf, Inf), 1,2))
d78 = identical(anet$val[[78]]$active,   # deact (-Inf, -Inf) in (a, a) spell
      matrix(c(Inf, Inf), 1,2))  
d79 = identical(anet$val[[79]]$active,   # deact (-Inf, Inf)  in (a, a) spell
      matrix(c(Inf, Inf), 1,2))
d80 = identical(anet$val[[80]]$active,   # deact (-Inf, L)    in (a, a) spell
      matrix(c(10,10),1,2))
d81 = identical(anet$val[[81]]$active,   # deact (-Inf, a)    in (a, a) spell
      matrix(c(10,10), 1,2)) 
d82 = identical(anet$val[[82]]$active,   # deact (-Inf, H)    in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d83 = identical(anet$val[[83]]$active,   # deact (L, Inf)     in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d84 = identical(anet$val[[84]]$active,   # deact (a, Inf)     in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d85 = identical(anet$val[[85]]$active,   # deact (H, Inf)     in (a, b) spell
      matrix(c(10,10),1,2))
d86 = identical(anet$val[[86]]$active,   # deact (L, a)       in (a, b) spell
      matrix(c(10, 10), 1,2))
d87 = identical(anet$val[[87]]$active,   # deact (a, a)       in (a, b) spell
      matrix(c(Inf,Inf), 1,2))
d88 = identical(anet$val[[88]]$active,   # deact (a, H)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))
d89 = identical(anet$val[[89]]$active,   # deact (L1, L2)     in (a, b) spell
      matrix(c(10,10),1,2))
d90 = identical(anet$val[[90]]$active,   # deact (L, H)       in (a, b) spell
      matrix(c(Inf, Inf), 1,2))    
d91 = identical(anet$val[[91]]$active,   # deact (H1, H2)     in (a, b) spell
      matrix(c(10,10),1,2))

               
d.tests = paste("d", seq(1,91), sep="")
d.results= sapply(d.tests, function(x){eval(parse(text=x))})
if(any(!d.results)){
  bad.tests = paste("d", which(!d.results), sep="", collapse=" ")
  stop(paste("deactivate.vertices is incorrectly activating vertices in tests",
             bad.tests))
}

                 
anet<-anet.copy                 
activate.vertices(anet, -Inf,  10, v=seq(1,23))
activate.vertices(anet,   20,  30, v=seq(1,23))
activate.vertices(anet,   30,  30, v=seq(1,23))
activate.vertices(anet,   40,  50, v=seq(1,23))
activate.vertices(anet,   60, Inf, v=seq(1,23))
deactivate.vertices(anet, c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,  20,  25,  30,  35, 20, 20, 20, 20, 30, 30, 30, 15, 15, 15,  0, 0),
                       c(Inf,  Inf,   20,   25,   30,   35,   55, Inf, Inf, Inf, Inf, 40, 45, 50, 55, 40, 45, 55, 35, 45, 55, 25, 55),
                   v = c(  1,    2,    3,    4,    5,    6,    7,   8,   9,  10,  11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23))
emat0 = matrix(c(Inf,Inf),1,2)
emat1 = matrix(c(20,30,40,60,30,30,50,Inf),4,2)
emat2 = matrix(c(25,30,40,60,30,30,50,Inf),4,2)
emat3 = matrix(c(40,60,50,Inf),2,2)
emat4 = matrix(c(60,Inf),1,2)
emat5 = matrix(c(-Inf,10),1,2)
emat6 = matrix(c(-Inf,20,10,25),2,2)
emat7 = matrix(c(-Inf,20,10,30),2,2)
emat8 = matrix(c(-Inf,40,60,10,50,Inf),3,2)
emat9 = matrix(c(-Inf,45,60,10,50,Inf),3,2)
emat10 = matrix(c(-Inf,60,10,Inf),2,2)
emat11 = matrix(c(-Inf,20,40,60,10,30,50,Inf),4,2)
emat12 = matrix(c(-Inf,20,45,60,10,30,50,Inf),4,2)
emat13 = matrix(c(-Inf,20,60,10,30,Inf),3,2)
emat14 = matrix(c(-Inf,60,55,Inf),2,2)
emat15 = matrix(c(-Inf,20,30,10,30,30),3,2)
emat16 = matrix(c(30,40,60,30,50,Inf),3,2)
emat17 = matrix(c(-Inf,25,30,40,60,0,30,30,50,Inf),5,2)
emat18 = matrix(c(-Inf,60,0,Inf),2,2)

e1 = identical(anet$val[[1]]$active, emat0)  # deact (Inf, Inf)   in set of spells
e2 = identical(anet$val[[2]]$active, emat0)  # deact (-Inf, Inf)  in set of spells
e3 = identical(anet$val[[3]]$active, emat1)  # deact (-Inf, a2)   in set of spells
e4 = identical(anet$val[[4]]$active, emat2)  # deact (-Inf, M2)   in set of spells
e5 = identical(anet$val[[5]]$active, emat16) # deact (-Inf, b2)   in set of spells
e6 = identical(anet$val[[6]]$active, emat3)  # deact (-Inf, G23)  in set of spells
e7 = identical(anet$val[[7]]$active, emat4)  # deact (-Inf, G34)  in set of spells
e8 = identical(anet$val[[8]]$active, emat5)  # deact (a2, Inf)    in set of spells
e9 = identical(anet$val[[9]]$active, emat6)  # deact (M2, Inf)    in set of spells
e10 = identical(anet$val[[10]]$active, emat7)  # deact (b2, Inf)    in set of spells
e11 = identical(anet$val[[11]]$active, emat15) # deact (G23, Inf)   in set of spells
e12 = identical(anet$val[[12]]$active, emat8)  # deact (a2, a3)     in set of spells
e13 = identical(anet$val[[13]]$active, emat9)  # deact (a2, M3)     in set of spells
e14 = identical(anet$val[[14]]$active, emat10) # deact (a2, b3)     in set of spells
e15= identical(anet$val[[15]]$active, emat10) # deact (a2, G34)    in set of spells
e16= identical(anet$val[[16]]$active, emat11) # deact (b2, a3)     in set of spells
e17= identical(anet$val[[17]]$active, emat12) # deact (b2, M3)     in set of spells
e18= identical(anet$val[[18]]$active, emat13) # deact (b2, G34)    in set of spells
e19 = identical(anet$val[[19]]$active, emat8)  # deact (G12, G23)   in set of spells
e20 = identical(anet$val[[20]]$active, emat9)  # deact (G12, M3)    in set of spells
e21 = identical(anet$val[[21]]$active, emat10) # deact (G12, G34)   in set of spells
e22 = identical(anet$val[[22]]$active, emat17) # deact (M1, M2)     in set of spells
e23 = identical(anet$val[[23]]$active, emat18) # deact (M1, G34)    in set of spells
               
e.tests = paste("e", seq(1,23), sep="")
e.results= sapply(e.tests, function(x){eval(parse(text=x))})
if(any(!e.results)){
  bad.tests = paste("e", which(!e.results), sep="", collapse=" ")
  stop(paste("deactivate.vertices is incorrectly activating vertices in tests",
             bad.tests))
}


# other tests for 'deactivate.vertices'
anet <- anet.copy
# default behavior
deactivate.vertices(anet, v=1)                         # no inputs, no activity matrix
f1 = identical(anet$val[[1]]$active,emat0)
activate.vertices(anet, 3, 10, v=2)                    # no inputs, activity matrix
deactivate.vertices(anet, v=2)                     
f2 = identical(anet$val[[2]]$active,emat0)
anet$val[[3]]$active <- matrix(c(-Inf,-Inf), 1,2) # no inputs, matrix of null spell
deactivate.vertices(anet, v=3) 
f3 = identical(anet$val[[3]]$active,emat0)
# using 'at' rather than 'onset' and 'terminus'
# ignored unless previously activated
activate.vertices(anet, at=10, v=seq(4, 9))
deactivate.vertices(anet, c(0, 10, 20), c(0, 10, 20), v=4:6)
deactivate.vertices(anet, at=c(0, 10, 20), v=7:9)
f4 = identical(anet$val[[4]]$active, anet$val[[7]]$active)
f5 = identical(anet$val[[5]]$active, anet$val[[8]]$active)
f6 = identical(anet$val[[6]]$active, anet$val[[9]]$active)
# or 'length' rather than 'terminus'
activate.vertices(anet, -Inf,  10, v=seq(25, 43))
activate.vertices(anet,   20,  30, v=seq(25, 43))
activate.vertices(anet,   40, Inf, v=seq(25, 43))
deactivate.vertices(anet, c(15, 15, 15, 15, 15, 20, 20, 20, 20),
                       c(18, 20, 25, 30, 35, 25, 30, 35, 40),
                     v=c(25, 26, 27, 28, 29, 30, 31, 32, 33))
deactivate.vertices(anet, c(15, 15, 15, 15, 15, 20, 20, 20, 20),
               length =c( 3,  5, 10, 15, 20,  5, 10, 15, 20),
                    v =c(35, 36, 37, 38, 39, 40, 41, 42, 43))
f7  = identical(anet$val[[25]]$active, anet$val[[35]]$active)
f8  = identical(anet$val[[26]]$active, anet$val[[36]]$active)
f9  = identical(anet$val[[27]]$active, anet$val[[37]]$active)
f10 = identical(anet$val[[28]]$active, anet$val[[38]]$active)
f11 = identical(anet$val[[29]]$active, anet$val[[39]]$active)
f12 = identical(anet$val[[30]]$active, anet$val[[40]]$active)
f13 = identical(anet$val[[31]]$active, anet$val[[41]]$active)
f14 = identical(anet$val[[32]]$active, anet$val[[42]]$active)
f15 = identical(anet$val[[33]]$active, anet$val[[43]]$active)


f.tests = paste("f", seq(1,15), sep="")
f.results= sapply(f.tests, function(x){eval(parse(text=x))})
if(any(!f.results)){
  bad.tests = paste("f", which(!f.results), sep="", collapse=" ")
  stop(paste("deactivate.vertices is incorrectly activating vertices in tests",
             bad.tests))
}

# testing deactivate edges associated with the deactivated vertex
# ---------------
# deactivate associated edges when deactivating a vertex
# all vertices at once
net <-network.initialize(3)
net[1,2]<-1;
net[2,3]<-1;
activate.edges(net,onset=1,terminus=Inf,e=1)
activate.edges(net,onset=2,terminus=3,e=2)
activate.vertices(net, onset=1, terminus=Inf)
deactivate.vertices(net, onset=2, terminus=3, deactivate.edges=TRUE)
dv.test <- as.data.frame(net)
if (!all(dv.test[,1:4] == matrix(c(1,2,1,2, 3,Inf,1,2), byrow=T, ncol=4))) {
  stop('deactivate.vertices did not perform as expected when deactivate.edges=TRUE')
}

# deactivate associated edges when deactivating a vertex
# one vertex at a time
net <-network.initialize(3)
net[1,2]<-1;
net[2,3]<-1;
activate.edges(net,onset=1,terminus=Inf,e=1)
activate.edges(net,onset=2,terminus=3,e=2)
activate.vertices(net, onset=1, terminus=Inf)
deactivate.vertices(net, onset=2, terminus=3, v=1, deactivate.edges=TRUE)
dv.test <- as.data.frame(net)
if (!all(dv.test[,1:4] == matrix(c(1,2,1,2, 3,Inf,1,2, 2,3,2,3), byrow=T, ncol=4))) {
  stop('deactivate.vertices did not perform as expected when deactivate.edges=TRUE')
}


# deactivate associated edges when deactivating a vertex
# no edge activity
net <-network.initialize(3)
net[1,2]<-1;
net[2,3]<-1;
activate.vertices(net, onset=1, terminus=Inf)
deactivate.vertices(net, onset=2, terminus=3, v=2, deactivate.edges=TRUE)
dv.test = as.data.frame(net)
if (!all(dv.test[,1:4] == matrix(c(-Inf,2,1,2, 3,Inf,1,2, -Inf,2,2,3, 3, Inf, 2, 3), byrow=TRUE, ncol=4))) {
  stop('deactivate.vertices did not perform as expected when deactivate.edges=TRUE')
}

# deactivate associated edges when deactivating a vertex
# no active edges (Inf, Inf) spells
net <-network.initialize(3)
net[1,2]<-1;
net[2,3]<-1;
activate.vertices(net, v=2, onset=1, terminus=Inf)
deactivate.edges(net, e=1, onset=-Inf, terminus=Inf)
deactivate.edges(net, e=2, onset=-Inf, terminus=Inf)
deactivate.vertices(net, onset=2, terminus=3, v=2, deactivate.edges=TRUE)
dv.test = as.data.frame(net)
if (nrow(dv.test[,1:4]) != 0) {
  stop('deactivate.vertices did not perform as expected when deactivate.edges=TRUE')
}

# deactivate associated edges when deactivating a vertex
# no active vertices. It will still run deactivate.edges on the vertices!
net <-network.initialize(3)
net[1,2]<-1;
net[2,3]<-1;
deactivate.vertices(net, onset=-Inf, terminus=Inf)
deactivate.vertices(net, onset=2, terminus=3, v=2, deactivate.edges=T)
dv.test = as.data.frame(net)
if (!all(dv.test[,1:4] == matrix(c(-Inf,2,1,2, 3,Inf,1,2, -Inf,2,2,3, 3, Inf, 2, 3), byrow=T, ncol=4))) {
  stop('deactivate.verteices did not perform as expected when deactivate.edges=TRUE')
}

expect_equal(network.size.active(deactivate.vertices(network.initialize(0)),onset=-Inf,terminus=Inf),0)

cat("ok\n")

#------------------- DELETE.ACTIVITY TESTS ------------------------
# Notes:
#  -- given the seriously limited time that I am working this
#     summer, this testing is minimal, at best. Feel free to add
#     to this if you have time.  --alc
#-----------------------------------------------------------------

cat("testing delete.*.activity ... ")

# a case with no edges
gnet1 = network.initialize(5)
delete.edge.activity(gnet1)
active1 = unlist(lapply(gnet1$val, "[[", "active"))
g1 = is.null(active1)
activate.vertices(gnet1, at=2)
delete.vertex.activity(gnet1, v=c(1,3,5))
active2 = unlist(lapply(gnet1$val, "[[", "active"))
g2 = (length(active2)==4 & all(active2==2))

# a case with missing edges
data(flo)
gnet3 = network(flo)
gnet3[9,16] = 0  # makes edge 39 NULL
activate.edges(gnet3, 2, 3)
delete.edge.activity(gnet3, e=1:37)
active3 = unlist(lapply(lapply(gnet3$mel, "[[", "atl"), "[[", "active"))
g3 = all(active3 == c(2,3,2,3))
delete.edge.activity(gnet3)
active4 = unlist(lapply(lapply(gnet3$mel, "[[", "atl"), "[[", "active"))
g4 = is.null(active4)

# a more typical case
gnet5 = network(flo)
activate.edges(gnet5, at = 1:40)
activate.vertices(gnet5, at = 101:117)
gnet6 = gnet7 = gnet5
delete.edge.activity(gnet5)
delete.vertex.activity(gnet5)
delete.vertex.activity(gnet6, v=1:8)
delete.edge.activity(gnet7, e=seq(2,40,2))
active5 = c(unlist(lapply(lapply(gnet5$mel, "[[", "atl"), "[[", "active")),
            unlist(lapply(gnet5$val, "[[", "active")))
active6 = c(unlist(lapply(lapply(gnet6$mel, "[[", "atl"), "[[", "active")),
            unlist(lapply(gnet6$val, "[[", "active")))
active7 = c(unlist(lapply(lapply(gnet7$mel, "[[", "atl"), "[[", "active")),
            unlist(lapply(gnet7$val, "[[", "active")))
g5 = (is.null(active5))
g6 = all(active6 == c(rep(1:40, each=2), rep(109:116, each=2))) 
g7 = all(active7 == c(rep(seq(1,39,2), each=2), rep(101:116, each=2)))

g.tests = paste("g", seq(1,7), sep="")
g.results= sapply(g.tests, function(x){eval(parse(text=x))})
if(any(!g.results)){
  bad.tests = paste("g", which(!g.results), sep="", collapse=" ")
  stop(paste("remove.activity is incorrectly removing activity matrices in tests",
             bad.tests))
}

expect_equal(network.size.active(delete.vertex.activity(network.initialize(0)),onset=-Inf,terminus=Inf),0)

expect_equal(network.edgecount.active(delete.edge.activity(network.initialize(0)),onset=-Inf,terminus=Inf),0)

# this triggered error for issue #298
net<-network.initialize(2)
activate.vertices(net,at=1)
activate.vertices(net,at=5)
deactivate.vertices(net,onset=-Inf,terminus=5)
expect_equal(as.numeric(get.vertex.activity(net)[[1]]),c(5,5))

# Assignment to networks as list elements
net <- list()
net[[1]] <- network.initialize(2)
activate.vertices(net[[1]],v=1,at=1)
activate.vertices(net[[1]],v=2,at=2)
expect_equal(
  c(
    network.vertex.names(network.extract(net[[1]],at=1)),
    network.vertex.names(network.extract(net[[1]],at=2))
  ),
  c(1,2)
)

cat("ok\n")

