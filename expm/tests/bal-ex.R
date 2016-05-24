library(expm)
source(system.file("test-tools.R", package= "expm"), keep.source=FALSE)## -> assertError()...

## A matrix with 'Inf'
mI <- rbind(0, c(-Inf, Inf, 0, 0), 0, 0)
bal3 <-
    list(dB = dgebal(mI, "B"), # = default
         dP = dgebal(mI, "P"),
         dN = dgebal(mI, "N"))
str(bal3)
stopifnot(identical(mI, bal3$dN$z),
          with(bal3, all.equal(dB, dP, tol=1e-14)),
          all.equal(bal3$dB$z, rbind(c(Inf,-Inf,0,0), 0,0,0), tol=1e-14),
          all.equal(bal3$dB$scale, c(1,1,3,4)))
assertError(dgebal(mI, "S"))# gave infinite loop



## Compare the two different "balance" pre-conditioning versions in Ward77:
set.seed(1)
mList <- lapply(integer(100), function(...) rSpMatrix(20, nnz=80))
re20 <- sapply(mList, function(M)
               relErr(expm(M, precond = "2bal"),
                      expm(M, precond = "1bal")))
re20 ## ahh.: zero or ~ 1e-13 ... good
table(re20 == 0)
summary(re20[re20 != 0])
## Pentium M (ubuntu)
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
## 2.593e-14 8.703e-14 1.282e-13 2.434e-13 4.177e-13 6.295e-13


demo(balanceTst) #-> the function definition and the first few examples

dm4. <- dgebal(m4)
storage.mode(m4) <- "integer"
stopifnot(identical(dm4., dgebal(m4)))

expm(m)
expm(m,"Pade") ## are different indeed {when bug still existed}
expm(m,"R_Pade")# same as Pade


## a non-empty ``non-balanced'' example  ---

expm.t.identity(m4, "Ward")

m6 <- zeroTrace(matrix(outer(2^(-8:9),c(-1,1)), 6,6)); m6
m6[lower.tri(m6)] <- 0 ## plus one non-zero
m6[4,2] <- 77
p <- c(6,4,5,2:1,3); m6 <- m6[p,p]
expm.t.identity(m6, "Ward") ##  difference; indeed
expm(m6) # is very different from
expm(m6,"R_Pade")

str(dm6 <- dgebalTst(m6))
## Now, that's interesting:
##
## 1.  'S' scales *more* (2 .. 5) than just (2:4 == i1:i2) !
##
## 2.  'B' has quite different scaling and it does (must!) obey rule
##     scale i1:i2 only
##
## 3. 'B'(oth) is better than  "P" and "S" separately:
##
kappa(eigen(m6)$vectors)#      597.5588
kappa(eigen(dm6$P$z)$vectors)# 597.5588
kappa(eigen(dm6$S$z)$vectors)#  42.58396
kappa(eigen(dm6$B$z)$vectors)#  22.20266


## An n=17 example where octave's expm() is wrong too
m17 <- matrix(c(10,0, 0, 2, 3,-1, 0, 0, 0, 0, 0, 4, 0, 5, 0, 0,-2,
                0, 0, 0, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 7, 0,
                0, 0,10, 0, 0,-4, 9, 0, 0, 0,-5, 0,-6, 0, 0, 0, 0,
                0, 0,-7, 0, 0, 0, 0, 0, 0,10, 0, 0, 0, 0, 0,11, 0,
                0, 0, 0, 0, 0, 0,12, 0, 0, 0, 0, 0,-8, 0, 0, 0, 0,
                0, 0,-9, 0, 0, 0, 0, 0, 0,-10,0,13,14,-11,-12,-13, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0,10, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0,-14,16,0,-10,0,17, 0, 0, 0, 0, 0, 0,
                0, 0,-16,0, 0,18,19, 0, 0, 0, 0, 0, 0, 0,20, 0, 21,
                22,0, 0, 0, 0, 0,-17,0, 0, 0,-10,-19,-20,0,0,0, 0,
                0,-21,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0,23,24, 0,25,26, 0, 0,27,-22,0,28,-23,0,-24,
                0,-25,0,29, 0, 0, 0, 0, 0, 0, 0,30,31, 0, 0, 0, 0,
                0, 0,-26,32,0, 0, 0, 0, 0,-27,0,33,34, 0, 0, 0, 0,
                0,-28,-29,0,0, 0,35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0,36,37, 0, 0, 0, 0, 0, 0, 0, 0,-10),
              17, 17)
str(dm17 <- dgebalTst(m17))
sapply(dm17[1:3], `[[`, "scale")

## The balancing was really rather harmful -- cond(V) *not* improved:
condX <- function(x) kappa(x, exact=TRUE)
condX(eigen(m17)$vectors)#      8.9e16
condX(eigen(dm17$P$z)$vectors)# 1.37e17
condX(eigen(dm17$S$z)$vectors)# 1.44e17
condX(eigen(dm17$B$z)$vectors)# 1.43e17 (very slightly smaller)
