## private communication from Gwen (gleday <gleday@few.vu.nl>) via Kurt Hornik
## produced an endless loop on some platforms.

library(quadprog)
load("bug.RData")
sol <- solve.QP(Dmat, as.vector(dvec), Amat, bvec, meq=meq)
print(lapply(sol, zapsmall))
