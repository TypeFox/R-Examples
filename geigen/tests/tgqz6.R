
library("geigen")
source("testqz.R")

set.seed(11)

n <- 5
A <- matrix(runif(n*n),nrow=n)+0i
B <- matrix(runif(n*n),nrow=n)+0i

# Test interface to zgges (QZ method)

z <- gqz(A, B,"-")
testqz(A,B,z)
gev <- gevalues(z)
all(which(Re(gev)<0) == seq_len(z$sdim))

z <- gqz(A, B,"+")
testqz(A,B,z)
gev <- gevalues(z)
all(which(Re(gev)>0) == seq_len(z$sdim))

z <- gqz(A, B,"S")
testqz(A,B,z)
gev <- gevalues(z)
all(which(abs(gev)<1) == seq_len(z$sdim))

z <- gqz(A, B,"B")
testqz(A,B,z)
gev <- gevalues(z)
all(which(abs(gev)>1) == seq_len(z$sdim))

z <- gqz(A, B,"R")
testqz(A,B,z) 
gev <- gevalues(z)
all(which(abs(Im(gev)) <= 100*.Machine$double.eps) == seq_len(z$sdim))
