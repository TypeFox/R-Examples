
library("geigen")
source("testqz.R")

set.seed(11)

n <- 5
A <- matrix(runif(n*n),nrow=n)+0i
B <- matrix(runif(n*n),nrow=n)+0i

B[5,] <- (B[4,]+B[3,])/2

# Test interface to zgges (QZ method)

z <- gqz(A, B,"N")
testqz(A,B,z)
gev <- gevalues(z)
z$sdim == 0

z <- gqz(A, B,"-")
testqz(A,B,z)
gev <- gevalues(z)
all(which(ifelse(is.finite(gev),Re(gev)<0,FALSE)) == seq_len(z$sdim))

z <- gqz(A, B,"+")
testqz(A,B,z)
gev <- gevalues(z)
all(which(ifelse(is.finite(gev),Re(gev)>0,FALSE)) == seq_len(z$sdim))

z <- gqz(A, B,"S")
testqz(A,B,z)
gev <- gevalues(z)
all(which(ifelse(is.finite(gev),abs(gev)<1,FALSE)) == seq_len(z$sdim))

z <- gqz(A, B,"B")
testqz(A,B,z)
gev <- gevalues(z)
all(which(ifelse(is.finite(gev),abs(gev)>1,FALSE)) == seq_len(z$sdim))

z <- gqz(A, B,"R")
testqz(A,B,z)
gev <- gevalues(z)
all(which(ifelse(is.finite(gev),abs(Im(gev)) <= 100*.Machine$double.eps,FALSE)) == seq_len(z$sdim))
