redistF <- function(y, d, Fdist) {
n <- length(d)
if ( length(y) != n ) stop("length of y and d must agree")
if ( any((d != 0) & (d !=1)) ) stop("d must be either 0 or 1")
if ( length(Fdist) != n ) stop("Fdist must have length n")
# Fdist must be a probability distribution
# the time vector is the ysort, and the prob mass is in Fdist

yorder <- order(y, -d)
ysort <- y[yorder]
dsort <- d[yorder]

WeightMat <- diag( rep(1, n) )

for (i in 1:(n-1)) 
     if ( dsort[i] == 0 ) {
           WeightMat[i, (i+1):n] <- Fdist[(i+1):n]/sum(Fdist[(i+1):n])
           WeightMat[i,i] <- 0
                          }

list(y=ysort, d=dsort, weight=WeightMat, ordY=yorder)
}
