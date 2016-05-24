bjtestII <- function(y, d, x, beta) {  ## depends on WKM( ), redistF( )
n <- length(y)
x <- as.matrix(x)
xdim <- dim(x)
if ( xdim[1] != n ) stop("check dim of x and y")
if ( length(beta) != xdim[2] ) stop("check dim of x and beta")

e <- y - as.vector( x %*% beta )
ordere <- order(e, -d)
esort <- e[ordere]
dsort <- d[ordere]
xsort <- as.matrix(x[ordere,])
dsort[length(dsort)] <- 1  #last one as uncensored always

## use KM as F (need to be an n vector prob)
temp0 <- WKM(esort,dsort, zc=1:n)
pKM <- temp0$jump
temp <- redistF(y=esort, d=dsort, Fdist=pKM) 
### what to do if there is tied obs.???
weight <- temp$weight/n   #the prob weight matrix

A <- matrix(0, ncol=xdim[2], nrow=n)
for (i in 1:n) if (dsort[i] == 1) { 
  A[i, ] <- t(as.matrix(weight[1:i,i])) %*% xsort[1:i, ]
}
####
Mi <- dsort*esort
Mi <- Mi*A
####

temp2 <- el.test(x=Mi, mu=rep(0, xdim[2]))
return(temp2)
}
