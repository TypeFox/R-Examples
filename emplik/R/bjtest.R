bjtest <- function(y, d, x, beta) {  ## depends on WKM( ), redistF( )
n <- length(y)
x <- as.matrix(x)
xdim <- dim(x)
if ( xdim[1] != n ) stop("check dim of x")
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
temp <- redistF( y=esort, d=dsort, Fdist=pKM ) 
### what to do if there is tied obs.???
weight <- temp$weight/n   #the prob weight matrix

A <- matrix(0, ncol=xdim[2], nrow=n)
for (i in 1:n) if (dsort[i] == 1) { 
  A[i, ] <- t(as.matrix(weight[1:i,i])) %*% xsort[1:i, ]
  A[i, ] <- A[i,]/pKM[i] 
}

myfun <- function(t, q) { t*q }

temp2 <- el.cen.EM2(x=esort, d=dsort, fun=myfun, mu=rep(0, xdim[2]), q=A)
pnew <- temp2$prob
logel1 <- temp0$logel
logel2 <- temp2$loglik 

list(prob=pnew, logel=logel1, logel2=logel2, "-2LLR"=2*(logel1-logel2))
}
