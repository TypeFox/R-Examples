bjtest1d <- function(y, d, x, beta) {  ## depends on WKM( ), redistF( )
n <- length(y)                         ## dimension of x must be n x 1.
if ( length(x) != n ) stop("check dim of x")
if ( length(beta) != 1 ) stop("check dim of beta")

e <- y - beta * x 
ordere <- order(e, -d)
esort <- e[ordere]
dsort <- d[ordere]
xsort <- x[ordere]
dsort[length(dsort)] <- 1  #last one as uncensored always

## use KM as F (need to be an n vector prob)
temp0 <- WKM(esort,dsort)
pold <- temp0$jump
temp <- redistF( y=esort, d=dsort, Fdist=pold ) 
weight <- temp$weight/n   #the prob weight matrix
####pold <- colSums(weight)  ##not needed, just let pold= WKM()$jump

A <- rep(0, n)
for (i in 1:n) if (dsort[i] == 1) { 
  A[i] <- sum( weight[1:i,i] * xsort[1:i] )/pold[i] 
}

AA <- A[ dsort == 1 ]
myfun <- function(t, q) { t*q }

temp2 <- el.cen.EM(x=esort, d=dsort, fun=myfun, mu=0, q=AA)
pnew <- temp2$prob
logel1 <- temp0$logel
logel2 <- temp2$loglik 

list(prob=pnew, logel=logel1, logel2=logel2, "-2LLR"=2*(logel1-logel2))
}
