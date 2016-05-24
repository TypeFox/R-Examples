RankRegTest <- function(y, d, x, beta, type="Gehan") {  ## depends on WKM( )
n <- length(y)                          ## dimension of x must be n x q.
x <- as.matrix(x)                       ## x must NOT including an intercept.
xdim <- dim(x)
if ( xdim[1] != n ) stop("check dim of x")
if ( length(beta) != xdim[2] ) stop("check dim of beta and x")
if(any((d!=0)&(d!=1)))
     stop("d must be 0(right-censored) or 1(uncensored)")
if (length(d) != n) stop("check the length of d")

e <- y - as.vector( x %*% beta ) 
ordere <- order(e, -d)
esort <- e[ordere]
dsort <- d[ordere]
xsort <- as.matrix(x[ordere,])
dsort[length(dsort)] <- 1       #last one as uncensored always

## compute KM  (need to be an n vector prob)
temp0 <- WKM(esort,dsort, zc=1:n)
pKM <- temp0$jump

##xbar <- rev(cumsum(rev(xsort)))/(n:1)

xbar <- xsort
####for(j in 1:(n-1)) xbar[j,] <- colMeans(xsort[j:n,])
for(j in 1:xdim[2]) xbar[,j] <- cumsumsurv(xsort[,j])/(n:1)  ## rev(cumsum(rev(xsort[,j])))/(n:1) 3/2015 MZ

if(type == "Gehan") {A <- (n:1) * (xsort - xbar)/pKM}
 else {if(type == "Logrank") A <- (xsort - xbar)/pKM
        else stop("type must be either Gehan or Logrank") }
####AA <- as.matrix(A[dsort == 1,])

myfun <- function(t, q) { q }

temp2 <- el.cen.EM2(x=esort, d=dsort, fun=myfun, mu=rep(0,xdim[2]), q=A)
pnew <- temp2$prob
logel1 <- temp0$logel
logel2 <- temp2$loglik 

list(prob=pnew, logel=logel1, logel2=logel2, "-2LLR"=2*(logel1-logel2))
} 
#
# Use empirical likelihood to test in the AFT model rank estimator.
#
