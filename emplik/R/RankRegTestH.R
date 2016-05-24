RankRegTestH <- function(y, d, x, beta, type="Gehan") {
     n <- length(y)          ## dimension of x must be n x q.
     x <- as.matrix(x)       ## x must NOT including an intercept.
     xdim <- dim(x)
     if( xdim[1] != n ) stop("check dim of x")
     if( length(beta) != xdim[2] ) stop("check dim of beta and x")
 
e <- y - as.vector( x %*% beta )
ordere <- order(e, -d)
esort <- e[ordere]
dsort <- d[ordere]
xsort <- as.matrix(x[ordere,])
dsort[length(dsort)] <- 1       #last one as uncensored always?
 
##xbar <- rev(cumsum(rev(xsort)))/(n:1)

xbar <- xsort
####for(j in 1:(n-1)) xbar[j,] <- colMeans(xsort[j:n,])
for(j in 1:xdim[2]) xbar[,j] <- cumsumsurv(xsort[,j])/(n:1)   ##  rev(cumsum(rev(xsort[,j])))/(n:1)  3/2015 MZ

if(type == "Gehan") {A <- (n:1)^2 * (xsort - xbar)}
 else {if(type == "Logrank") A <- (n:1) * (xsort - xbar)
        else stop("type must be either Gehan or Logrank") }

## A1 <- (n:1)^2 * (xsort -xbar)   #  for over -determine case
## A2 <- (n:1) * (xsort - xbar)    #  both constraints
## A <- cbind(A1, A2)              #
## un-comment the above 3 lines to get over determined estimator
AA <- as.matrix(A[dsort == 1,])

myfun <- function(t, A){ return(A) }
myfun2 <- function(t){ matrix(0, ncol=ncol(AA), nrow=length(t)) }
x20 <- 1:30        ## fake data for sample two
d20 <- rep(1, 30)  ## fake data for sample two

temp2 <- emplikHs.test22(x1=esort, d1=dsort, x2=x20, d2=d20, theta=rep(0,ncol(AA)), fun1=myfun, fun2=myfun2, A=AA)

Samp1EL <- temp2$"-2LLR(sample1)"

list(loglikH0=temp2$"Llik(sample1)",  "-2LLR"=Samp1EL)
}
