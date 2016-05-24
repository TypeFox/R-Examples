##' Random generation noisy blockwise matrices
##'
##' Function to draw data.
##'
##' @param n number of rows and columns.
##' @param mu symetric matrix to the means.
##' @param sigma variance of the variables.
##' @param type represent the spacing between two change-point: "Eq" for a homogenous spacement, 
##'"NEq" for an arithmetic spacement and "NEqbis" for a decreasing arithmetic spacement.
##'
##' @name rblockdata
##' @rdname rblockdata
##'
##' @examples
##' ## model parameters 
##' n <- 100 
##' K <- 5
##' mu <- suppressWarnings(matrix(rep(c(1,0),ceiling(K**2/2)), K,K))
##' Y <- rblockdata(n,mu,sigma=.5)
##'
##' @export rblockdata
rblockdata <- function(n,mu,sigma,type=c("Eq","NEq","NEqbis")) {

    type <- match.arg(type)
    K <- ncol(mu)
    
    if (type == "Eq") {
        bt <- rep(n/K,K)
        empl <- c(0,cumsum(bt))
    } else if (type == "NEqbis") {
        empl <- cumsum(floor(n/sum(1:K)*(1:K)))
        empl[K] <- n
        empl <- c(0,empl)
        bt <- empl[2:(K+1)]-empl[1:K]        
    } else if (type == "NEqbis") {
        empl <- cumsum(floor(n/sum(K:1)*(K:1)))
        empl[K] <- n
        empl <- c(0,empl)
        bt <- empl[2:(K+1)]-empl[1:K]
    }

    Y <- matrix(0,n,n)
    betastar <- Matrix(0,n,n,sparse=TRUE)
    for (k in 1:K){
        for (l in k:K){
            Y[(empl[k]+1):empl[k+1],(empl[l]+1):empl[l+1]]=
                matrix(c(rnorm(bt[k]*bt[l],mu[k,l],sd=sigma)),nrow=bt[k])
            if (k!=l){
                Y[(empl[l]+1):empl[l+1],(empl[k]+1):empl[k+1]]=
                    Y[(empl[k]+1):empl[k+1],(empl[l]+1):empl[l+1]]
            }
            if ((k==1)&&(l==1)){
                betastar[empl[k]+1,empl[l]+1]=mu[1,1]
            } else{
                betastar[empl[k]+1,empl[l]+1]=mu[k,l]-mu[k,l-1]
                if (k!=l){
                    betastar[empl[l]+1,empl[k]+1]=betastar[empl[k]+1,empl[l]+1]
                }
            }            
        }
    }
    return(list(Y = as.matrix(forceSymmetric(Y)),
                betastar = betastar, bt = bt, empl = empl))
}

## ##' MODEL:  Y = U + E with U = T B T'
## ##' with T lower triangular and
## ##' B block wise such that each block has the pattern
## ##'  (bij  ... 0 )
## ##'  ( 0   ... 0 )
## ##'  ( .       . )
## ##'  ( .       . )
## ##'  ( 0       0 )
## ##' t.star defined the blocs in B (changes)
## rHiCdata <- function(n, t.star = c(1,floor(n/4),floor(n/2),floor(3*n/4))) {

##     ##  B matrix
##     B <- Matrix(0,n,n)

##     ## define the bloc changes
##     nb.blocs <- length(t.star)

##     dim_mu_tr <- nb.blocs*(nb.blocs-1)/2
##     mat_mu=matrix(0,nb.blocs,nb.blocs)
##     mat_mu[lower.tri(mat_mu,diag=FALSE)]=floor(runif(dim_mu_tr,5,20))
##     mat_mu=mat_mu+t(mat_mu)
##     diag(mat_mu)=c(5,10,15,20)
##     for (i in 1:length(t.star)) {
##         for(j in 1:length(t.star)) {
##             B[t.star[i],t.star[j]]=mat_mu[i,j]
##         }
##     }
    
##     ## Noise matrix (symmetric)
##     E <- Matrix(0,n,n)
##     E[upper.tri(E, diag=TRUE)] <- rnorm((n+n*(n-1)/2),0,1)
##     E <- forceSymmetric(E)
    
##     ## Triangular low (with bandSparse : extremely small amount of memory)
##     ## T <- bandSparse(n,k=(-n+1):0,diagonals=matrix(rep(1,n*(n+1),2),n,n))
##     ## T <- bandSparse(n,k=(-n+1):0) * 1
##     T <- Matrix(0,n,n,sparse=TRUE)
##     T[lower.tri(T,diag=TRUE)] <- 1
    
##     ## Observed Matrix
##     Y <- T %*% B %*% t(T) + E

##     return(list(Y=as.matrix(Y),B=B))
## }

