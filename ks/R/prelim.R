##############################################################################
## Parse variable name
############################################################################

parse.name <-function(x)
{
  if (is.vector(x))
  {
    d <- 1
    x.names <- deparse(substitute(x))
  }
  else
  {  
    d <- ncol(x)
    x.names <- colnames(x)
    if (is.null(x.names))
    {
      x.names <- strsplit(deparse(substitute(x)), "\\[")[[1]][1]
      x.names <- paste(x.names, "[, ", 1:d,"]",sep="") 
    }
  }
  return(x.names)
}



#############################################################################
## Basic vectors and matrices and their operations
#############################################################################
  
## Vec operator 

vec <- function(x, byrow=FALSE)
{
  if (is.vector(x)) return (x)
  
  if (byrow) x <- t(x)
  d <- ncol(x)
  vecx <- vector()
  for (j in 1:d)
    vecx <- c(vecx, x[,j])
  
  return(vecx)           
}

## Vech operator

vech <- function(x)
{
  if (is.vector(x))
  {
    if (length(x)==1)
      return (x)
    else
      stop("vech undefined for vectors")
  }
  else if (is.matrix(x))
  {  
    d <- ncol(x)
    if (d!=nrow(x)) ##if (!isSymmetric(x))
      stop("vech only defined for square matrices")
    
    vechx <- vector()
    for (j in 1:d)
      vechx <- c(vechx, x[j:d,j])
    return(vechx)
  }
  
}

  
## Inverse vec operator 
    
invvec <- function(x, ncol, nrow, byrow=FALSE)
{
  if (length(x)==1)
    return(x)
  
  d <- sqrt(length(x))
  if (missing(ncol) | missing(nrow))
  {
    ncol <- d; nrow <- d
    if (round(d) != d)
      stop("Need to specify nrow and ncol for non-square matrices")
  }
  
  invvecx <- matrix(0, nrow = nrow, ncol = ncol)
  if (byrow)
    for (j in 1:nrow)
      invvecx[j,] <- x[c(1:ncol) + (j-1)*ncol]
  else
    for (j in 1:ncol)
      invvecx[,j] <- x[c(1:nrow) + (j-1)*nrow]
  
  return(invvecx)
}

## Inverse vech operator 

invvech <- function(x)
{
  if (length(x)==1)
    return(x)
  
  d <- (-1 + sqrt(8*length(x) + 1))/2
  if (round(d) != d)
    stop("Number of elements in x will not form a square matrix")
  invvechx <- matrix(0, nrow=d, ncol=d)

  for (j in 1:d)
    invvechx[j:d,j] <- x[1:(d-j+1)+ (j-1)*(d - 1/2*(j-2))]
  
  invvechx <- invvechx + t(invvechx) - diag(diag(invvechx))
  
  return(invvechx)
}

## Trace of matrix

tr <- function(A)
{
  count <- 0
  if (is.vector(A)) return (A[1])
  if (nrow(A)!=ncol(A))
    stop("Not square matrix")

  else 
     for (i in 1:nrow(A))
       count <- count + A[i,i]
  return(count)
}


## Elementary vector 
    
elem <- function(i, d)
{
  elem.vec <- rep(0, d)
  elem.vec[i] <- 1
  
  return(elem.vec)
}      

## Commutation matrix (taken from MCMCglmmm library)

comm <- function(m,n){
  K<-matrix(0,m*n, m*n)
  H<-matrix(0,m,n)
  for(i in 1:m){
    for(j in 1:n){ 
      H[i,j]<-1
      K<-K+kronecker(H,t(H))
      H[i,j]<-0
    }
  }
  return(K)
}

###############################################################################
## Duplication matrix
## Taken from Felipe Osorio http://www.ime.usp.br/~osorio/files/dupl.q
###############################################################################

dupl <- function(order, ret.q = FALSE)
{
    ## call
    cl <- match.call()
    time1 <- proc.time()
    if (!is.integer(order))
        order <- as.integer(order)
    n <- order - 1
    
    ## initial duplication matrix
    d1 <- matrix(0, nrow = 1, ncol = 1)
    d1[1,1] <- 1
    if (!is.integer(d1))
        storage.mode(d1) <- "integer"
    
    ## recursive formula
    if (n > 0){
    	for (k in 1:n){
    	    drow <- 2*k + 1 + nrow(d1)
    	    dcol <- k + 1 + ncol(d1)
    	    d2 <- matrix(0, nrow = drow, ncol=dcol)
    	    storage.mode(d2) <- "integer"
    	    d2[1,1] <- 1
    	    d2[2:(k+1),2:(k+1)] <- diag(k)
    	    d2[(k+2):(2*k+1),2:(k+1)] <- diag(k)
    	    d2[(2*k+2):drow,(k+2):dcol] <- d1
    	    ## permutation matrix
    	    q <- permute.mat(k)
    	    ## new duplication matrix
    	    d2 <- q %*% d2
    	    storage.mode(d2) <- "integer"
    	    d1 <- d2
    	}
    }
    else {
    	d2 <- q <- d1
    }
    
    ## results
    obj <- list(call=cl, order=order, d=d2)
    if (ret.q)
        obj$q <- q
    obj$time <- proc.time() - time1
    obj
}


###############################################################################
## Pre-scaling
## Parameters
## x - data points
##
## Returns
## Pre-scaled x values
###############################################################################

pre.scale <- function(x, mean.centred=FALSE)
{
  S <- diag(diag(var(x)))
  Sinv12 <- matrix.sqrt(chol2inv(chol(S)))

  if (mean.centred) x.scaled <- sweep(x, 2, apply(x, 2, mean))
  else x.scaled <- x  
  x.scaled <- x.scaled %*% Sinv12
  
  return (x.scaled)
}

###############################################################################
## Pre-sphering
## Parameters
## x - data points
##
## Returns
## Pre-sphered x values
###############################################################################

pre.sphere <- function(x, mean.centred=FALSE)
{
  S <- var(x)
  Sinv12 <- matrix.sqrt(chol2inv(chol(S)))

  if (mean.centred) x.sphered <- sweep(x, 2, apply(x, 2, mean))
  else x.sphered <- x  
  x.sphered <- x.sphered %*% Sinv12

  return (x.sphered)
}

##############################################################################
## Boolean functions
###############################################################################

is.even <- function(x)
{
  y <- x[x>0] %%2
  return(identical(y, rep(0, length(y))))
}

is.diagonal <- function(x)
{
  return(identical(diag(diag(x)),x))
}

###############################################################################
## Finds row index matrix
## Parameters
## x - data points
##
## Returns
## i  - if r==mat[i,]
## NA - otherwise
###############################################################################

which.mat <- function(r, mat)
{
  ind <- numeric()
  
  for (i in 1:nrow(mat))
    if (identical(r, mat[i,])) ind <- c(ind,i)

  return(ind)  
}



###################################################################
## Permutation functions
###################################################################

####################################################################
## Exactly the same function as combinat:::permn
####################################################################

permn.ks <- function (x, fun = NULL, ...) 
{
    if (is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x) 
      x <- seq(x)
    n <- length(x)
    nofun <- is.null(fun)
    out <- vector("list", gamma(n + 1))
    p <- ip <- seqn <- 1:n
    d <- rep(-1, n)
    d[1] <- 0
    m <- n + 1
    p <- c(m, p, m)
    i <- 1
    use <- -c(1, n + 2)
    while (m != 1) {
        out[[i]] <- if (nofun) 
            x[p[use]]
        else fun(x[p[use]], ...)
        i <- i + 1
        m <- n
        chk <- (p[ip + d + 1] > seqn)
        m <- max(seqn[!chk])
        if (m < n) 
            d[(m + 1):n] <- -d[(m + 1):n]
        index1 <- ip[m] + 1
        index2 <- p[index1] <- p[index1 + d[m]]
        p[index1 + d[m]] <- m
        tmp <- ip[index2]
        ip[index2] <- ip[m]
        ip[m] <- tmp
    }
    out
}

##########################################################################
## Permutations with repetitions of the first d naturals (1:d) taking 
## k elements at a time. There are d^k of them, each having length k 
## => We arrange them into a matrix of order d^k times k
## Each row represents one permutation
## Second version: filling in the matrix comlumn-wise (slightly faster)
##########################################################################

perm.rep<-function(d,r)
{
    if(r==0){PM<-1}
    if(r>0){
    PM<-matrix(nrow=d^r,ncol=r)
    for(pow in 0:(r-1)){
        t2<-d^pow
        p1<-1
        while(p1<=d^r){
            for(al in 1:d){for(p2 in 1:t2){
                PM[p1,r-pow]<-al
                p1<-p1+1}}}}
    }
    return(PM)
}

###############################################################################
## Permute a list of values
##
## Same function as EXPAND.GRID (base package), modified to take 
## list as an argument and returns a matrix 
###############################################################################

permute <- function (args) 
{
  nargs <- length(args)
  if (!nargs) 
    return(as.data.frame(list()))
  if (nargs == 1 && is.list(a1 <- args[[1]])) 
    nargs <- length(args <- a1)
  if (nargs <= 1) 
    return(as.data.frame(if (nargs == 0 || is.null(args[[1]])) list() else args, 
                         optional = TRUE))
  cargs <- args
  rep.fac <- 1
  orep <- prod(sapply(args, length))
  
  for (i in 1:nargs) {
    x <- args[[i]]
    nx <- length(x)
    orep <- orep/nx
    cargs[[i]] <- rep(rep(x, rep(rep.fac, nx)), orep)
    rep.fac <- rep.fac * nx
  }
  do.call("cbind", cargs)
} 

permute.mat <- function(order)
{
    m <- as.integer(order)
    m <- m + 1
    eye <- diag(m)
    u1 <- eye[1:m,1]
    u2 <- eye[1:m,2:m]
    q1 <- kronecker(eye, u1)
    q2 <- kronecker(eye, u2)
    q <- matrix(c(q1, q2), nrow = nrow(q2), ncol = ncol(q1) + ncol(q2))
    if (!is.integer(q))
        storage.mode(q) <- "integer"
    q
}

##########################################################################
### pinv.all generates all the permutations PR_{d,r} as described in
### Appendix B of Chacon and Duong (2014)
##########################################################################
    
pinv.all<-function(d,r){
    i<-1:d^r
    n<-i-1
    dpow<-d^(0:r)
    n.mat<-matrix(rep(n,r+1),byrow=FALSE,nrow=d^r,ncol=r+1)
    dpow.mat<-matrix(rep(dpow,d^r),byrow=TRUE,nrow=d^r,ncol=r+1)
    ndf.mat<-floor(n.mat/dpow.mat)
    ans<-ndf.mat[,r:1]-d*ndf.mat[,(r+1):2]    
    return(ans+1)
    } 


##############################################################################
## Block indices for double sums
##############################################################################

block.indices <- function(nx, ny, d, r=0, diff=FALSE, block.limit=1e6, npergroup)
{
  if (missing(npergroup)) 
  { 
    if (diff) npergroup <- max(c(block.limit %/% (nx*d^r), 1))
    else npergroup <- max(c(block.limit %/% nx,1))
  }
  nseq <- seq(1, ny, by=npergroup)
  if (tail(nseq,n=1) <= ny) nseq <- c(nseq, ny+1)
  if (length(nseq)==1) nseq <- c(1, ny+1)
  return(nseq)
}

block.indices2 <- function(nx, ny, block.limit=1e6, npergroup)
{
  if (missing(npergroup)) npergroup <- max(c(block.limit %/% nx,1))
  nseq <- seq(1, ny, by=npergroup)
  if (tail(nseq,n=1) <= ny) nseq <- c(nseq, ny+1)
  if (length(nseq)==1) nseq <- c(1, ny+1)
  return(nseq)
}


####################################################################
## Differences for double sums calculations
####################################################################

differences <- function(x, y, upper=FALSE, ff=FALSE, Kpow=0)
{
  if (missing(y)) y <- x
  if (is.vector(x)) x <- t(as.matrix(x))
  if (is.vector(y)) y <- t(as.matrix(y))

  nx <- nrow(x)
  ny <- nrow(y)
  d <- ncol(x)

  if (ff) difs <- ff(init=0, dim=c(nx*ny,d))
  else difs <- matrix(ncol=d,nrow=nx*ny)

  for (j in 1:d)
  {
    difs[,j] <- rep(x[,j], times=ny) - rep(y[1:ny,j], each=nx)
    ## jth column of difs contains all the differences X_{ij}-Y_{kj}
  }
 
  if (upper)
  {
    ind.remove <- numeric()
    for (j in 1:(nx-1))
      ind.remove <- c(ind.remove, (j*nx+1):(j*nx+j))
      
    return(difs[-ind.remove,])
  }
  else
    return(difs)
}


##### Odd factorial

OF<-function(m){factorial(m)/(2^(m/2)*factorial(m/2))} 

###############################################################################
## Matrix square root - taken from Stephen Lake 
## http://www5.biostat.wustl.edu/s-news/s-news-archive/200109/msg00067.html
###############################################################################

matrix.sqrt <- function(A)
{
  if (length(A)==1)
    return(sqrt(A))
  sva <- svd(A)
  if (min(sva$d)>=0)
    Asqrt <- sva$u %*% diag(sqrt(sva$d)) %*% t(sva$v)
  else
    stop("Matrix square root is not defined")
  return(Asqrt)
}

###############################################################################
## Matrix power
###############################################################################

matrix.pow <- function(A, n)
{
  if (nrow(A)!=ncol(A)) stop("A must be a square matrix")
  if (floor(n)!=n) stop("n must be an integer")
  if (n==0) return(diag(ncol(A)))
  if (n < 0) return(matrix.pow(A=chol2inv(chol(A)), n=-n))
        
  # trap non-integer n and return an error
  if (n == 1) return(A)
  result <- diag(1, ncol(A))
  while (n > 0) {
    if (n %% 2 != 0) {
      result <- result %*% A
      n <- n - 1
    }
    A <- A %*% A
    n <- n / 2
  }
  return(result)
}


##########################################################################
### Kmat computes the commutation matrix of orders m,n
##########################################################################    
    
Kmat<-function(m,n){ 
    K<-matrix(0,ncol=m*n,nrow=m*n)
    i<-1:m;j<-1:n
    rows<-rowSums(expand.grid((i-1)*n,j))
    cols<-rowSums(expand.grid(i,(j-1)*m))    
    positions<-cbind(rows,cols)
    K[positions]<-1
    return(K)
}

##########################################################################
### mat.Kprod computes row-wise Kronecker products of matrices
########################################################################## 

mat.Kprod<-function(U,V){ #### Returns a matrix with rows U[i,]%x%V[i,]
  n1<-nrow(U)

  n2<-nrow(V)
  if(n1!=n2)stop("U and V must have the same number of vectors")
  p<-ncol(U)
  q<-ncol(V)
  onep<-rep(1,p)
  oneq<-rep(1,q)
  P<-(U%x%t(oneq))*(t(onep)%x%V)
  return(P)
}


##########################################################################
## Kpow computes the Kronecker power of a matrix A
##########################################################################

Kpow<-function(A,pow){    
  if(floor(pow)!=pow)stop("pow must be an integer")
  Apow<-A
  if(pow==0){Apow<-1}
  if(pow>1){
    for(i in 2:pow) Apow<-Apow%x%A
  }
  return(Apow)
} 


#### Kronecker sum

Ksum <- function(A,B)
{
  AB <- numeric()
  for (i in 1:nrow(A))
    for (j in 1:nrow(B))
      AB <- rbind(AB, A[i,] + B[j,])

  return(AB)
}

#### Returns a matrix with the pow-th Kronecker power of A[i,] in the i-th row

mat.Kpow<-function(A,pow){ 
  Apow<-A
  if(pow==0){Apow<-matrix(1,nrow=nrow(A), ncol=1)}  
  if(pow>1){
    for(i in 2:pow) Apow<-mat.Kprod(Apow,A)
  }
  return(Apow)
}


#### Vector of all r-th partial derivatives of the normal density at x=0, i.e., D^{\otimes r)\phi(0)

DrL0 <- function(d,r)
{
  v <- as.vector(Kpow(A=vec(diag(d)),pow=r/2))
  DL0<-(-1)^(r/2)*(2*pi)^(-d/2)*OF(r)*matrix(Sdrv(d=d, r=r, v=v), ncol=1) 
  return(DL0)
}



#########################################################################
### Wrapper functions for Chacon & Duong (2014) 
##########################################################################

Sdr<-function(d, r, type="recursive"){
  type1 <- match.arg(type, c("recursive", "direct"))
  Sdr.mat <- do.call(paste("Sdr", type1, sep="."), list(d=d, r=r))
  return(Sdr.mat)
}

Sdrv<-function(d, r, v, type="recursive"){
  type1 <- match.arg(type, c("recursive", "direct"))
  v <- as.vector(v)
  Sdrvec <- do.call(paste("Sdrv", type1, sep="."), list(d=d, r=r, v=v))
  return(Sdrvec)
}


##########################################################################
## Symmetriser matrix
##########################################################################


############################################################################
### Sdr.direct computes the symmetrizer matrix S_{d,r} based on Equation (4)
### as described in Section 3 of Chacon and Duong (2014)
############################################################################

Sdr.direct<-function(d,r){ 
    S<-matrix(0,ncol=d^r,nrow=d^r)
    per<-permn.ks(r)
    per.rep<-pinv.all(d,r)
    nper<-factorial(r)
    nper.rep<-d^r
    per<-matrix(unlist(per), byrow=TRUE, ncol=r, nrow=nper)    
    pow<-0:(r-1)
    dpow<-d^pow
    
    if(nper.rep<=nper){
    dpow.mat<-matrix(rep(dpow,nper),byrow=TRUE,ncol=r,nrow=nper)
    for(i in 1:nper.rep){      ## Loop over no. perms with reps (d^r)
        pinvi<-per.rep[i,]
        sigpinvi<-matrix(pinvi[per],byrow=FALSE,nrow=nrow(per),ncol=ncol(per))
        psigpinvi<-drop(1+rowSums((sigpinvi-1)*dpow.mat))
        S[i,]<-tabulate(psigpinvi,nbins=d^r)
    }}
    
    if(nper<nper.rep){
    dpow.mat<-matrix(rep(dpow,nper.rep),byrow=TRUE,ncol=r,nrow=nper.rep)
    for(s in 1:nper){          ## Loop over no. perms (r!)
        sig<-per[s,]
        sigpinv<-per.rep[,sig]
        psigpinv<-drop(1+rowSums((sigpinv-1)*dpow.mat))
        locations<-cbind(1:d^r,psigpinv)
        S[locations]<-S[locations]+1
    }}    
    return(S/nper)
}
     
############################################################################
### Sdr.recursive computes the symmetrizer matrix S_{d,r} based on 
### the recursive approach detailed in Algorithm 1 in Section 3 of 
### Chacon and Duong (2014)
############################################################################     
     
Sdr.recursive<-function(d,r){
    S<-diag(d)
    if(r==0)S<-1
    if(r>=2){
        Id<-diag(d)
        T<-Id 
        A<-Kmat(d,d)        
        for(j in 2:r){
            T<-((j-1)/j)*(A%*%(T%x%Id)%*%A)+A/j
            S<-(S%x%Id)%*%T
            if(j<r){A<-Id%x%A}
            }}
    return(S)
}



##########################################################################
### Symmetrizer matrix applied to a vector
##########################################################################
    

############################################################################
### Sdrv.direct computes the result of multiplying the symmetrizer matrix
### S_{d,r} by a vector v of length d^r, based on Equation (4) as described
### in Section 4 of Chacon and Duong (2014)
############################################################################ 

Sdrv.direct<-function(d,r,v){
    Sv<-rep(0,d^r)
    per<-permn.ks(r)
    per.rep<-pinv.all(d,r)
    nper<-factorial(r)
    nper.rep<-d^r
    per<-matrix(unlist(per), byrow=TRUE, ncol=r, nrow=nper)    
    pow<-0:(r-1)
    dpow<-d^pow
    
    if(nper.rep<=nper){
    dpow.mat<-matrix(rep(dpow,nper),byrow=TRUE,ncol=r,nrow=nper)
    for(i in 1:nper.rep){      ## Loop over no. perms with reps (d^r)
        pinvi<-per.rep[i,]
        sigpinvi<-matrix(pinvi[per],byrow=FALSE,nrow=nrow(per),ncol=ncol(per))
        psigpinvi<-drop(1+rowSums((sigpinvi-1)*dpow.mat))
        Sv[i]<-sum(tabulate(psigpinvi,nbins=d^r)*v)
    }}
    
    if(nper<nper.rep){
    dpow.mat<-matrix(rep(dpow,nper.rep),byrow=TRUE,ncol=r,nrow=nper.rep)
    for(s in 1:nper){          ## Loop over no. perms (r!)
        sig<-per[s,]
        sigpinv<-per.rep[,sig]
        psigpinv<-drop(1+rowSums((sigpinv-1)*dpow.mat))
        Sv<-Sv+v[psigpinv]
    }}    
    return(Sv/nper)
}
   

############################################################################
### Sdrv.recursive computes the result of multiplying the symmetrizer matrix
### S_{d,r} by a vector v of length d^r, based on the recursive Algorithm 2
### as described in Section 4 of Chacon and Duong (2014)
############################################################################ 

Sdrv.recursive<-function(d,r,v){
    if((!is.vector(v))&(!is.matrix(v))){stop("v must be a vector or a matrix")}
    if(is.vector(v)){v<-matrix(v,nrow=1)}
    n<-nrow(v)   
    if(ncol(v)!=d^r){stop("Length of the vector(s) must equal d^r")}
    
    if((r==0)|(r==1)){w<-v}
    
    else{
    per.rep<-pinv.all(d,r)
    nper.rep<-d^r
    dpow<-d^((r-1):0)
    dpow.mat<-matrix(rep(dpow,d^r),ncol=r,nrow=d^r,byrow=TRUE) 

    w<-v
    for(p in 2:r){
        Tv<-matrix(0,nrow=n,ncol=d^r)
        for(j in 1:p){
        per.rep.tau<-per.rep;per.rep.tau[,j]<-per.rep[,p];per.rep.tau[,p]<-per.rep[,j]
        positions<-1+rowSums((per.rep.tau-1)*dpow.mat)
        Tv[,positions]<-Tv[,positions]+w
    }
    w<-Tv/p}
    }
    return(drop(w))
}




#############################################################################
## Lp norm between two functions (grid based)
#############################################################################

Lpdiff <- function(f1, f2, p=2, index=1)
{
  if (is.vector(f1$H)) d <- 1 else d <- ncol(f1$H)
  f.diff <- f1
  if (is.list(f1$estimate)) f.diff$estimate <- (abs(f1$estimate[[index]] - f2$estimate[[index]]))^p
  else f.diff$estimate <- (abs(f1$estimate - f2$estimate))^p
  if (d==1)
  {
    delta <- diff(f.diff$eval.points)
    riemann.sum <- sum(c(delta[1], delta)*f.diff$estimate)
  }
  else if (d>1)
  {
    delta <- sapply(f.diff$eval.points, diff)
    delta <- rbind(head(delta, n=1), delta)
    if (d==2) riemann.sum <- sum(outer(delta[,1], delta[,2]) * f.diff$estimate)
    else if (d==3) riemann.sum <- sum(outer(outer(delta[,1], delta[,2]), delta[,3]) * f.diff$estimate)
  }  
  
  return(riemann.sum)
}
