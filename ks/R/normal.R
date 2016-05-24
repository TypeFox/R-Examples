
###############################################################################
## Univariate mixture normal densities
###############################################################################

rnorm.mixt <- function(n=100, mus=0, sigmas=1, props=1, mixt.label=FALSE)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one")

  ### single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
  {
    if (mixt.label)
      rand <- cbind(rnorm(n=n, mean=mus, sd=sigmas), rep(1, n))
    else
      rand <- rnorm(n=n, mean=mus, sd=sigmas)
  }
  ### multiple component mixture
  else
  {
    k <- length(props)
    n.samp <- sample(1:k, n, replace=TRUE, prob=props) 
    n.prop <- numeric(0)

    ## compute number taken from each mixture
    for (i in 1:k)
      n.prop <- c(n.prop, sum(n.samp == i))
    
    rand <- numeric(0)
    
    for (i in 1:k) ##for (i in as.numeric(rownames(n.prop)))
    {
      ## compute random sample from normal mixture component
      if (n.prop[i] > 0)
        if (mixt.label)
          rand <- rbind(rand, cbind(rnorm(n=n.prop[i], mean=mus[i], sd=sigmas[i]), rep(i, n.prop[i])))
        else
          rand <- c(rand, rnorm(n=n.prop[i], mean=mus[i], sd=sigmas[i]))
    }
  }
  if (mixt.label)
    return(rand[sample(n),])
  else
    return(rand[sample(n)])
}


dnorm.mixt <- function(x, mus=0, sigmas=1, props=1)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one")

  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
    dens <- dnorm(x, mean=mus[1], sd=sigmas[1])

  ## multiple component mixture
  else   
  {   
    k <- length(props)
    dens <- 0

    ## sum of each normal density value from each component at x  
    for (i in 1:k)
      dens <- dens + props[i]*dnorm(x, mean=mus[i], sd=sigmas[i])
  }
  
  return(dens)
}   


###############################################################################
## Partial derivatives of the univariate normal (mean 0) 
## 
## Parameters
## x - points to evaluate at
## sigma - std deviation
## r - derivative index 
#
## Returns
## r-th derivative at x
###############################################################################

dnorm.deriv <- function(x, mu=0, sigma=1, deriv.order=0)
{
  r <- deriv.order
  phi <- dnorm(x, mean=mu, sd=sigma) 
  x <- (x - mu)
  
  if (r==0)
    return(phi)
  else if (r==1)
    derivt <- -x/sigma^2*phi
  else if (r==2)
    derivt <- (x^2-sigma^2)/sigma^4*phi
  else if (r==3)
    derivt <- -(x^3 - 3*x*sigma^2)/sigma^6*phi
  else if (r==4)
    derivt <- (x^4 - 6*x^2*sigma^2 + 3*sigma^4)/sigma^8*phi
  else if (r==5)
    derivt <- -(x^5 - 10*x^3*sigma^2 + 15*x*sigma^4)/sigma^10*phi
  else if (r==6)
    derivt <- (x^6 - 15*x^4*sigma^2 + 45*x^2*sigma^4 - 15*sigma^6)/sigma^12*phi
  else if (r==7)
    derivt <- -(x^7 - 21*x^5*sigma^2 + 105*x^3*sigma^4 - 105*x*sigma^6)/sigma^14*phi
  else if (r==8)
    derivt <- (x^8 - 28*x^6*sigma^2 + 210*x^4*sigma^4 - 420*x^2*sigma^6 + 105*sigma^8)/sigma^16*phi
  else if (r==9)
    derivt <- -(x^9 - 36*x^7*sigma^2 + 378*x^5*sigma^4 - 1260*x^3*sigma^6 + 945*x*sigma^8)/sigma^18*phi
  else if (r==10)
    derivt <- (x^10 - 45*x^8*sigma^2 + 630*x^6*sigma^4 - 3150*x^4*sigma^6 + 4725*x^2*sigma^8 - 945*sigma^10)/sigma^20*phi
  
  if (r > 10)
    stop ("Up to 10th order derivatives only")
    
  return(derivt)
}

###############################################################################
## Double sum  of K(X_i - X_j) used in density derivative estimation
#
## Parameters
## x - points to evaluate
## Sigma - variance matrix
## inc - 0 - exclude diagonals
##     - 1 - include diagonals
#
## Returns
## Double sum at x
###############################################################################


dnorm.deriv.sum <- function(x, sigma, deriv.order, inc=1, binned=FALSE, bin.par, kfe=FALSE)
{
  r <- deriv.order
  n <- length(x)
  if (binned)
  {
    if (missing(bin.par)) bin.par <- binning(x, h=sigma, supp=4+r) 
    est <- kdde.binned(x=x, H=sigma^2, h=sigma, deriv.order=r, bin.par=bin.par)$estimate
    sumval <- sum(bin.par$counts*est*n)
    if (inc == 0) 
      sumval <- sumval - n*dnorm.deriv(x=0, mu=0, sigma=sigma, deriv.order=r)
  }
  else
  {
    sumval <- 0
    for (i in 1:n)
      sumval <- sumval + sum(dnorm.deriv(x=x[i] - x, mu=0, sigma=sigma, deriv.order=r)) 
    if (inc == 0) 
      sumval <- sumval - n*dnorm.deriv(x=0, mu=0, sigma=sigma, deriv.order=r)
  }

  if (kfe)
    if (inc==1) sumval <- sumval/n^2
    else sumval <- sumval/(n*(n-1))
  
  return(sumval)
  
}

dnorm.deriv.mixt <- function(x, mus=0, sigmas=1, props=1, deriv.order=0)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one")

  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
    dens <- dnorm.deriv(x, mu=mus[1], sigma=sigmas[1], deriv.order=deriv.order)

  ## multiple component mixture
  else   
  {   
    k <- length(props)
    dens <- 0

    ## sum of each normal density value from each component at x  
    for (i in 1:k)
      dens <- dens + props[i]*dnorm.deriv(x=x, mu=mus[i], sigma=sigmas[i], deriv.order=deriv.order)
  }
  
  return(dens)
}   

###############################################################################
# Multivariate normal densities and derivatives
###############################################################################


###############################################################################
## Multivariate normal mixture - random sample
## 
## Parameters
## n - number of samples
## mus - matrix of means (each row is a vector of means from each component
##       density)
## Sigmas - matrix of covariance matrices (every d rows is a covariance matrix 
##          from each component density) 
## props - vector of mixing proportions 
## 
## Returns
## Vector of n observations from the normal mixture 
###############################################################################

rmvnorm.mixt <- function(n=100, mus=c(0,0), Sigmas=diag(2), props=1, mixt.label=FALSE)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one")
  
  #if (is.vector(Sigmas))
  ##  return(rnorm.mixt(n=n, mus=mus, sigmas=Sigmas, props=props))
  
  ### single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
   if (mixt.label)
     rand <- cbind(rmvnorm(n=n, mean=mus, sigma=Sigmas), rep(1, n))
   else
     rand <- cbind(rmvnorm(n=n, mean=mus, sigma=Sigmas))
    
  ### multiple component mixture
  else
  {
    k <- length(props)
    d <- ncol(Sigmas)
    n.samp <- sample(1:k, n, replace=TRUE, prob=props) 
    n.prop <- numeric(0)

    ## compute number taken from each mixture
    for (i in 1:k)
      n.prop <- c(n.prop, sum(n.samp == i))
    
    rand <- numeric(0)
    
    for (i in 1:k)
    {
      ## compute random sample from normal mixture component
      if (n.prop[i] > 0)
      {       
        if (mixt.label)
          rand <- rbind(rand, cbind(rmvnorm(n=n.prop[i], mean=mus[i,], sigma=Sigmas[((i-1)*d+1) : (i*d),]), rep(i, n.prop[i])))
        else
          rand <- rbind(rand, rmvnorm(n=n.prop[i], mean=mus[i,], sigma=Sigmas[((i-1)*d+1) : (i*d),]))
      }    
    }
  }

  return(rand[sample(n),])
}


###############################################################################
## Multivariate normal mixture - density values
## 
## Parameters
## x - points to compute density at 
## mus - matrix of means
## Sigmas - matrix of covariance matrices 
## props - vector of mixing proportions 
## 
## Returns
## Density values from the normal mixture (at x)
###############################################################################

dmvnorm.mixt <- function(x, mus, Sigmas, props=1)
{  
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one")

  if (is.vector(x)) d <- length(x)
  else d <- ncol(x)
  
  if (missing(mus)) mus <- rep(0,d)
  if (missing(Sigmas)) Sigmas <- diag(d)
     
  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
  {
    if (is.matrix(mus)) mus <- mus[1,]
    dens <- dmvnorm(x=x, mean=mus, sigma=Sigmas[1:d,])
  }
  ## multiple component mixture
  else   
  {   
    k <- length(props)
    dens <- 0
    ## sum of each normal density value from each component at x  
    for (i in 1:k)
      dens <- dens + props[i]*dmvnorm(x=x, mean=mus[i,], sigma=Sigmas[((i-1)*d+1):(i*d),])
  }
  
  return(dens)
}   

##########################################################################
### Computation of the r-th derivative vector of the Gaussian density
##########################################################################

dmvnorm.deriv <- function(x, mu, Sigma, deriv.order=0, deriv.vec=TRUE, add.index=FALSE, only.index=FALSE, type="unique")
{
  type1 <- match.arg(type, c("recursive", "direct", "unique"))

  r <- deriv.order
  if(length(r)>1) stop("deriv.order should be a non-negative integer")
  sumr <- sum(r)
  
  if (missing(x)) d <- ncol(Sigma)
  else
  {
    if (is.vector(x))
      x <- t(as.matrix(x))
    
    if (is.data.frame(x)) x <- as.matrix(x)
    d <- ncol(x)
    n <- nrow(x)
  }
  if (add.index | only.index | !deriv.vec)
  {    
      ## matrix of derivative indices
      ind.mat <- 0
      sumr.counter <- sumr
      if (sumr>=1) ind.mat <- diag(d)
      {
          while (sumr.counter >1)
          {
              ind.mat <- Ksum(diag(d), ind.mat)
              sumr.counter <- sumr.counter - 1
          }
      }
      ind.mat.minimal <- unique(ind.mat)
      ind.mat.minimal.logical <- !duplicated(ind.mat)
      
      if (only.index)
          if (deriv.vec) return (ind.mat)
          else  return(ind.mat.minimal)
  }
  
  if (missing(mu)) mu <- rep(0,d)
  if (missing(Sigma)) Sigma <- diag(d)
 
  x.centred <- sweep(x, 2, mu)
 
  dens <- do.call(paste("dmvnorm.deriv", type1, sep="."), list(x=x.centred, Sigma=Sigma, deriv.order=r))
  if (is.vector(dens) & r>0) dens <- matrix(dens, nrow=1)
  
  if (!deriv.vec & r>0)
  {
    ind.select <- numeric()
    for (i in 1:nrow(ind.mat.minimal))
      ind.select <- c(ind.select, head(which.mat(ind.mat.minimal[i,], ind.mat), n=1))
    
    dens <- dens[,ind.select]
    ind.mat <- ind.mat.minimal
  }
  if (add.index) return(list(deriv=dens, deriv.ind=ind.mat))
  else return(dens)
}


############################################################################
### dmvnorm.deriv.direct computes the vector derivative of the Gaussian
### density phi_Sigma(x) on the basis of Equation (1) and Algotihm 2, as
### described in Chacon and Duong (2014)
############################################################################

dmvnorm.deriv.direct<-function(x,Sigma,deriv.order=0){ 
  if(is.vector(x)){x<-matrix(x,nrow=1)}
  d<-ncol(Sigma)
  n<-nrow(x)
  r<-deriv.order
  
  Sigmainv<-chol2inv(chol(Sigma))
  Hermx<-matrix(0,nrow=n,ncol=d^r)
  for(i in 1:n)
  {
    for (j in 0:floor(r/2))
      Hermx[i,]<-Hermx[i,] + (-1)^j/(factorial(j)*factorial(r-2*j)*2^j)*
      Kpow(Sigmainv%*%x[i,], r-2*j)%x%Kpow(vec(Sigmainv),j)
    Hermx[i,]<-Sdrv.recursive(d=d,r=r,v=Hermx[i,])
  }
  dens<-(-1)^r*factorial(r)*Hermx*dmvnorm(x,mean=rep(0,d), sigma=Sigma)
  return(drop(dens))
}

############################################################################
### dmvnorm.deriv.recursive computes the vector derivative of the Gaussian
### density phi_Sigma(x) on the basis of Equation (7) and Algorithm 2 as
### described in Section 5 of Chacon and Duong (2014)
############################################################################

dmvnorm.deriv.recursive<-function(x,Sigma,deriv.order=0){ 
    if(is.vector(x)){x<-matrix(x,nrow=1)}
    d<-ncol(Sigma)
    n<-nrow(x)
    r<-deriv.order
    G<-Sigma
    Ginv<-chol2inv(chol(G))
    vGinv<-vec(Ginv)
    nvGinv<-matrix(rep(vGinv,n),byrow=TRUE,ncol=length(vGinv),nrow=n)
    arg<-matrix(x%*%Ginv, nrow=n)
    
    hmold0 <- matrix(1,nrow=n,ncol=1)
    hmold1 <- arg
    hmnew <- hmold0

    if(r==1){hmnew<-hmold1}
    if(r>=2){for(i in 2:r){
        hmnew<-mat.Kprod(hmold1,arg)-(i-1)*mat.Kprod(nvGinv,hmold0)
        hmnew<-matrix(Sdrv.recursive(d=d,r=i,v=hmnew), nrow=n)
        hmold0<-hmold1
        hmold1<-hmnew        
        }    
    }
    dens<-dmvnorm(x,mean=rep(0,d),sigma=Sigma)
    result<-matrix(rep(dens,d^r),byrow=FALSE,nrow=n,ncol=d^r)*hmnew*(-1)^r
    return(drop(result))
}

###############################################################################
### dmvnorm.deriv.unique computes the whole vector derivative of the Gaussian
### density phi_Sigma(x) from its unique coordinates, based on Algorithm 3 as
### described in Section 5 of Chacon and Duong (2014)
###############################################################################

dmvnorm.deriv.unique<-function(x,Sigma,deriv.order=0){
    if(is.vector(x)){x<-matrix(x,nrow=1)}
    d<-ncol(x)
    n<-nrow(x)
    r<-deriv.order 
    G<-Sigma
    Ginv<-chol2inv(chol(G))
    arg<-x%*%Ginv
    
    hmold0 <- matrix(1,nrow=n,ncol=1)
    hmold1 <- arg
    hmnew <- hmold0  
    udind0<-matrix(rep(0,d),nrow=1,ncol=d)
    udind1<-diag(d)
     
    if(r==1){hmnew<-hmold1}
    if(r>=2){for(i in 2:r){
        Ndi1<-ncol(hmold1)
        Ndi0<-ncol(hmold0)
        hmnew<-numeric()
        
        for(j in 1:d){        
        nrecj<-choose(d-j+i-1,i-1)        
        hmnew.aux<-arg[,j]*hmold1[,Ndi1-(nrecj:1)+1]
              
        for(k in j:d){        
            udind0.aux<-matrix(udind1[Ndi1-(nrecj:1)+1,],ncol=d,byrow=FALSE)
            udind0.aux[,k]<-udind0.aux[,k]-1
            valid.udind0<-as.logical(apply(udind0.aux>=0,1,min))
            enlarged.hmold0<-matrix(0,ncol=nrow(udind0.aux),nrow=n)
            for(l in 1:nrow(udind0.aux)){
                if(valid.udind0[l]){
                    pos<-which(rowSums((udind0-matrix(rep(udind0.aux[l,],nrow(udind0)),
                         nrow=nrow(udind0),byrow=TRUE))^2)==0)
                    enlarged.hmold0[,l]<-hmold0[,pos]}
                }
                hmnew.aux<-hmnew.aux-Ginv[j,k]*matrix(rep(udind1[Ndi1-(nrecj:1)+1,k],n),
                           nrow=n,byrow=TRUE)*enlarged.hmold0
            }
            hmnew<-cbind(hmnew,hmnew.aux)
        }
        hmold0 <- hmold1
        hmold1 <- hmnew

        ##Compute the unique i-th derivative multi-indexes 
        nudind1<-nrow(udind1)
        udindnew<-numeric()
        for(j in 1:d){
            Ndj1i<-choose(d+i-1-j,i-1)
            udind.aux<-matrix(udind1[nudind1-(Ndj1i:1)+1,],ncol=d,byrow=FALSE)
            udind.aux[,j]<-udind.aux[,j]+1
            udindnew<-rbind(udindnew,udind.aux)
        }
        udind0<-udind1
        udind1<-udindnew
        }}
        
    if(r==0) result<-dmvnorm(x,mean=rep(0,d),sigma=Sigma)
    if(r==1) result<-(-1)*matrix(rep(dmvnorm(x,mean=rep(0,d),sigma=Sigma),d),
                     nrow=n,byrow=FALSE)*hmnew
    if(r>=2){
        per<-pinv.all(d=d,r=r)
        dind<-numeric()
        udind<-udind1
        dind.base<-rep(0,d^r)        
        udind.base<-rep(0,choose(d+r-1,r))
                
        for(i in 1:d){
            dind<-cbind(dind,rowSums(per==i)) ## Matrix of derivative indices
            dind.base<-dind.base+dind[,i]*(r+1)^(d-i) ## Transform each row to base r+1           
            udind.base<-udind.base+udind[,i]*(r+1)^(d-i) ## Transform each row to base r+1                       
        }        
        
        dlabs<-match(dind.base,udind.base)
        
        deriv.vector<-hmnew[,dlabs]*matrix(rep((-1)^rowSums(dind),n),nrow=n,byrow=TRUE)
       
        result<-matrix(rep(dmvnorm(x,mean=rep(0,d),sigma=Sigma),ncol(deriv.vector)),
                nrow=n,byrow=FALSE)*deriv.vector
    }
    return(drop(result))
}


dmvnorm.deriv.mixt <- function(x, mus, Sigmas, props, deriv.order, deriv.vec=TRUE, add.index=FALSE, only.index=FALSE)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one")

  if (is.vector(x)) d <- length(x)
  else d <- ncol(x)
  
  if (missing(mus)) mus <- rep(0,d)
  if (missing(Sigmas)) Sigmas <- diag(d)
  r <- deriv.order
  sumr <- sum(r)
 
  if (only.index | add.index) 
    ind.mat <- dmvnorm.deriv(x=x, mu=mus[1,], Sigma=Sigmas[1:d,], deriv.order=r, only.index=TRUE)
  if (only.index)  
    if (deriv.vec) return (ind.mat)
    else return(unique(ind.mat))
  
  ## derivatives 
  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
  {
    if (is.matrix(mus)) mus <- mus[1,]
    dens <- dmvnorm.deriv(x=x, mu=mus, Sigma=Sigmas[1:d,], deriv.order=sumr)
  }
  ## multiple component mixture
  else   
  {   
    k <- length(props)
    dens <- 0
    ## sum of each normal density value from each component at x  
    for (i in 1:k)
      dens <- dens + props[i]*dmvnorm.deriv(x=x, mu=mus[i,], Sigma=Sigmas[((i-1)*d+1):(i*d),], deriv.order=sumr)  
  }

  if (!deriv.vec)
  { 
    dens <- dens[,!duplicated(ind.mat)]
    ind.mat <- unique(ind.mat)
  }
 
  if (add.index) return(list(deriv=dens, deriv.ind=ind.mat))
  else return(deriv=dens)
}



###############################################################################
## Double sum  of K(X_i - X_j) used in density derivative estimation
#
## Parameters
## x - points to evaluate
## Sigma - variance matrix
## inc - 0 - exclude diagonals
##     - 1 - include diagonals
#
## Returns
## Double sum at x
##############################################################################



dmvnorm.deriv.sum <- function(x, Sigma, deriv.order=0, inc=1, binned=FALSE, bin.par, bgridsize, kfe=FALSE, deriv.vec=TRUE, add.index=FALSE, verbose=FALSE)
{
  r <- deriv.order
  d <- ncol(x)
  n <- nrow(x)
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
 
  if (binned)
  {
    d <- ncol(Sigma)
    n <- nrow(x)
    ##if (is.diagonal(Sigma))
    ##{
      if (missing(bin.par)) bin.par <- binning(x, H=diag(diag(Sigma)), bgridsize=bgridsize)  
      est <- kdde.binned(x=x, bin.par=bin.par, H=Sigma, deriv.order=r, verbose=verbose)$estimate
      
      if (r>0)
      {
        sumval <- rep(0, length(est))
        for (j in 1:length(est)) sumval[j] <- sum(bin.par$counts * n * est[[j]])
      }
      else
        sumval <- sum(bin.par$counts * n * est)
    ##}
    ## transformation approach from Jose E. Chacon 06/12/2010
    if (0)
    {
      Sigmainv12 <- matrix.sqrt(chol2inv(chol(Sigma)))
      y <- x %*% Sigmainv12
      if (missing(bin.par)) bin.par <- binning(x=y, H=diag(d), bgridsize=bgridsize)
     
      est <- kdde.binned(x=y, bin.par=bin.par, H=diag(d), deriv.order=r, verbose=verbose)$estimate
      if (r>0)
      {
        sumval <- rep(0, length(est))
        for (j in 1:length(est)) sumval[j] <- sum(bin.par$counts *n*est[[j]]) 
      }
      else
        sumval <- sum(bin.par$counts * n * est)

      sumval <- det(Sigmainv12) * sumval  %*% Kpow(Sigmainv12, pow=r)
    }
  }
  ## exact computation 
  else
  {
    if (verbose) pb <- txtProgressBar() 
    if (r==0)
    {
      n.seq <- block.indices(n, n, d=d, r=r, diff=TRUE)
      sumval <- 0
      for (i in 1:(length(n.seq)-1))
      {  
        difs <- differences(x=x, y=x[n.seq[i]:(n.seq[i+1]-1),])
        sumval <- sumval + sum(dmvnorm.deriv(x=difs, mu=rep(0,d), Sigma=Sigma, deriv.order=r, deriv.vec=deriv.vec))
        if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1)) 
      }
    }
    else
    {
      ## only with r>0 
      ## original recursive code from Jose E. Chacon 03/2012
      if (2*floor(r/2)!=r){sumval <- rep(0,d^r)}
      else
      {
        Sigmainv<-chol2inv(chol(Sigma))
        per <- perm.rep(d=d,r=r)
        dind <- numeric()
        for(i in 1:d) dind <- cbind(dind,rowSums(per==i)) ###Matrix of derivative indices
        
        udind <- unique(dind)
        nudind <- nrow(udind)        
        dlabs <- numeric(nrow(dind))
        for (i in 1:nrow(udind))
          dlabs <- dlabs+i*(rowSums((dind-matrix(rep(udind[i,],nrow(dind)),nrow=nrow(dind),byrow=TRUE))^2)==0)
              
        result<-rep(0,nudind)
        
        ndif <- n*(n-1)/2
        dif.ind <- numeric()
        ##for(k in 2:n) dif.ind <- rbind(dif.ind,cbind(1:(k-1),rep(k,k-1)))
        
        M <- 1e6
        max.loop.size <- ceiling(M/nudind) ### Inside the loop we need to store a matrix of order max.loop.size x nudind <= M
        nblocks <- ceiling(ndif/max.loop.size)
        blength <- c(rep(max.loop.size,nblocks-1),ndif-max.loop.size*(nblocks-1)) #length of each of the blocks, the last one could be smaller
        
        tri.num <- (1:n)*((1:n)-1)/2
        for(kk in 1:nblocks){
          b <- blength[kk]
          if (verbose) setTxtProgressBar(pb, kk/nblocks) 
          kkm <- (kk-1)*max.loop.size
         
          tri.ind <- findInterval(kkm:(kkm+b-1), tri.num)
          dif.ind.block <- cbind(kkm:(kkm+b-1) - tri.num[tri.ind]+1, tri.ind+1)
                             
          ##difs.block <- x[dif.ind[((kk-1)*max.loop.size+1):((kk-1)*max.loop.size+b),1],]-x[dif.ind[((kk-1)*max.loop.size+1):((kk-1)*max.loop.size+b),2],]
          difs.block <- x[dif.ind.block[,1],]-x[dif.ind.block[,2],]
          
          arg <- difs.block %*% Sigmainv    
          narg <- nrow(arg) 
          
          hmold0 <- matrix(rep(1,narg),ncol=1,nrow=narg)
          hmold1 <- arg
          hmnew <- hmold0
          
          udind0 <- matrix(rep(0,d),nrow=1,ncol=d)
          udind1 <- diag(d)
          
          if (r==1){hmnew<-hmold1}
          if (r >= 2)
          { 
            for (i in 2:r)
            {     	     
              Ndi1 <- ncol(hmold1)
              ##Ndi0 <- ncol(hmold0)
              hmnew <- numeric()
              
              for(j in 1:d)
              {        
                nrecj <- choose(d-j+i-1,i-1)        
                hmnew.aux <- arg[,j]*hmold1[,Ndi1-(nrecj:1)+1]
                
                for(k in j:d)
                {        
                  udind0.aux <- matrix(udind1[Ndi1-(nrecj:1)+1,],ncol=d,byrow=FALSE)
                  udind0.aux[,k] <- udind0.aux[,k]-1
                  valid.udind0 <- as.logical(apply(udind0.aux>=0,1,min))
                  enlarged.hmold0 <- matrix(0,ncol=nrow(udind0.aux),nrow=narg)
                  for(ell in 1:nrow(udind0.aux))
                  {
                    if(valid.udind0[ell]){
                      pos <- which(rowSums((udind0-matrix(rep(udind0.aux[ell,],nrow(udind0)),nrow=nrow(udind0),byrow=TRUE))^2)==0)
                      enlarged.hmold0[,ell]<-hmold0[,pos]}
                  }
                  ##In enlarged.hmold0 we put the vector hmold0 in those positions not having a -1 in any of the derivative order after subtracting e_k
                  ##The remaining positions are zeroes
                  ##Surely this could be done in a more efficient way
                  
                  hmnew.aux <- hmnew.aux-Sigmainv[j,k]*matrix(rep(udind1[Ndi1-(nrecj:1)+1,k],narg),nrow=narg,byrow=TRUE)*enlarged.hmold0
                }
            		hmnew<-cbind(hmnew,hmnew.aux)
              }
              hmold0 <- hmold1
              hmold1 <- hmnew
              
              ##Compute the unique i-th derivative multi-indexes 
              nudind1<-nrow(udind1)
              udindnew<-numeric()
              
              for(j in 1:d)
              {
                Ndj1i <- choose(d+i-1-j,i-1)
                udind.aux <- matrix(udind1[nudind1-(Ndj1i:1)+1,],ncol=d,byrow=FALSE)
                udind.aux[,j] <- udind.aux[,j]+1
                udindnew <- rbind(udindnew,udind.aux)
              }
              udind0 <- udind1
              udind1 <- udindnew
            }
          }
          hmnew<-hmnew*matrix(rep((-1)^rowSums(udind),narg),nrow=narg,byrow=TRUE)
          phi <- dmvnorm(difs.block, mean=rep(0,d), sigma=Sigma)
          phi <- matrix(rep(phi,nudind),ncol=nudind,byrow=FALSE)
          result <- result+drop(colSums(phi*hmnew))
        }
        
        result <- result[dlabs]    
        hm0 <- dmvnorm.deriv(x=rep(0,d),Sigma=Sigma,deriv.order=deriv.order)
        sumval <- 2*result+n*hm0
        if (verbose) close(pb)
        }   
    }        
    if (verbose) close(pb)
  }
  
  if (inc==0)
    sumval <- sumval - n*dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=Sigma, deriv.order=r)

  sumval <- drop(sumval)  
  if (kfe) { if (inc==1) sumval <- sumval/n^2 else sumval <- sumval/(n*(n-1)) }
  
  if (add.index)
  {
    ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=r, only.index=TRUE)
    if (deriv.vec) return(list(sum=sumval, deriv.ind=ind.mat))
    else return(list(sum=sumval, deriv.ind=unique(ind.mat)))
  }
  else return(sum=sumval)
}

## Single partial derivative of the multivariate normal with scalar variance matrix sigma^2 I_d  
## Code by Jose  Chacon 04/09/2007


dmvnorm.deriv.scalar <- function(x, mu, sigma, deriv.order, binned=FALSE)
{
  r <- deriv.order
  d <- ncol(x)

  sderiv <- sum(r)
  arg <- x/sigma
  darg <- dmvnorm(arg, mean=mu)/(sigma^(sderiv+d))
  for (j in 1:d)
  {
    hmold0 <- 1
    hmold1 <- arg[,j]
    hmnew <- 1
    if (r[j] ==1){hmnew<-hmold1}
    if (r[j] >= 2) ## Multiply by the corresponding Hermite polynomial, coordinate-wise, using Fact C.1.4 in W&J (1995) and Willink (2005, p.273)
      for (i in (2:r[j]))
      {
        hmnew <- arg[,j] * hmold1 - (i - 1) * hmold0
        hmold0 <- hmold1
        hmold1 <- hmnew
      }
    darg <- hmnew * darg
  }
  
  val <- darg*(-1)^sderiv
  return(val)
}


dmvnorm.deriv.scalar.sum <- function(x, sigma, deriv.order=0, inc=1, kfe=FALSE, binned=FALSE, bin.par, verbose=FALSE)
{
  r <- deriv.order
  d <- ncol(x)
  n <- nrow(x)

  if (binned)
  {
    if (missing(bin.par)) bin.par <- binning(x, H=diag(d)*sigma^2)  
    n <- sum(bin.par$counts)

    ind.mat <- dmvnorm.deriv(x=rep(0,d), Sigma=diag(d), deriv.order=sum(r), deriv.vec=TRUE, only.index=TRUE)
    fhatr <- kdde.binned(bin.par=bin.par, H=sigma^2*diag(d), deriv.order=sum(r), deriv.vec=TRUE, w=rep(1,n), deriv.index=which.mat(r=r, ind.mat)[1])
    sumval <- sum(bin.par$counts * n * fhatr$est[[1]])
  }
  else
  {
    if (verbose) pb <- txtProgressBar() 
    
    n.seq <- block.indices(n, n, d=d, r=r, diff=TRUE)
    sumval <- 0
    for (i in 1:(length(n.seq)-1))
    {  
      difs <- differences(x=x, y=x[n.seq[i]:(n.seq[i+1]-1),])
      sumval <- sumval + sum(dmvnorm.deriv.scalar(x=difs, mu=rep(0,d), sigma=sigma, deriv.order=r))
      if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1)) 
    }
  }
  if (verbose) close(pb)
  
  if (inc==0)
    sumval <- sumval - n*dmvnorm.deriv.scalar(x=t(as.matrix(rep(0,d))), mu=rep(0,d), sigma=sigma, deriv.order=r)
  
  if (kfe)
    if (inc==1) sumval <- sumval/n^2
    else sumval <- sumval/(n*(n-1))
  
  return(sumval)
}

##########################################################################
## Normal scale psi functionals
##########################################################################

psins.1d <- function(r, sigma)
{
  if (r %% 2 ==0)
    psins <- (-1)^(r/2)*factorial(r)/((2*sigma)^(r+1)*factorial(r/2)*pi^(1/2))
  else
    psins <- 0
    
  return(psins)  
}


psins <- function(r, Sigma, deriv.vec=length(r)==1)
{
  d <- ncol(Sigma)
  if (deriv.vec)
  {
    dens <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), deriv.order=r, Sigma=2*Sigma, add.index=FALSE)
    return(drop(dens)) ##dens$deriv
  }
  else
  {
     dens <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), deriv.order=sum(r), Sigma=2*Sigma, add.index=TRUE)
     if (!is.vector(dens$deriv.ind))
     {
       i <- head(which.mat(r, dens$deriv.ind),n=1)
       dens <- dens$deriv[1,i]
     }
     else
       dens <- dens$deriv
    return(dens)
  }
}


##########################################################################
### Vector moments of the normal distribution
##########################################################################

mur <- function(r, A, mu, Sigma, type="unique")
{
  type1 <- match.arg(type, c("direct", "recursive", "unique"))
  mur.val <- do.call(paste("mur", type1, sep="."), list(r=r, A=A, mu=mu, Sigma=Sigma))
  return(mur.val)
}


#############################################################################
### mur.direct computes the vector moment E[X^{\otimes r}] for a random
### vector with N(mu,Sigma) distribution, on the basis of Equation (8) in
### Section 6 of Chacon and Duong (2014)
#############################################################################

mur.direct<-function(r,mu,Sigma){
    d<-ncol(Sigma)
    result<-as.vector(Kpow(mu,r))
    vS<-vec(Sigma)
    vSj<-1
    if(r>=2){
        for(j in 1:floor(r/2)){
            vSj<-as.vector(vSj%x%vS)
            cj<-prod(r:(r-2*j+1))/(prod(1:j)*2^j)
            mur2j<-as.vector(Kpow(mu,r-2*j))
            result<-result+cj*as.vector(mur2j%x%vSj)
            }
        }
    return(drop(Sdrv.recursive(d=d,r=r,v=result)))
}



#############################################################################
### mur.recursive computes the vector moment E[X^{\otimes r}] for a random
### vector with N(mu,Sigma) distribution, on the basis of Equation (9) in
### Section 6 of Chacon and Duong (2014), using Equation (7) in Section 5
### to obtain the Hermite polynomial
#############################################################################

mur.recursive<-function(r,mu,Sigma){ 
    d<-ncol(Sigma)
    G<- -Sigma
    vG<-vec(G)
    arg<-mu
    hmold0 <- 1
    hmold1 <- arg
    hmnew <- hmold0

    if(r==1){hmnew<-hmold1}
    if(r>=2){for(i in 2:r){
        hmnew<-as.vector(arg%x%hmold1-(i-1)*(vG%x%hmold0))
        hmold0<-hmold1
        hmold1<-hmnew        
        }    
    }
    return(drop(Sdrv.recursive(d=d,r=r,v=hmnew)))   
}

###############################################################################
### mur.unique computes the vector moment E[X^{\otimes r}] for a random vector 
### with N(mu,Sigma) distribution, on the basis of Equation (9) in Section 6
### of Chacon and Duong (2014), using Algorithm 3 in Section 5, based on the
### unique partial derivatives, to obtain the Hermite polynomial
###############################################################################

mur.unique<-function(r,mu,Sigma){
    d<-ncol(Sigma)
    G<- -Sigma
    arg<-mu
    hmold0 <- 1
    hmold1 <- arg
    hmnew <- hmold0
    udind0<-matrix(rep(0,d),nrow=1,ncol=d)
    udind1<-diag(d)
    
    if(r==1){hmnew<-hmold1}
    if(r>=2){for(i in 2:r){
        Ndi1<-length(hmold1)
        Ndi0<-length(hmold0)
        hmnew<-numeric()
        
        for(j in 1:d){        
        nrecj<-choose(d-j+i-1,i-1)        
        hmnew.aux<-arg[j]*hmold1[Ndi1-(nrecj:1)+1]
              
        for(k in j:d){        
            udind0.aux<-matrix(udind1[Ndi1-(nrecj:1)+1,],ncol=d,byrow=FALSE)
            udind0.aux[,k]<-udind0.aux[,k]-1
            valid.udind0<-as.logical(apply(udind0.aux>=0,1,min))
            enlarged.hmold0<-rep(0,nrow(udind0.aux))
            for(l in 1:nrow(udind0.aux)){
                if(valid.udind0[l]){
                    pos<-which(rowSums((udind0-matrix(rep(udind0.aux[l,],nrow(udind0)),
                         nrow=nrow(udind0),byrow=TRUE))^2)==0)
                    enlarged.hmold0[l]<-hmold0[pos]}
                }
            hmnew.aux<-hmnew.aux-G[j,k]*udind1[Ndi1-(nrecj:1)+1,k]*enlarged.hmold0
            }
        hmnew<-c(hmnew,hmnew.aux)
        }
        hmold0 <- hmold1
        hmold1 <- hmnew

        ## Compute the unique i-th derivative multi-indexes 
        nudind1<-nrow(udind1)
        udindnew<-numeric()
        for(j in 1:d){
            Ndj1i<-choose(d+i-1-j,i-1)
            udind.aux<-matrix(udind1[nudind1-(Ndj1i:1)+1,],ncol=d,byrow=FALSE)
            udind.aux[,j]<-udind.aux[,j]+1
            udindnew<-rbind(udindnew,udind.aux)
        }
        udind0<-udind1
        udind1<-udindnew
        }}
        
    if(r==0){ result<-1 }
    if(r==1){ result<-(-1)*hmnew } 
    if(r>=2){
        per<-pinv.all(d=d,r=r)
        
        dind<-numeric()
        udind<-udind1
        dind.base<-rep(0,d^r)        
        udind.base<-rep(0,choose(d+r-1,r))
                
        for(i in 1:d){
            dind<-cbind(dind,rowSums(per==i)) ## Matrix of derivative indices
            dind.base<-dind.base+dind[,i]*(r+1)^(d-i) ## Transform each row to base r+1           
            udind.base<-udind.base+udind[,i]*(r+1)^(d-i) ## Transform each row to base r+1                       
        }        
        
        dlabs<-match(dind.base,udind.base)
   
    result<-hmnew[dlabs]*(-1)^rowSums(dind)
    }
    return(drop(result*(-1)^r))
}


##########################################################################
### Moments of quadratic forms in normal variables
##########################################################################

nur <- function(r, A, mu, Sigma, type="cumulant")
{
  type1 <- match.arg(type, c("direct", "recursive", "unique", "cumulant"))
  nur.val <- do.call(paste("nur", type1, sep="."), list(r=r, A=A, mu=mu, Sigma=Sigma))
  return(nur.val)
}

nurs <- function(r, s, A, B, mu, Sigma, type="cumulant")
{
  type1 <- match.arg(type, c("direct", "recursive", "unique", "cumulant"))
  nur.val <- do.call(paste("nurs", type1, sep="."), list(r=r, s=s, A=A, B=B, mu=mu, Sigma=Sigma))
  return(nur.val)
}


#############################################################################
### nur.direct computes the moment E[(X^T AX)^r] of the quadratic form
### X^T AX where X is a random vector with N(mu,Sigma) distribution, using
### Equation (10) in Section 6 of Chacon and Duong (2014), and the direct
### implementation mur.direct of the normal moments
#############################################################################

nur.direct<-function(r,A,mu,Sigma){
    vA<-vec(A)
    result<-drop(Kpow(t(vA),r)%*%mur.direct(2*r,mu,Sigma))
    return(result)
}

#############################################################################
### nur.recursive computes the moment E[(X^T AX)^r] of the quadratic form
### X^T AX where X is a random vector with N(mu,Sigma) distribution, using
### Equation (10) in Section 6 of Chacon and Duong (2014), and the recursive
### implementation mur.recursive of the normal moments
#############################################################################

nur.recursive<-function(r,A,mu,Sigma){
    vA<-vec(A)
    result<-drop(Kpow(t(vA),r)%*%mur.recursive(2*r,mu,Sigma))
    return(result)
}

#############################################################################
### nur.unique computes the moment E[(X^T AX)^r] of the quadratic form
### X^T AX where X is a random vector with N(mu,Sigma) distribution, using
### Equation (10) in Section 6 of Chacon and Duong (2014), and the function
### mur.unique to compute the normal moments from its unique coordinates
#############################################################################

nur.unique<-function(r,A,mu,Sigma){
    vA<-vec(A)
    result<-sum(Kpow(vA,r)*mur.unique(2*r,mu,Sigma))
    return(result)
}

#############################################################################
### nur.cumulant computes the moment E[(X^T AX)^r] of the quadratic form
### X^T AX where X is a random vector with N(mu,Sigma) distribution, using
### the recursive formula relating moments and cumulants
#############################################################################

nur.cumulant<-function(r,A,mu,Sigma){
    if(r==0){result<-1}
    if(r==1){result<-sum(diag(A%*%Sigma))+t(mu)%*%A%*%mu}
    if(r>=2){
        ASigma<-A%*%Sigma
        AS<-ASigma
        Amu<-A%*%mu
        kappas<-sum(diag(ASigma)+mu*Amu)
        nus<-kappas
        for(k in 2:r){
            knew<-k*t(mu)%*%ASigma%*%Amu
            ASigma<-ASigma%*%AS
            knew<-(knew+sum(diag(ASigma)))*factorial(k-1)*2^(k-1)
            nnew<-knew+sum(choose(k-1,1:(k-1))*nus*rev(kappas))
            kappas<-c(kappas,knew)
            nus<-c(nus,nnew)
            }
        result<-nnew
        }
    return(drop(result))
}

#############################################################################
### nurs.direct computes the joint moment E[(X^T AX)^r (X^T BX)^s] of the
### quadratic forms X^T AX and X^T BX, where X is a random vector with
### N(mu,Sigma) distribution, using Equation (10) in Section 6 of Chacon and
### Duong (2014), and the direct implementation mur.direct of the
### normal moments
#############################################################################

nurs.direct<-function(r,s,A,B,mu,Sigma){
    vA<-vec(A)
    vB<-vec(B)    
    result<-(Kpow(t(vA),r)%x%Kpow(t(vB),s))%*%mur.direct(2*r+2*s,mu,Sigma)
    return(drop(result))
}

#############################################################################
### nurs.recursive computes the joint moment E[(X^T AX)^r (X^T BX)^s] of the
### quadratic forms X^T AX and X^T BX, where X is a random vector with
### N(mu,Sigma) distribution, using Equation (10) in Section 6 of Chacon and
### Duong (2014), and the recursive implementation mur.recursive of the
### normal moments
#############################################################################

nurs.recursive<-function(r,s,A,B,mu,Sigma){
    vA<-vec(A)
    vB<-vec(B)    
    result<-(Kpow(t(vA),r)%x%Kpow(t(vB),s))%*%mur.recursive(2*r+2*s,mu,Sigma)
    return(drop(result))
}

#############################################################################
### nurs.unique computes the joint moment E[(X^T AX)^r (X^T BX)^s] of the
### quadratic forms X^T AX and X^T BX, where X is a random vector with
### N(mu,Sigma) distribution, using Equation (10) in Section 6 of Chacon and
### Duong (2014), and the function mur.unique to compute the normal moments
### from its unique coordinates
#############################################################################

nurs.unique<-function(r,s,A,B,mu,Sigma){
    vA<-vec(A)
    vB<-vec(B)    
    result<-drop((Kpow(t(vA),r)%x%Kpow(t(vB),s))%*%mur.unique(2*r+2*s,mu,Sigma))
    return(drop(result))
}

#############################################################################
### nurs.cumulant computes the joint moment E[(X^T AX)^r (X^T BX)^s] of the
### quadratic forms X^T AX and X^T BX, where X is a random vector with
### N(mu,Sigma) distribution, using the recursive formula (11) in Section 6
### of Chacon and Duong (2014), relating moments and cumulants. The cumulants
### are computed using the function kappars, which is based on Theorem 3
#############################################################################

kappars<-function(r,s,A,B,mu,Sigma){ 
    d<-ncol(A)
    if(r+s>1 & r>0 & s>0){ind<-multicool::allPerm(multicool::initMC(c(rep(1,r),rep(2,s))))}
    if(r+s==1 | r==0 | s==0){ind<-matrix(c(rep(1,r),rep(2,s)),nrow=1)}
    if(r+s==0){return(0)}
    nper<-nrow(ind)
   
    result<-0
    Dmat<-solve(Sigma)%*%mu%*%t(mu)
    ASigma<-A%*%Sigma
    BSigma<-B%*%Sigma
    Id<-diag(d)
    
    for(i in 1:nper){
        product<-Id
        for(j in 1:(r+s)){
            if(ind[i,j]==1){product<-product%*%ASigma}
            else if(ind[i,j]==2){product<-product%*%BSigma}
          }
        result<-result+sum(diag(product%*%(Id/(r+s)+Dmat))) 
    }
    
    result<-result*factorial(r)*factorial(s)*2^(r+s-1)
    return(drop(result))
}


nurs.cumulant<-function(r,s,A,B,mu,Sigma){ 
    if(r==0 & s>0){nurs<-nur.cumulant(s,B,mu,Sigma)}
    if(r>0 & s==0){nurs<-nur.cumulant(r,A,mu,Sigma)}
    if(r==0 & s==0){nurs<-1}
    if((r>0)&(s>0)){
    
    K<-matrix(0,nrow=r+1,ncol=s)
    for(i in 0:r){for(j in 1:s){
        K[i+1,j]<-kappars(i,j,A,B,mu,Sigma)
        }}

    N<-matrix(0,nrow=r+1,ncol=s)
    for(i in 0:r){
        N[i+1,1]<-nur.cumulant(r=i,A=A,mu=mu,Sigma=Sigma)
        }
    if(s>1){
    for(j in 1:(s-1)){for(i in 0:r){
        Choose<-outer(choose(i,0:i),choose(j-1,0:(j-1)))
        N[i+1,j+1]<-sum(Choose*N[1:(i+1),1:j]*K[(i:0)+1,j:1])
    }}
    }
    Choose<-outer(choose(r,0:r),choose(s-1,0:(s-1)))  
    nurs<-sum(Choose*N[1:(r+1),1:s]*K[(r:0)+1,s:1])
    }
return(nurs)
}


##########################################################################
### V-statistics with multivariate Gaussian derivatives kernel
##########################################################################

Qr <- function(x, y, Sigma, deriv.order=0, inc=1, type="cumulant", verbose=FALSE)
{
  if (missing(y)) y <- x
  type1 <- match.arg(type, c("direct", "cumulant"))
  if (type1=="direct")
    Qr.val <- Qr.direct(x=x, y=y, Sigma=Sigma, r=deriv.order, inc=inc, verbose=verbose)
  else
    Qr.val <- Qr.cumulant(x=x, y=y, Sigma=Sigma, r=deriv.order, inc=inc, verbose=verbose)
  return(Qr.val)
}


############################################################################
### Qr.direct computes the V-statistic using the direct approach, as
### described in Section 6 of Chacon and Duong (2014)
############################################################################

Qr.direct <- function(x, y, Sigma, r=0, inc=1, binned=FALSE, bin.par, bgridsize, verbose=FALSE)
{
  if (is.vector(x)) x <- matrix(x, nrow=1)
  d <- ncol(x)

  eta <- drop(Kpow(vec(diag(d)), r/2) %*% kfe(x=x, G=Sigma, deriv.order=r, inc=inc, verbose=verbose, binned=binned, bin.par=bin.par, add.index=FALSE))
  ##if (inc==0) eta <- eta/(nrow(x)*(nrow(x)-1))
  ##if (inc==1) eta <- eta/(nrow(x)*nrow(x))
  return(eta)
}



############################################################################
### Qr.cumulant computes the V-statistic using the relationship with nur
### shown in Theorem 4 of Chacon and Duong (2014)
############################################################################

Qr.cumulant <- function(x, y, Sigma, r=0, inc=1, verbose=FALSE)
{
  if (is.vector(x)) x <- matrix(x, nrow=1)
  d <- ncol(x)
  r <- r/2
  if (missing(y)) y <- x
  if (is.vector(y)) y <- matrix(y, nrow=1)
  nx <- as.numeric(nrow(x))
  ny <- as.numeric(nrow(y))
  G <- Sigma
  Ginv <- chol2inv(chol(G))
  G2inv <- Ginv%*%Ginv
  G3inv <- G2inv%*%Ginv
  trGinv <- sum(diag(Ginv))
  trG2inv <- sum(diag(G2inv))   
  detG <- det(G)

  ## indices for separating into blocks for double sum calculation
  n.seq <- block.indices(nx, ny, d=d, r=r, diff=FALSE)
  if (verbose) pb <- txtProgressBar()
  
  if (r==0)
  {
    xG <- x%*%Ginv
    a <- rowSums(xG*x)
    eta <- 0
    for (i in 1:(length(n.seq)-1))
    {
      nytemp <- n.seq[i+1] - n.seq[i]
      ytemp <- matrix(y[n.seq[i]:(n.seq[i+1]-1),], ncol=d)
      aytemp <- rowSums((ytemp %*% Ginv) *ytemp)
      M <- a%*%t(rep(1,nytemp)) + rep(1, nx)%*%t(aytemp) - 2*(xG%*%t(ytemp))
      em2 <- exp(-M/2)
      eta <- eta + (2*pi)^(-d/2)*detG^(-1/2)*sum(em2)
      if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
    }
  } 
  else if (r==1)
  {
    xG <- x%*%Ginv
    xG2 <- x%*%G2inv
    a <- rowSums(xG*x)
    a2 <- rowSums(xG2*x)
    
    eta <- 0
    for (i in 1:(length(n.seq)-1))
    {
      nytemp <- n.seq[i+1] - n.seq[i]
      ytemp <- matrix(y[n.seq[i]:(n.seq[i+1]-1),], nrow=nytemp)
      aytemp <- rowSums((ytemp %*% Ginv) *ytemp)
      aytemp2 <- rowSums((ytemp %*% G2inv) *ytemp)
      M  <- a%*%t(rep(1,nytemp))+rep(1,nx)%*%t(aytemp)-2*(xG%*%t(ytemp))
      M2 <- a2%*%t(rep(1,nytemp))+rep(1,nx)%*%t(aytemp2)-2*(xG2%*%t(ytemp))
      eta <- eta + (2*pi)^(-d/2)*detG^(-1/2)*sum(exp(-M/2)*(M2-trGinv))
      if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
    }
  }
  else if (r==2)
  {
    xG <- x%*%Ginv
    xG2 <- x%*%G2inv
    xG3 <- x%*%G3inv
    a <- rowSums(xG*x)
    a2 <- rowSums(xG2*x)
    a3 <- rowSums(xG3*x)
    
    eta <- 0
    for (i in 1:(length(n.seq)-1))
    {
      nytemp <- n.seq[i+1] - n.seq[i]
      ytemp <- matrix(y[n.seq[i]:(n.seq[i+1]-1),], ncol=d)
      aytemp <- rowSums((ytemp %*% Ginv) *ytemp)
      aytemp2 <- rowSums((ytemp %*% G2inv) *ytemp)
      aytemp3 <- rowSums((ytemp %*% G3inv) *ytemp)
      M  <- a%*%t(rep(1,nytemp))+rep(1,nx)%*%t(aytemp)-2*(xG%*%t(ytemp))
      M2 <- a2%*%t(rep(1,nytemp))+rep(1,nx)%*%t(aytemp2)-2*(xG2%*%t(ytemp))
      M3 <- a3%*%t(rep(1,nytemp))+rep(1,nx)%*%t(aytemp3)-2*(xG3%*%t(ytemp))
      eta <- eta + (2*pi)^(-d/2)*detG^(-1/2)*sum(exp(-M/2)*(2*trG2inv-4*M3
             +(-trGinv+M2)^2))
      if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
    }
  }
  else if (r>2)
  {
    xG <- x%*%Ginv
    a <- rowSums(xG*x)
    eta <- 0
    for (i in 1:(length(n.seq)-1))
    {
      nytemp <- n.seq[i+1] - n.seq[i]
      ytemp <- matrix(y[n.seq[i]:(n.seq[i+1]-1),], ncol=d)
      aytemp <- rowSums((ytemp %*% Ginv) *ytemp)
      M <- a %*% t(rep(1,nytemp)) + rep(1,nx)%*%t(aytemp) - 2*(xG%*%t(ytemp))
      edv2 <- exp(-M/2)
      
      P0<-Ginv
      kappas <- matrix(nrow=as.numeric(nx*nytemp), ncol=r)
      for (j in 1:r)
      {
        Gi1inv <- P0%*%Ginv
        trGi0inv <- sum(diag(P0))    
        xGi1inv <- x%*%Gi1inv
        xGi1invx <- rowSums(xGi1inv*x)
        aytemp <- rowSums((ytemp %*% Gi1inv) *ytemp)
        dvi1 <- xGi1invx%*%t(rep(1,nytemp))+rep(1,nx)%*%t(aytemp)-2*(xGi1inv%*%t(ytemp))
        kappas[,j] <- (-2)^(j-1)*factorial(j-1)*(-trGi0inv+j*dvi1)
        P0 <- Gi1inv
      }
      
      nus <- matrix(nrow=as.numeric(nx*nytemp), ncol=r+1)        
      nus[,1] <- 1        
      for (j in 1:r)
      {
        js<-0:(j-1)
        if (j==1) nus[,2] <- kappas[,1]
        else nus[,j+1] <- rowSums(kappas[,j:1]*nus[,1:j]/matrix(rep(factorial(js)*
             factorial(rev(js)),nx*nytemp),nrow=nx*nytemp,byrow=TRUE))*factorial(j-1)
      }
      eta <- eta + (2*pi)^(-d/2)*detG^(-1/2)*sum(edv2*nus[,r+1])
      if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
    }
  }
  if (verbose) close(pb)
  if (inc==0) eta <- (eta - (-1)^r*nx*nur.cumulant(r=r, A=Ginv, mu=rep(0,d), Sigma=diag(d))*(2*pi)^(-d/2)*detG^(-1/2))/(nx*(ny-1))
  if (inc==1) eta <- eta/(nx*ny)  
  return(eta)
}


Qr.1d <- function(x, y, sigma, deriv.order=0, inc=1, verbose=FALSE)
{
  d <- 1
  r <- deriv.order/2
  if (missing(y)) y <- x
  nx <- length(x)
  ny <- length(y)
  g <- sigma
  
  n.seq <- block.indices(nx, ny, d=1, r=0, diff=FALSE)
  eta <- 0
  if (verbose) pb <- txtProgressBar() 
  
  if (r==0)
  {
    a <- x^2
    for (i in 1:(length(n.seq)-1))
    {
      if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
      nytemp <- n.seq[i+1] - n.seq[i]
      ytemp <- y[n.seq[i]:(n.seq[i+1]-1)]
      aytemp <- ytemp^2
      M <- a %*%t(rep(1,nytemp)) + rep(1, nx)%*%t(aytemp) - 2*(x %*% t(ytemp))
      em2 <- exp(-M/(2*g^2))
      eta <- eta + (2*pi)^(-d/2)*g^(-1)*sum(em2)
    }
  }
  else if (r>0)
  {
    a <- x^2
    for (i in 1:(length(n.seq)-1))
    {
      if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
      nytemp <- n.seq[i+1] - n.seq[i]
      ytemp <- y[n.seq[i]:(n.seq[i+1]-1)]
      aytemp <- ytemp^2 
      M <- a %*% t(rep(1,nytemp)) + rep(1,nx)%*%t(aytemp) - 2*(x %*%t(ytemp))
      edv2 <- exp(-M/(2*g^2))

      kappas <- matrix(nrow=as.numeric(nx*nytemp), ncol=r)
      for (i in 1:r)
      {
        aytemp <- ytemp^2
        dvi1 <- (a %*% t(rep(1,nytemp)) + rep(1,nx) %*% t(aytemp) - 2*(x%*%t(ytemp)))/g^(2*(i+1))
        kappas[,i] <- (-2)^(i-1)*factorial(i-1)*(-g^(-2*i)+i*dvi1)
      }
      
      nus <- matrix(nrow=as.numeric(nx*nytemp), ncol=r+1)        
      nus[,1] <- 1        
      for (j in 1:r)
      {
        js<-0:(j-1)
        if (j==1) nus[,2] <- kappas[,1]
        else nus[,j+1] <- rowSums(kappas[,j:1]*nus[,1:j]/matrix(rep(factorial(js)*factorial(rev(js)),nx*nytemp),nrow=nx*nytemp,byrow=TRUE))*factorial(j-1)
      }
      eta <- eta + (2*pi)^(-d/2)*g^(-1)*sum(edv2*nus[,r+1])
    }
  }
  if (verbose) close(pb)
  if (inc==0) eta <- (eta - nx*dnorm.deriv(x=0, mu=0, sigma=g, deriv.order=deriv.order))/(nx*(ny-1))
  if (inc==1) eta <- eta/(nx*ny) 
  
  return(eta)

}



###############################################################################
## Compute moments of multivariate normal mixture
###############################################################################

moments.mixt <- function (mus, Sigmas, props)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one")
  d <- ncol(Sigmas)
  k <- length(props)
  mn <- rep(0, d)
  va <- matrix(0, nrow=d, ncol=d)
  for (i in 1:k)
  {
    mn <- mn + props[i] * mus[i,]
    va <- va + props[i] * (Sigmas[((i-1)*d+1):(i*d),] + mus[i,] %*% t(mus[i,]))
  } 
  va <- va + mn %*% t(mn)
  return( list(mean=mn, var=va))
}


###############################################################################
## Creates plots of mixture density functions
#
## Parameters
## mus - means
## Sigmas - variances
## props - vector of proportions of each mixture component 
## dfs - degrees of freedom
## dist - "normal" - normal mixture
##      - "t" - t mixture
## ...
###############################################################################


plotmixt <- function(mus, sigmas, Sigmas, props, dfs, dist="normal", draw=TRUE, deriv.order=0, which.deriv.ind=1, binned=TRUE, ...)
{
  ## locally set random seed not to interfere with global random number generators
  if (!exists(".Random.seed")) rnorm(1)
  old.seed <- .Random.seed
  on.exit( { .Random.seed <<- old.seed } )
  set.seed(8192)
  
  if (!missing(sigmas)) plotmixt.1d(mus=mus, sigmas=sigmas, props=props, dfs=dfs, dist=dist, draw=draw, deriv.order=deriv.order, which.deriv.ind=which.deriv.ind, ...)
  else if (ncol(Sigmas)==2)
    plotmixt.2d(mus=mus, Sigmas=Sigmas, props=props, dfs=dfs, dist=dist, draw=draw, deriv.order=deriv.order, which.deriv.ind=which.deriv.ind,binned=binned, ...)
  else if (ncol(Sigmas)==3)
    plotmixt.3d(mus=mus, Sigmas=Sigmas, props=props, dfs=dfs, dist=dist, draw=draw, deriv.order=deriv.order, which.deriv.ind=which.deriv.ind, binned=binned, ...)
}

plotmixt.1d <- function(mus, sigmas, props, dfs, dist="normal", xlim, ylim, gridsize, draw=TRUE, deriv.order, which.deriv.ind, ...)
{
  dist1 <- match.arg(dist, c("normal", "t")) 
  maxsigmas <- 4*max(sigmas)

  if (missing(xlim)) xlim <- c(min(mus) - maxsigmas, max(mus) + maxsigmas)  
  if (missing(gridsize)) gridsize <- default.gridsize(1)

  x <- seq(xlim[1]-0.1*abs(diff(xlim)), xlim[2]+0.1*abs(diff(xlim)), length=gridsize)
  if (dist1=="normal")
  {
    if (deriv.order<=0) dens <- dnorm.mixt(x=x, mus=mus, sigmas=sigmas, props=props)
    else  dens <- dnorm.deriv.mixt(x=x, mus=mus, sigmas=sigmas, props=props, deriv.order=deriv.order)
  }
  else if (dist1=="t") stop("1-d t mixture not yet implemented")

  fhat <- list()
  fhat$x <- x
  fhat$eval.points <- x
  fhat$estimate <- dens
  fhat$H <- diag(1)
  fhat$h <- 1
  fhat$gridtype <- "linear"
  fhat$gridded <- TRUE
  fhat$binned <- FALSE
  fhat$names <- parse.name(x)
  fhat$w <- rep(1, length(x))
  if (deriv.order>0)
  {
    fhat$deriv.order <- deriv.order
    fhat$deriv.ind <- which.deriv.ind
  }
  class(fhat) <- "kdde"
  
  if (draw) plot(fhat, xlim=xlim, ...)
  invisible(fhat)
}



plotmixt.2d <- function(mus, Sigmas, props, dfs, dist="normal", xlim, ylim, gridsize, nrand=1e4, draw=TRUE, binned, deriv.order, which.deriv.ind, display="slice", ...)
{
  dist1 <- match.arg(dist, c("normal", "t"))
  disp1 <- match.arg(display, c("slice", "image", "persp", "filled.contour", "filled.contour2"))  
  maxSigmas <- 4*max(Sigmas)

  if (is.vector(mus)) mus <- as.matrix(t(mus))
  if (missing(xlim)) xlim <- c(min(mus[,1]) - maxSigmas, max(mus[,1]) + maxSigmas)
  if (missing(ylim)) ylim <- c(min(mus[,2]) - maxSigmas, max(mus[,2]) + maxSigmas)
  if (missing(gridsize)) gridsize <- default.gridsize(2)

  x <- seq(xlim[1], xlim[2], length=gridsize[1])
  y <- seq(ylim[1], ylim[2], length=gridsize[2])
  xy <- permute(list(x, y))
  d <- ncol(Sigmas)
 
  if (dist1=="normal")
  {
    if (deriv.order<=0) dens <- dmvnorm.mixt(xy, mus=mus, Sigmas=Sigmas, props=props)
    else  dens <- dmvnorm.deriv.mixt(xy, mus=mus, Sigmas=Sigmas, props=props, deriv.order=deriv.order)
  }
  else if (dist1=="t")
  {
    if (deriv.order>0) stop("deriv.order>0 for t mixture not yet implemented")
    dens <- dmvt.mixt(xy, mus=mus, Sigmas=Sigmas, props=props, dfs=dfs)
  }
  
  if (deriv.order<=0) dens.mat <- matrix(dens, ncol=length(x), byrow=FALSE)
  else
  {
    dens.mat <- list()
    for (i in 1:ncol(dens))
      dens.mat[[i]] <- matrix(dens[,i], ncol=length(x), byrow=FALSE)
  }
  
  if (dist1=="normal")
    x.rand <- rmvnorm.mixt(n=nrand, mus=mus, Sigmas=Sigmas, props=props)
  else if (dist1=="t")
    x.rand <- rmvt.mixt(n=nrand, mus=mus, Sigmas=Sigmas, props=props, dfs=dfs)
  
  H <- Hns(x=x.rand, deriv.order=deriv.order)
  if (binned) H <- diag(diag(H))
  fhat.rand <- kdde(x=x.rand, H=H, deriv.order=deriv.order, binned=binned)
  fhat <- fhat.rand 
  fhat$x <- x.rand
  fhat$eval.points <- list(x,y)
  fhat$estimate <- dens.mat
  fhat$names <- c("x", "y")

  if (deriv.order>0)
  {
    deriv.ind <- dmvnorm.deriv.mixt(xy, mus=mus, Sigmas=Sigmas, props=props, add.index=TRUE, only.index=TRUE, deriv.order=deriv.order, deriv.vec=TRUE)
    fhat$deriv.order <- deriv.order
    fhat$deriv.ind <- deriv.ind
    class(fhat) <- "kdde"
    if (draw)
    {
      plot(fhat, which.deriv.ind=which.deriv.ind, xlim=xlim, ylim=ylim, ...)
    }
  }
  else
  {
    if (draw)
    {
      if (disp1=="persp") plot(fhat, display=display, ...)
      else plot(fhat, xlim=xlim, ylim=ylim, display=display, ...)
    }
  }
  invisible(fhat)
}


plotmixt.3d <- function(mus, Sigmas, props, dfs, dist="normal", xlim, ylim, zlim, gridsize, nrand=1e4, draw=TRUE, binned, deriv.order, which.deriv.ind, ...)
{
  d <- 3
  dist1 <- match.arg(dist, c("normal", "t")) 
  maxsd <- sqrt(apply(Sigmas, 2, max))

  if (is.vector(mus)) mus <- as.matrix(t(mus))
  if (missing(xlim)) xlim <- c(min(mus[,1]) - 4*maxsd[1], max(mus[,1]) + 4*maxsd[1])
  if (missing(ylim)) ylim <- c(min(mus[,2]) - 4*maxsd[2], max(mus[,2]) + 4*maxsd[2])
  if (missing(zlim)) zlim <- c(min(mus[,3]) - 4*maxsd[3], max(mus[,3]) + 4*maxsd[3])
  if (missing(gridsize)) gridsize <- default.gridsize(3)
  
  x <- seq(xlim[1], xlim[2], length=gridsize[1])
  y <- seq(ylim[1], ylim[2], length=gridsize[2])
  z <- seq(zlim[1], zlim[2], length=gridsize[3])
  xy <- permute(list(x,y))

  if (deriv.order>0)
  {
    if (dist1=="t")
      stop("deriv.order>0 for t mixture not yet implemented")
    else if (dist1=="normal")
      deriv.ind <- dmvnorm.deriv.mixt(cbind(xy,z[1]), mus=mus, Sigmas=Sigmas, props=props, deriv.order=deriv.order, add.index=TRUE, only.index=TRUE)
  }
     
  if (deriv.order<=0) dens.array <- array(0, dim=gridsize)
  else { dens.array <- list(); for (i in 1:nrow(deriv.ind)) dens.array <- c(dens.array, list(array(0, dim=gridsize))) }

  
  for (i in 1:length(z))
  {
    if (dist1=="normal")
    {
      if (deriv.order<=0)
        dens <- dmvnorm.mixt(cbind(xy, z[i]), mus=mus, Sigmas=Sigmas, props=props)
      else
        dens <- dmvnorm.deriv.mixt(cbind(xy, z[i]), mus=mus, Sigmas=Sigmas, props=props, deriv.order=deriv.order)
    }
    else if (dist1=="t")
      dens <- dmvt.mixt(cbind(xy, z[i]), mus=mus, Sigmas=Sigmas, dfs=dfs, props=props)

    if (deriv.order<=0)
    {
      dens.mat <- matrix(dens, ncol=length(x), byrow=FALSE)
      dens.array[,,i] <- dens.mat
    }
    else
    {
      for (j in 1:ncol(dens))
      {
        dens.mat <- matrix(dens[,j], ncol=length(x), byrow=FALSE)
        dens.array[[j]][,,i] <- dens.mat
      }
    }
  
  }

  if (dist1=="normal")
    x.rand <- rmvnorm.mixt(n=nrand, mus=mus, Sigmas=Sigmas, props=props)
  else if (dist1=="t")
    x.rand <- rmvt.mixt(n=nrand, mus=mus, Sigmas=Sigmas, props=props, dfs=dfs)

  H <- Hns(x=x.rand, deriv.order=deriv.order)
  if (binned) H <- diag(diag(H))
  fhat.rand <- kdde(x=x.rand, H=H, deriv.order=deriv.order, binned=binned)
    
  fhat <- fhat.rand
  fhat$x <- head(x.rand, n=100)
  fhat$eval.points <- list(x,y,z)
  fhat$estimate <- dens.array
  fhat$names <- c("x", "y", "z")
  fhat$H <- H
  fhat$w <- rep(1,nrow(fhat$x))
  
  if (deriv.order>0)
  {
    deriv.ind <- dmvnorm.deriv.mixt(xy, mus=mus, Sigmas=Sigmas, props=props, add.index=TRUE, only.index=TRUE, deriv.order=deriv.order, deriv.vec=TRUE)
    fhat$deriv.order <- deriv.order
    fhat$deriv.ind <- deriv.ind
    class(fhat) <- "kdde"
    if (draw) plot(fhat, which.deriv.ind=which.deriv.ind, xlim=xlim, ylim=ylim, zlim=zlim, ...)
  }
  else
  {
    if (draw) plot(fhat, xlim=xlim, ylim=ylim, zlim=zlim, ...)
  }
 
  invisible(fhat)
}







###############################################################################
## Multivariate t mixture - density values
##
## Parameters
## x - points to compute density at    
## mus - vector of means 
## Sigmas - dispersion matrices
## dfs - degrees of freedom
## props - vector of mixing proportions
##
## Returns
## Value of multivariate t mixture density at x
###############################################################################


dmvt.mixt <- function(x, mus, Sigmas, dfs, props)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one")
  else if (length(dfs) != length(props))
    stop("Length of df and mixing proportions vectors not equal")
  
  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
    dens <- dmvt(x, delta=mus, sigma=Sigmas, df=dfs, log=FALSE)
  
  ## multiple component mixture
  else   
  {   
    if (is.vector(mus)) d <- length(mus)
    else d <- ncol(mus)
    k <- length(props)
    dens <- 0      
    for (i in 1:k)
      dens <- dens+props[i]*dmvt(x,delta=mus[i,],sigma=Sigmas[((i-1)*d+1):(i*d),], df=dfs[i], log=FALSE)
  }
  
  return(dens)
}


###############################################################################
## Multivariate t mixture - random sample
## 
## Parameters
## n - number of samples
## mus - means 
## Sigmas - matrix of dispersion matrices
## dfs - vector of degrees of freedom
## props - vector of mixing proportions 
## 
## Returns
## Vector of n observations from the t mixture
###############################################################################

rmvt.mixt <- function(n=100, mus=c(0,0), Sigmas=diag(2), dfs=7, props=1)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))  
    stop("Proportions don't sum to one")
  else if (length(dfs) != length(props))
    stop("Length of df and mixing proportions vectors not equal")  

  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
  {
    rand <- rmvt(n=n, sigma=Sigmas, df=dfs)
    for (i in 1:length(mus))
      rand[,i] <- rand[,i] + mus[i]
  }
  
  ## multiple component mixture
  else
  {
    k <- length(props)
    d <- ncol(Sigmas)
    n.samp <- sample(1:k, n, replace=TRUE, prob=props) 
    n.prop <- numeric(0)

    ## compute number to be drawn from each component 
    for (i in 1:k)
      n.prop <- c(n.prop, sum(n.samp == i))

    ## generate random samples from each component
    rand <- numeric(0)  
    for (i in 1:k)
    {
      if (n.prop[i] > 0)
      {  
        rand.temp<-rmvt(n=n.prop[i],sigma=Sigmas[((i-1)*d+1):(i*d),],df=dfs[i])
        for (j in 1:length(mus[k,]))
          rand.temp[,j] <- rand.temp[,j] + mus[i,j]
       
        rand <- rbind(rand, rand.temp)
      }
    }
  }
  
  return(rand[sample(n),])
}


