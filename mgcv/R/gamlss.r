## (c) Simon N. Wood (2013,2014) distributed under GPL2
## Code for the gamlss families.
## idea is that there are standard functions converting
## derivatives w.r.t. mu to derivatives w.r.t. eta, given 
## given the links and derivatives. 
## Then there are standard routines to take the family 
## specific derivatives and the model matrices, and convert 
## these to the required gradient, hessian, etc...

## Storage convections:
## 1. A single model matrix is stored, along with a single param vector.
##    an index object associates columns with particular gamlss predictors.
## 2. Distribution specific derivatives are stored in d1l-d4l.  

## Need to somehow record block starts... 
## idea is that if n blocks are stored using loops with the 
## given l >= k >= j >= i structure then the block for index
## i,j,k,l starts at i4[i,j,k,l]*n+1, given symmetry over the indices. 

trind.generator <- function(K=2) {
## Generates index arrays for 'upper triangular' storage up to order 4
## Suppose you fill an array using code like...
## m = 1
## for (i in 1:K) for (j in i:K) for (k in j:K) for (l in k:K) {
##   a[,m] <- something; m <- m+1 }
## ... and do this because actually the same 'something' would 
## be stored for any permutation of the indices i,j,k,l.
## Clearly in storage we have the restriction l>=k>=j>=i,
## but for access we want no restriction on the indices.
## i4[i,j,k,l] produces the appropriate m for unrestricted 
## indices. i3 and i2 do the same for 3d and 2d arrays.
  i4 <- array(0,dim=c(K,K,K,K))
  m.start <- 1
  m <- m.start
  for (i in 1:K) for (j in i:K) for (k in j:K) for (l in k:K) {
    i4[i,j,k,l] <- i4[i,j,l,k] <- i4[i,k,l,j] <- i4[i,k,j,l] <- i4[i,l,j,k] <- 
    i4[i,l,k,j] <- 
    i4[j,i,k,l] <- i4[j,i,l,k] <- i4[j,k,l,i] <- i4[j,k,i,l] <- i4[j,l,i,k] <- 
    i4[j,l,k,i] <- 
    i4[k,j,i,l] <- i4[k,j,l,i] <- i4[k,i,l,j] <- i4[k,i,j,l] <- i4[k,l,j,i] <- 
    i4[k,l,i,j] <- 
    i4[l,j,k,i] <- i4[l,j,i,k] <- i4[l,k,i,j] <- i4[l,k,j,i] <- i4[l,i,j,k] <- 
    i4[l,i,k,j] <- m
    m <- m + 1
  }

  i3 <- array(0,dim=c(K,K,K))
  m <- m.start
  for (j in 1:K) for (k in j:K) for (l in k:K) {
    i3[j,k,l] <- i3[j,l,k] <- i3[k,l,j] <- i3[k,j,l] <- i3[l,j,k] <- 
    i3[l,k,j] <- m
    m <- m + 1
  }

  i2 <- array(0,dim=c(K,K))
  m <- m.start
  for (k in 1:K) for (l in k:K) {
    i2[k,l] <- i2[l,k] <- m
    m <- m + 1
  }
  
  list(i2=i2,i3=i3,i4=i4)
} ## trind.generator

gamlss.etamu <- function(l1,l2,l3=NULL,l4=NULL,ig1,g2,g3=NULL,g4=NULL,i2,i3=NULL,i4=NULL,deriv=0) {
## lj is the array of jth order derivatives of l
## gj[,k] contains the jth derivatives for the link of the kth lp
## ig1 is one over first deriv of link
## kth parameter. This routine transforms derivatives 
## w.r.t. the parameters (mu_1..mu_K) to derivatives
## w.r.t. the linear predictors (eta_1.. eta_K)
## i2, i3 and i4 are the upper triangular indexing arrays
## e.g. l4[,i4[i,j,l,m]] contains the partial w.r.t.  
## params indexed by i,j,l,m with no restriction on
## the index values except that they are in 1..K 
  K <- ncol(l1) ## number of parameters of distribution
  d1 <- l1
  for (i in 1:K) { ## first derivative loop
    d1[,i] <- l1[,i]*ig1[,i]
  }

  ##n <- length(ig1[,1])

  k <- 0
  d2 <- l2
  for (i in 1:K) for (j in i:K) {
    ## obtain the order of differentiation associated 
    ## with the i,j derivatives...
    ord <- rep(1,2);k <- k+1
    if (i==j) {ord[1] <- ord[1] + 1; ord[2] <- 0 }
    ## l2[,k] is derivative to transform
    mo <- max(ord)
    if (mo==2) { ## pure 2nd derivative transform
      d2[,k] <- (l2[,k] - l1[,i]*g2[,i]*ig1[,i])*ig1[,i]^2   
    } else { ## all first derivative
      d2[,k] <- l2[,k]*ig1[,i]*ig1[,j]
    }
  } ## 2nd order transform done


  k <- 0
  d3 <- l3
  if (deriv>0) for (i in 1:K) for (j in i:K) for (l in j:K) {
    ## obtain the order of differentiation associated 
    ## with the i,j,l derivatives...
    ord <- rep(1,3);k <- k+1
    if (i==j) {ord[1] <- ord[1] + 1; ord[2] <- 0 }
    if (i==l) {ord[1] <- ord[1] + 1; ord[3] <- 0 }
    if (ord[2]) {
      if (j==l) {ord[2] <- ord[2] + 1; ord[3] <- 0 }
    }
    ii <- c(i,j,l)
    ## l3[,k] is derivative to transform
    mo <- max(ord)
    if (mo==3) { ## pure 3rd derivative transform
      d3[,k] <- (l3[,k] - 3*l2[,i2[i,i]]*g2[,i]*ig1[,i] +
                l1[,i]*(3*g2[,i]^2*ig1[,i]^2 - g3[,i]*ig1[,i]))*ig1[,i]^3          
    } else if (mo==1) { ## all first derivative
      d3[,k] <- l3[,k]*ig1[,i]*ig1[,j]*ig1[,l]
    } else { ## 2,1 deriv
      k1 <- ii[ord==1] ## index of order 1 deriv
      k2 <- ii[ord==2] ## index of order 2 part
      d3[,k] <- (l3[,k] - l2[,i2[k2,k1]]*g2[,k2]*ig1[,k2])*
                ig1[,k1]*ig1[,k2]^2
    } 
  } ## 3rd order transform done
  
  k <- 0
  d4 <- l4
  if (deriv>2) for (i in 1:K) for (j in i:K) for (l in j:K) for (m in l:K) {
    ## obtain the order of differentiation associated 
    ## with the i,j,l & m derivatives...
    ord <- rep(1,4);k <- k+1
    if (i==j) {ord[1] <- ord[1] + 1; ord[2] <- 0 }
    if (i==l) {ord[1] <- ord[1] + 1; ord[3] <- 0 }
    if (i==m) {ord[1] <- ord[1] + 1; ord[4] <- 0 }
    if (ord[2]) {
      if (j==l) {ord[2] <- ord[2] + 1; ord[3] <- 0 }
      if (j==m) {ord[2] <- ord[2] + 1; ord[4] <- 0 }
    }
    if (ord[3]&&l==m) { ord[3] <- ord[3] + 1; ord[4] <- 0 }
    ii <- c(i,j,l,m)
    ## l4[,k] is derivative to transform
    mo <- max(ord)
    if (mo==4) { ## pure 4th derivative transform
    d4[,k] <-  (l4[,k] - 6*l3[,i3[i,i,i]]*g2[,i]*ig1[,i] + 
        l2[,i2[i,i]]*(15*g2[,i]^2*ig1[,i]^2 - 4*g3[,i]*ig1[,i]) - 
        l1[,i]*(15*g2[,i]^3*ig1[,i]^3 - 10*g2[,i]*g3[,i]*ig1[,i]^2 
         + g4[,i]*ig1[,i]))*ig1[,i]^4    
    } else if (mo==1) { ## all first derivative
      d4[,k] <- l4[,k]*ig1[,i]*ig1[,j]*ig1[,l]*ig1[,m]
    } else if (mo==3) { ## 3,1 deriv
      k1 <- ii[ord==1] ## index of order 1 deriv
      k3 <- ii[ord==3] ## index of order 3 part
      d4[,k] <- (l4[,k] - 3*l3[,i3[k3,k3,k1]]*g2[,k3]*ig1[,k3] +
        l2[,i2[k3,k1]]*(3*g2[,k3]^2*ig1[,k3]^2 - g3[,k3]*ig1[,k3])         
      )*ig1[,k1]*ig1[,k3]^3
    } else { 
      if (sum(ord==2)==2) { ## 2,2
        k2a <- (ii[ord==2])[1];k2b <- (ii[ord==2])[2]
        d4[,k] <- (l4[,k] - l3[,i3[k2a,k2b,k2b]]*g2[,k2a]*ig1[,k2a]
          -l3[,i3[k2a,k2a,k2b]]*g2[,k2b]*ig1[,k2b] + 
           l2[,i2[k2a,k2b]]*g2[,k2a]*g2[,k2b]*ig1[,k2a]*ig1[,k2b]
        )*ig1[,k2a]^2*ig1[,k2b]^2
      } else { ## 2,1,1
        k2 <- ii[ord==2] ## index of order 2 derivative
        k1a <- (ii[ord==1])[1];k1b <- (ii[ord==1])[2]
        d4[,k] <- (l4[,k] - l3[,i3[k2,k1a,k1b]]*g2[,k2]*ig1[,k2] 
                   )*ig1[,k1a]*ig1[,k1b]*ig1[,k2]^2
      }
    }
  } ## 4th order transform done

  list(l1=d1,l2=d2,l3=d3,l4=d4)
} # gamlss.etamu


gamlss.gH0 <- function(X,jj,l1,l2,i2,l3=0,i3=0,l4=0,i4=0,d1b=0,d2b=0,deriv=0,fh=NULL,D=NULL) {
## X[,jj[[i]]] is the ith model matrix.
## lj contains jth derivatives of the likelihood for each datum,
## columns are w.r.t. different combinations of parameters.
## ij is the symmetric array indexer for the jth order derivs...
## e.g. l4[,i4[i,j,l,m]] contains derivatives with 
## respect to parameters indexed by i,j,l,m
## d1b and d2b are first and second derivatives of beta w.r.t. sps.
## fh is a factorization of the penalized hessian, while D contains the corresponding
##    Diagonal pre-conditioning weights.
## deriv: 0 - just grad and Hess
##        1 - diagonal of first deriv of Hess
##        2 - first deriv of Hess
##        3 - everything.
  K <- length(jj)
  p <- ncol(X);n <- nrow(X)
  trHid2H <- d1H <- d2H <- NULL ## defaults

  ## the gradient...
  lb <- rep(0,p)
  for (i in 1:K) { ## first derivative loop
    lb[jj[[i]]] <- colSums(l1[,i]*X[,jj[[i]],drop=FALSE])
  }
  
  ## the Hessian...
  lbb <- matrix(0,p,p)
  for (i in 1:K) for (j in i:K) {
    lbb[jj[[i]],jj[[j]]] <- t(X[,jj[[i]],drop=FALSE])%*%(l2[,i2[i,j]]*X[,jj[[j]],drop=FALSE])
    lbb[jj[[j]],jj[[i]]] <- t(lbb[jj[[i]],jj[[j]]])
  } 

  if (deriv>0) {
    ## the first derivative of the Hessian, using d1b
    ## the first derivates of the coefficients wrt the sps
    m <- ncol(d1b) ## number of smoothing parameters
    ## stack the derivatives of the various linear predictors on top
    ## of each other...
    d1eta <- matrix(0,n*K,m)
    ind <- 1:n
    for (i in 1:K) { 
      d1eta[ind,] <- X[,jj[[i]],drop=FALSE]%*%d1b[jj[[i]],]
      ind <- ind + n
    }
  }

  if (deriv==1) { 
    d1H <- matrix(0,p,m) ## only store diagonals of d1H
    for (l in 1:m) {
      for (i in 1:K) {
        v <- rep(0,n);ind <- 1:n
        for (q in 1:K) { 
          v <- v + l3[,i3[i,i,q]] * d1eta[ind,l]
          ind <- ind + n
        }
        d1H[jj[[i]],l] <- colSums(X[,jj[[i]],drop=FALSE]*(v*X[,jj[[i]],drop=FALSE]))
      } 
    }
  } ## if deriv==1

  if (deriv>1) {
    d1H <- list()
    for (l in 1:m) {
      d1H[[l]] <- matrix(0,p,p)
      for (i in 1:K) for (j in i:K) {
        v <- rep(0,n);ind <- 1:n
        for (q in 1:K) { 
          v <- v + l3[,i3[i,j,q]] * d1eta[ind,l]
          ind <- ind + n
        }
        ## d1H[[l]][jj[[j]],jj[[i]]] <- 
        d1H[[l]][jj[[i]],jj[[j]]] <- t(X[,jj[[i]],drop=FALSE])%*%(v*X[,jj[[j]],drop=FALSE])
        d1H[[l]][jj[[j]],jj[[i]]] <- t(d1H[[l]][jj[[i]],jj[[j]]])
      } 
    }
  } ## if deriv>1

  if (deriv>2) {
    ## need tr(Hp^{-1} d^2H/drho_k drho_j)
    ## First form the expanded model matrix...
    VX <- Xe <- matrix(0,K*n,ncol(X))
    ind <- 1:n 
    for (i in 1:K) { 
      Xe[ind,jj[[i]]] <- X[,jj[[i]]]
      ind <- ind + n
    }
    ## Now form Hp^{-1} Xe'...
    if (is.list(fh)) { ## then the supplied factor is an eigen-decomposition
       d <- fh$values;d[d>0] <- 1/d[d>0];d[d<=0] <- 0
       Xe <- t(D*((fh$vectors%*%(d*t(fh$vectors)))%*%(D*t(Xe))))
    } else { ## the supplied factor is a choleski factor
       ipiv <- piv <- attr(fh,"pivot");ipiv[piv] <- 1:p
       Xe <- t(D*(backsolve(fh,forwardsolve(t(fh),(D*t(Xe))[piv,]))[ipiv,]))
    }
    ## now compute the required trace terms
    d2eta <- matrix(0,n*K,ncol(d2b))
    ind <- 1:n
    for (i in 1:K) { 
      d2eta[ind,] <- X[,jj[[i]],drop=FALSE]%*%d2b[jj[[i]],]
      ind <- ind + n
    }
    trHid2H <- rep(0,ncol(d2b))
    kk <- 0 ## counter for second derivatives
    for (k in 1:m) for (l in k:m) { ## looping over smoothing parameters...
      kk <- kk + 1
      for (i in 1:K) for (j in 1:K) {
        v <- rep(0,n);ind <- 1:n
        for (q in 1:K) { ## accumulate the diagonal matrix for X_i'diag(v)X_j
          v <- v + d2eta[ind,kk]*l3[,i3[i,j,q]]
          ins <- 1:n
          for (s in 1:K) { 
            v <- v + d1eta[ind,k]*d1eta[ins,l]*l4[,i4[i,j,q,s]]
            ins <- ins + n
          }
          ind <- ind + n
        }
        if (i==j) {
          rind <- 1:n + (i-1)*n
          VX[rind,jj[[i]]] <- v * X[,jj[[i]]]
        } else {
          rind1 <- 1:n + (i-1)*n
          rind2 <- 1:n + (j-1)*n
          VX[rind2,jj[[i]]] <- v * X[,jj[[i]]]
          VX[rind1,jj[[j]]] <- v * X[,jj[[j]]]
        }
      }
      trHid2H[kk] <- sum(Xe*VX)
    }
  } ## if deriv>2

  list(lb=lb,lbb=lbb,d1H=d1H,d2H=d2H,trHid2H=trHid2H)
} ## end of gamlss.gH0



gamlss.gH <- function(X,jj,l1,l2,i2,l3=0,i3=0,l4=0,i4=0,d1b=0,d2b=0,deriv=0,fh=NULL,D=NULL) {
## X[,jj[[i]]] is the ith model matrix.
## lj contains jth derivatives of the likelihood for each datum,
## columns are w.r.t. different combinations of parameters.
## ij is the symmetric array indexer for the jth order derivs...
## e.g. l4[,i4[i,j,l,m]] contains derivatives with 
## respect to parameters indexed by i,j,l,m
## d1b and d2b are first and second derivatives of beta w.r.t. sps.
## fh is a factorization of the penalized hessian, while D contains the corresponding
##    Diagonal pre-conditioning weights.
## deriv: 0 - just grad and Hess
##        1 - diagonal of first deriv of Hess
##        2 - first deriv of Hess
##        3 - everything.
  K <- length(jj)
  p <- ncol(X);n <- nrow(X)
  trHid2H <- d1H <- d2H <- NULL ## defaults

  ## the gradient...
  lb <- rep(0,p)
  for (i in 1:K) { ## first derivative loop
    lb[jj[[i]]] <- lb[jj[[i]]] + colSums(l1[,i]*X[,jj[[i]],drop=FALSE]) ## !
  }
  
  ## the Hessian...
  lbb <- matrix(0,p,p)
  for (i in 1:K) for (j in i:K) {
    A <- t(X[,jj[[i]],drop=FALSE])%*%(l2[,i2[i,j]]*X[,jj[[j]],drop=FALSE])
    lbb[jj[[i]],jj[[j]]] <- lbb[jj[[i]],jj[[j]]] + A 
    if (j>i) lbb[jj[[j]],jj[[i]]] <- lbb[jj[[j]],jj[[i]]] + t(A) 
  } 

  if (deriv>0) {
    ## the first derivative of the Hessian, using d1b
    ## the first derivates of the coefficients wrt the sps
    m <- ncol(d1b) ## number of smoothing parameters
    ## stack the derivatives of the various linear predictors on top
    ## of each other...
    d1eta <- matrix(0,n*K,m)
    ind <- 1:n
    for (i in 1:K) { 
      d1eta[ind,] <- X[,jj[[i]],drop=FALSE]%*%d1b[jj[[i]],]
      ind <- ind + n
    }
  }

  if (deriv==1) { 
    d1H <- matrix(0,p,m) ## only store diagonals of d1H
    for (l in 1:m) {
      for (i in 1:K) {
        v <- rep(0,n);ind <- 1:n
        for (q in 1:K) { 
          v <- v + l3[,i3[i,i,q]] * d1eta[ind,l]
          ind <- ind + n
        }
        d1H[jj[[i]],l] <-  d1H[jj[[i]],l] + colSums(X[,jj[[i]],drop=FALSE]*(v*X[,jj[[i]],drop=FALSE])) 
      } 
    }
  } ## if deriv==1

  if (deriv>1) {
    d1H <- list()
    for (l in 1:m) {
      d1H[[l]] <- matrix(0,p,p)
      for (i in 1:K) for (j in i:K) {
        v <- rep(0,n);ind <- 1:n
        for (q in 1:K) { 
          v <- v + l3[,i3[i,j,q]] * d1eta[ind,l]
          ind <- ind + n
        }
        ## d1H[[l]][jj[[j]],jj[[i]]] <- 
        A <- t(X[,jj[[i]],drop=FALSE])%*%(v*X[,jj[[j]],drop=FALSE])
        d1H[[l]][jj[[i]],jj[[j]]] <- d1H[[l]][jj[[i]],jj[[j]]] + A
        if (j>i) d1H[[l]][jj[[j]],jj[[i]]] <- d1H[[l]][jj[[j]],jj[[i]]] + t(A)
      } 
    }
  } ## if deriv>1

  if (deriv>2) {
    ## need tr(Hp^{-1} d^2H/drho_k drho_j)
    ## First form the expanded model matrix...
    VX <- Xe <- matrix(0,K*n,ncol(X))
    ind <- 1:n 
    for (i in 1:K) { 
      Xe[ind,jj[[i]]] <- X[,jj[[i]]]
      ind <- ind + n
    }
    ## Now form Hp^{-1} Xe'...
    if (is.list(fh)) { ## then the supplied factor is an eigen-decomposition
       d <- fh$values;d[d>0] <- 1/d[d>0];d[d<=0] <- 0
       Xe <- t(D*((fh$vectors%*%(d*t(fh$vectors)))%*%(D*t(Xe))))
    } else { ## the supplied factor is a choleski factor
       ipiv <- piv <- attr(fh,"pivot");ipiv[piv] <- 1:p
       Xe <- t(D*(backsolve(fh,forwardsolve(t(fh),(D*t(Xe))[piv,]))[ipiv,]))
    }
    ## now compute the required trace terms
    d2eta <- matrix(0,n*K,ncol(d2b))
    ind <- 1:n
    for (i in 1:K) { 
      d2eta[ind,] <- X[,jj[[i]],drop=FALSE]%*%d2b[jj[[i]],]
      ind <- ind + n
    }
    trHid2H <- rep(0,ncol(d2b))
    kk <- 0 ## counter for second derivatives
    for (k in 1:m) for (l in k:m) { ## looping over smoothing parameters...
      kk <- kk + 1
      for (i in 1:K) for (j in 1:K) {
        v <- rep(0,n);ind <- 1:n
        for (q in 1:K) { ## accumulate the diagonal matrix for X_i'diag(v)X_j
          v <- v + d2eta[ind,kk]*l3[,i3[i,j,q]]
          ins <- 1:n
          for (s in 1:K) { 
            v <- v + d1eta[ind,k]*d1eta[ins,l]*l4[,i4[i,j,q,s]]
            ins <- ins + n
          }
          ind <- ind + n
        }
        if (i==j) {
          rind <- 1:n + (i-1)*n
          VX[rind,jj[[i]]] <- v * X[,jj[[i]]]
        } else {
          rind1 <- 1:n + (i-1)*n
          rind2 <- 1:n + (j-1)*n
          VX[rind2,jj[[i]]] <- v * X[,jj[[i]]]
          VX[rind1,jj[[j]]] <- v * X[,jj[[j]]]
        }
      }
      trHid2H[kk] <- sum(Xe*VX)
    }
  } ## if deriv>2

  list(lb=lb,lbb=lbb,d1H=d1H,d2H=d2H,trHid2H=trHid2H)
} ## end of gamlss.gH


gaulss <- function(link=list("identity","logb"),b=0.01) {
## Extended family for Gaussian location scale model...
## so mu is mu1 and tau=1/sig is mu2
## tau = 1/(b + exp(eta)) eta = log(1/tau - b)
## 1. get derivatives wrt mu, tau
## 2. get required link derivatives and tri indices.
## 3. transform derivs to derivs wrt eta (gamlss.etamu).
## 4. get the grad and Hessian etc for the model
##    via a call to gamlss.gH  
## the first derivatives of the log likelihood w.r.t
## the first and second parameters...
  
  ## first deal with links and their derivatives...
  if (length(link)!=2) stop("gaulss requires 2 links specified as character strings")
  okLinks <- list(c("inverse", "log", "identity","sqrt"),"logb")
  stats <- list()
  if (link[[1]] %in% okLinks[[1]]) stats[[1]] <- make.link(link[[1]]) else 
  stop(link[[1]]," link not available for mu parameter of gaulss")
  fam <- structure(list(link=link[[1]],canonical="none",linkfun=stats[[1]]$linkfun,
           mu.eta=stats[[1]]$mu.eta),
           class="family")
  fam <- fix.family.link(fam)
  stats[[1]]$d2link <- fam$d2link
  stats[[1]]$d3link <- fam$d3link
  stats[[1]]$d4link <- fam$d4link
  if (link[[2]] %in% okLinks[[2]]) { ## creating the logb link
    stats[[2]] <- list()
    stats[[2]]$valideta <- function(eta) TRUE 
    stats[[2]]$link = link[[2]]
    stats[[2]]$linkfun <- eval(parse(text=paste("function(mu) log(1/mu -",b,")")))
    stats[[2]]$linkinv <- eval(parse(text=paste("function(eta) 1/(exp(eta) +",b,")")))
    stats[[2]]$mu.eta <- eval(parse(text=
                         paste("function(eta) { ee <- exp(eta); -ee/(ee +",b,")^2 }")))
    stats[[2]]$d2link <-  eval(parse(text=
    paste("function(mu) { mub <- 1 - mu *",b,";(2*mub-1)/(mub*mu)^2}" )))
    stats[[2]]$d3link <-  eval(parse(text=
    paste("function(mu) { mub <- 1 - mu *",b,";((1-mub)*mub*6-2)/(mub*mu)^3}" )))
    stats[[2]]$d4link <-  eval(parse(text=
    paste("function(mu) { mub <- 1 - mu *",b,";(((24*mub-36)*mub+24)*mub-6)/(mub*mu)^4}")))
  } else stop(link[[2]]," link not available for precision parameter of gaulss")
  
  residuals <- function(object,type=c("deviance","pearson","response")) {
      type <- match.arg(type)
      rsd <- object$y-object$fitted[,1]
      if (type=="response") return(rsd) else
      return((rsd*object$fitted[,2])) ## (y-mu)/sigma 
    }
  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    object$null.deviance <- sum(((object$y-mean(object$y))*object$fitted[,2])^2)
  })

  ll <- function(y,X,coef,wt,family,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
  ## function defining the gamlss Gaussian model log lik. 
  ## N(mu,sigma^2) parameterized in terms of mu and log(sigma)
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.
    jj <- attr(X,"lpi") ## extract linear predictor index
    eta <- X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]]
    mu <- family$linfo[[1]]$linkinv(eta)
    eta1 <- X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]] 
    tau <-  family$linfo[[2]]$linkinv(eta1) ## tau = 1/sig here
    
    n <- length(y)
    l1 <- matrix(0,n,2)
    ymu <- y-mu;ymu2 <- ymu^2;tau2 <- tau^2
 
    l <- sum(-.5 * ymu2 * tau2 - .5 * log(2*pi) + log(tau))

    if (deriv>0) {

      l1[,1] <- tau2*ymu
      l1[,2] <- 1/tau - tau*ymu2  

      ## the second derivatives
    
      l2 <- matrix(0,n,3)
      ## order mm,ms,ss
      l2[,1] <- -tau2
      l2[,2] <- 2*l1[,1]/tau
      l2[,3] <- -ymu2 - 1/tau2

      ## need some link derivatives for derivative transform
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),family$linfo[[2]]$d2link(tau))
    }
 
    l3 <- l4 <- g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
      ## the third derivatives
      ## order mmm,mms,mss,sss
      l3 <- matrix(0,n,4) 
      ## l3[,1] <- 0
      l3[,2] <- -2*tau
      l3[,3] <- 2*ymu
      l3[,4] <- 2/tau^3 
      g3 <- cbind(family$linfo[[1]]$d3link(mu),family$linfo[[2]]$d3link(tau))
    }

    if (deriv>3) {
      ## the fourth derivatives
      ## order mmmm,mmms,mmss,msss,ssss
      l4 <- matrix(0,n,5) 
      ## l4[,1] <- 0
      ## l4[,2] <- 0
      l4[,3] <- -2
      #l4[,4] <- 0
      l4[,5] <- -6/tau2^2 
      g4 <- cbind(family$linfo[[1]]$d4link(mu),family$linfo[[2]]$d4link(tau))
    }
    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3
      i4 <- family$tri$i4
   
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)

      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                      d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
    } else ret <- list()
    ret$l <- l; ret
  } ## end ll gaulss

  initialize <- expression({
  ## idea is to regress g(y) on model matrix for mean, and then 
  ## to regress the corresponding log absolute residuals on 
  ## the model matrix for log(sigma) - may be called in both
  ## gam.fit5 and initial.spg... note that appropriate E scaling
  ## for full calculation may be inappropriate for initialization 
  ## which is basically penalizing something different here.
  ## best we can do here is to use E only as a regularizer.
      n <- rep(1, nobs)
      ## should E be used unscaled or not?..
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      if (is.null(start)) {
        jj <- attr(x,"lpi")
        start <- rep(0,ncol(x))
        yt1 <- if (family$link[[1]]=="identity") y else 
               family$linfo[[1]]$linkfun(abs(y)+max(y)*1e-7)
        x1 <- x[,jj[[1]],drop=FALSE]
        e1 <- E[,jj[[1]],drop=FALSE] ## square root of total penalty
        #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
        if (use.unscaled) {
          qrx <- qr(rbind(x1,e1))
          x1 <- rbind(x1,e1)
          startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
          startji[!is.finite(startji)] <- 0       
        } else startji <- pen.reg(x1,e1,yt1)
        start[jj[[1]]] <- startji
        lres1 <- log(abs(y-family$linfo[[1]]$linkinv(x[,jj[[1]],drop=FALSE]%*%start[jj[[1]]])))
        x1 <-  x[,jj[[2]],drop=FALSE];e1 <- E[,jj[[2]],drop=FALSE]
        #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
        if (use.unscaled) {
          x1 <- rbind(x1,e1)
          startji <- qr.coef(qr(x1),c(lres1,rep(0,nrow(E))))   
          startji[!is.finite(startji)] <- 0
        } else startji <- pen.reg(x1,e1,lres1)
        start[jj[[2]]] <- startji
      }
  }) ## initialize gaulss

  structure(list(family="gaulss",ll=ll,link=paste(link),nlp=2,
    tri = trind.generator(2), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,
    linfo = stats, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    available.derivs = 2 ## can use full Newton here
    ),class = c("general.family","extended.family","family"))
} ## end gaulss


multinom <- function(K=1) {
## general family for multinomial logistic regression model...
## accepts no links as parameterization directly in terms of 
## linear predictor. 
## 1. get derivatives wrt mu, tau
## 2. get required link derivatives and tri indices.
## 3. transform derivs to derivs wrt eta (gamlss.etamu).
## 4. get the grad and Hessian etc for the model
##    via a call to gamlss.gH  
## the first derivatives of the log likelihood w.r.t
## the first and second parameters...
 
  if (K<1) stop("number of categories must be at least 2") 
   stats <- list()
  
  for (i in 1:K) {
    stats[[i]] <- make.link("identity")
    fam <- structure(list(link="identity",canonical="none",linkfun=stats[[i]]$linkfun,
           mu.eta=stats[[i]]$mu.eta),
           class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  } 

  residuals <- function(object,type=c("deviance")) {
  ## Deviance residuals where sign depends on whether classification correct (+ve)
  ## or not (-ve)...
      type <- match.arg(type)
      ## get category probabilities...
      p <- object$family$predict(object$family,eta=object$linear.predictors)[[1]]
      ## now get most probable category for each observation
      pc <- apply(p,1,function(x) which(max(x)==x)[1])-1 
      n <- length(pc)
      ## +ve sign if class correct, -ve otherwise
      sgn <- rep(-1,n); sgn[pc==object$y] <- 1
      ## now get the deviance...
      sgn*sqrt(-2*log(pmax(.Machine$double.eps,p[1:n + object$y*n]))) 
  } ## residuals

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family$data
  ## contains any family specific extra information. 
  ## if se = FALSE returns one item list containing matrix otherwise 
  ## list of two matrices "fit" and "se.fit"... 

    if (is.null(eta)) { 
      lpi <- attr(X,"lpi") 
      if (is.null(lpi)) {
        lpi <- list(1:ncol(X))
      } 
      K <- length(lpi) ## number of linear predictors
      eta <- matrix(0,nrow(X),K)
      if (se) { 
        ve <- matrix(0,nrow(X),K) ## variance of eta
        ce <- matrix(0,nrow(X),K*(K-1)/2) ## covariance of eta_i eta_j
      } 
      for (i in 1:K) { 
        Xi <- X[,lpi[[i]],drop=FALSE]
        eta[,i] <- Xi%*%beta[lpi[[i]]] ## ith linear predictor
        if (se) { ## variance and covariances for kth l.p.
          ve[,i] <- drop(pmax(0,rowSums((Xi%*%Vb[lpi[[i]],lpi[[i]]])*Xi)))
          ii <- 0
          if (i<K) for (j in (i+1):K) {
            ii <- ii + 1
            ce[,ii] <- drop(pmax(0,rowSums((Xi%*%Vb[lpi[[i]],lpi[[j]]])*X[,lpi[[j]]])))
          }
        }
      }
    } else { 
      se <- FALSE
    }
    gamma <- cbind(1,exp(eta))
    beta <- rowSums(gamma)
    gamma <- gamma/beta ## category probabilities
    vp <- gamma*0
    if (se) { ## need to loop to find se of probabilities...
      for (j in 1:(K+1)) {
        ## get dp_j/deta_k...
        if (j==1) dp <- -gamma[,-1,drop=FALSE]/beta else { 
          dp <- -gamma[,j]*gamma[,-1,drop=FALSE]
          dp[,j-1] <- gamma[,j]*(1-gamma[,j]) 
        }
        ## now compute variance... 
        vp[,j] <- rowSums(dp^2*ve)
        ii <- 0
        for (i in 1:K) if (i<K) for (k in (i+1):K) {
          ii <- ii + 1
          vp[,j] <- vp[,j] + 2 * dp[,i]*dp[,k]*ce[,ii] 
        }
        vp[,j] <- sqrt(pmax(0,vp[,j])) ## transform to se
      }
      return(list(fit=gamma,se.fit=vp))
    } ## if se
    list(fit=gamma)
  } ## multinom predict

  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    multinom <- list()
    object$y <- round(object$y) 
    multinom$nj <- tabulate(object$y+1) ## count each class membership
    multinom$n <- sum(multinom$nj) ## total number
    multinom$K <- length(multinom$nj)-1 ## number of linear predictors
    ## compute exp(eta) for categories 1..K
    multinom$gamma <- c(1,solve(diag(multinom$n/multinom$nj[-1],multinom$K)-
                        matrix(1,multinom$K,multinom$K),rep(1,multinom$K)))
    multinom$gamma <- log(multinom$gamma/sum(multinom$gamma))
    object$null.deviance <- -2*sum(multinom$gamma[object$y+1])
  })

  ll <- function(y,X,coef,wt,family,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL,eta=NULL) {
  ## Function defining the logistic multimomial model log lik. 
  ## Assumption is that coding runs from 0..K, with 0 class having no l.p.
  ## argument eta is for debugging only, and allows direct FD testing of the 
  ## derivatives w.r.t. eta. 
  ## ... this matches binary log reg case... 
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.
    n <- length(y)
    if (is.null(eta)) {
      return.l <- FALSE
      jj <- attr(X,"lpi") ## extract linear predictor index
      K <- length(jj) ## number of linear predictors 
      eta <- matrix(1,n,K+1) ## linear predictor matrix (dummy 1's in first column)
      for (i in 1:K) eta[,i+1] <- X[,jj[[i]],drop=FALSE]%*%coef[jj[[i]]]
    } else { l2 <- 0;K <- ncol(eta);eta <- cbind(1,eta); return.l <- TRUE}
 
    if (K!=family$nlp) stop("number of linear predictors doesn't match")
    y <- round(y) ## just in case
    if (min(y)<0||max(y)>K) stop("response not in 0 to number of predictors + 1")
    
    ee <- exp(eta[,-1,drop=FALSE])
    beta <- 1 + rowSums(ee); alpha <- log(beta)
    
    l0 <- eta[1:n+y*n] - alpha ## log likelihood
    l <- sum(l0)    

    l1 <- matrix(0,n,K) ## first deriv matrix
 
    if (deriv>0) {
      for (i in 1:K) l1[,i] <- ee[,i]/beta ## alpha1
    
      ## the second derivatives...
    
      l2 <- matrix(0,n,K*(K+1)/2)
      ii <- 0; b2 <- beta^2
      for (i in 1:K) for (j in i:K) {
        ii <- ii + 1 ## column index
        l2[,ii] <- if (i==j) -l1[,i] + ee[,i]^2/b2 else (ee[,i]*ee[,j])/b2
      }

      ## finish first derivatives...
      for (i in 1:K) l1[,i] <- as.numeric(y==i) - l1[,i] 

    } ## if (deriv>0)
 
    l3 <- l4 <- 0 ## defaults
    tri <- family$tri ## indices to facilitate access to earlier results
    
    if (deriv>1) { ## the third derivatives...
      l3 <- matrix(0,n,(K*(K+3)+2)*K/6)
      ii <- 0; b3 <- b2 * beta
      for (i in 1:K) for (j in i:K) for (k in j:K) {
        ii <- ii + 1 ## column index
        if (i==j&&j==k) { ## all same
           l3[,ii] <- l2[,tri$i2[i,i]] + 2*ee[,i]^2/b2 - 2*ee[,i]^3/b3
        } else if (i!=j&&j!=k&i!=k) { ## all different
           l3[,ii] <- -2*(ee[,i]*ee[,j]*ee[,k])/b3
        } else { ## two same one different
           kk <- if (i==j) k else j ## get indices for differing pair
           l3[,ii] <- l2[,tri$i2[i,kk]] - 2*(ee[,i]*ee[,j]*ee[,k])/b3
        }
      }
    } ## if (deriv>1)

    if (deriv>3) { ## the fourth derivatives...
      l4 <- matrix(0,n,(6+K*11+K^2*6+K^3)*K/24)
      ii <- 0; b4 <- b3 * beta
      for (i in 1:K) for (j in i:K) for (k in j:K) for (l in k:K) {
        ii <- ii + 1 ## column index
        uni <- unique(c(i,j,k,l));
        nun <- length(uni) ## number of unique indices
        if (nun==1) { ## all equal
          l4[,ii] <- l3[,tri$i3[i,i,i]] + 4*ee[,i]^2/b2 - 10*ee[,i]^3/b3 + 6*ee[,i]^4/b4
        } else if (nun==4) { ## all unequal
          l4[,ii] <- 6*ee[,i]*ee[,j]*ee[,k]*ee[,l]/b4
        } else if (nun==3) { ## 2 same 2 different
          l4[,ii] <- l3[,tri$i3[uni[1],uni[2],uni[3]]]  +6*ee[,i]*ee[,j]*ee[,k]*ee[,l]/b4
        } else if (sum(uni[1]==c(i,j,k,l))==2) { ## 2 unique (2 of each)
          l4[,ii] <- l3[,tri$i3[uni[1],uni[2],uni[2]]] - 2 * ee[,uni[1]]^2*ee[,uni[2]]/b3 + 
                     6*ee[,i]*ee[,j]*ee[,k]*ee[,l]/b4
        } else { ## 3 of one 1 of the other
          if (sum(uni[1]==c(i,j,k,l))==1) uni <- uni[2:1] ## first index is triple repeat index
          l4[,ii] <- l3[,tri$i3[uni[1],uni[1],uni[2]]] - 4 * ee[,uni[1]]^2*ee[,uni[2]]/b3 + 
                     6*ee[,i]*ee[,j]*ee[,k]*ee[,l]/b4
        }
      }
    } ## if deriv>3

    if (return.l) return(list(l=l0,l1=l1,l2=l2,l3=l3,l4=l4)) ## for testing...

    if (deriv) {
      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,l1,l2,tri$i2,l3=l3,i3=tri$i3,l4=l4,i4=tri$i4,
                      d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
    } else ret <- list()
    ret$l <- l; ret
  } ## end ll multinom

  rd <- function(mu,wt,scale) {
    ## simulate data given fitted linear predictor matrix in mu 
    p <- exp(cbind(0,mu))
    p <- p/rowSums(p)
    cp <- t(apply(p,1,cumsum))
    apply(cp,1,function(x) min(which(x>runif(1))))-1    
  } ## rd

  initialize <- expression({
  ## Binarize each category and lm on 6*y-3 by category.
 
      n <- rep(1, nobs)
      ## should E be used unscaled or not?..
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      if (is.null(start)) {
        jj <- attr(x,"lpi")
        start <- rep(0,ncol(x))
        for (k in 1:length(jj)) { ## loop over the linear predictors      
          yt1 <- 6*as.numeric(y==k)-3
          x1 <- x[,jj[[k]],drop=FALSE]
          e1 <- E[,jj[[k]],drop=FALSE] ## square root of total penalty
          if (use.unscaled) {
            qrx <- qr(rbind(x1,e1))
            x1 <- rbind(x1,e1)
            startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
            startji[!is.finite(startji)] <- 0       
          } else startji <- pen.reg(x1,e1,yt1)
          start[jj[[k]]] <- startji ## copy coefficients back into overall start coef vector
        } ## lp loop
      }
  }) ## initialize multinom

  structure(list(family="multinom",ll=ll,link=NULL,#paste(link),
    nlp=round(K),rd=rd,
    tri = trind.generator(K), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,predict=predict,
    linfo = stats, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    available.derivs = 2 ## can use full Newton here
    ),class = c("general.family","extended.family","family"))
} ## end multinom



pen.reg <- function(x,e,y) {
## get coefficients of penalized regression of y on matrix x
## where e is a square root penalty. Idea is to use e mainly for 
## regularization, so that edf is close to rank of x.  
  if (sum(abs(e))==0) { ## no penalization - easy
    b <- qr.coef(qr(x),y);b[!is.finite(b)] <- 0
    return(b)
  }
  ## need to adjust degree of penalization, so best to QR
  ## the x matrix up front...
  qrx <- qr(x,LAPACK=TRUE)
  R <- qr.R(qrx)
  r <- ncol(R)
  rr <- Rrank(R) ## rank of R/X
  R[,qrx$pivot] <- R ## unpivot
  Qy <- qr.qty(qrx,y)[1:ncol(R)]  
  ## now we want estimates with penalty weight low enough 
  ## EDF is k * rr where k is somewhere in e.g. (.7,.9)
  k <- .01 * norm(R)/norm(e)
  qrr <- qr(rbind(R,e*k));
  edf <- sum(qr.Q(qrr)[1:r,]^2) 
  while (edf > .9*rr) { ## increase penalization
    k <- k*10
    qrr <- qr(rbind(R,e*k));
    edf <- sum(qr.Q(qrr)[1:r,]^2)
  } 
  while (edf<.7*rr) { ## reduce penalization
    k <- k/20
    qrr <- qr(rbind(R,e*k));
    edf <- sum(qr.Q(qrr)[1:r,]^2)
  } 
  b <- qr.coef(qrr,c(Qy,rep(0,nrow(e))));b[!is.finite(b)] <- 0
  b
} ## pen.reg

## code for zero inflated Poisson models


#log1ex <- function(x) {
## evaluate log(1+exp(x)) accurately and avoiding overflow
#  y <- x
#  big <- -log(.Machine$double.eps)+5 ## exp(big) overwhelms 1
#  ind <- x > big
#  y[ind] <- x[ind] ## log(1+exp(x)) = x to machine precision
#  ## compute number below which log(1+exp(x)) = exp(x) to
#  ## machine precision... 
#  small <- log(sqrt(.Machine$double.eps))
#  ind1 <- x < small
#  y[ind1] <- exp(x[ind1])
#  ind <- !ind&!ind1 ## the moderate size elements 
#  y[ind] <- log(1+exp(x[ind]))
#  y 
#}

#logist <- function(x) {
## overflow proof logistic
#  ind <- x > 0; y <- x
#  y[ind] <- 1/(exp(-x[ind])+1) 
#  ex <- exp(x[!ind])
#  y[!ind] <- ex/(1+ex)
#  y
#}

l1ee <- function(x) {
## log(1-exp(-exp(x)))...
  ind <- x < log(.Machine$double.eps)/3
  ex <- exp(x);exi <- ex[ind]
  l <- log(1-exp(-ex))
  l[ind] <- log(exi-exi^2/2+exi^3/6)
  ind <- x < -log(.Machine$double.xmax)
  l[ind] <- x[ind]
  l
}

lee1 <- function(x) {
## log(exp(exp(x))-1)...
  ind <- x < log(.Machine$double.eps)/3
  ex <- exp(x);exi <- ex[ind]
  l <- log(exp(ex)-1)
  l[ind] <- log(exi+exi^2/2+exi^3/6)
  ind <- x < -log(.Machine$double.xmax)
  l[ind] <- x[ind]
  ind <- x > log(log(.Machine$double.xmax))
  l[ind] <- ex[ind]
  l
}

ldg <- function(g,deriv=4) {
  alpha <- function(g) {
    ind <- g > log(.Machine$double.eps)/3
    eg <- exp(g)
    g[ind] <- eg[ind]/(1-exp(-eg[ind]))
    g[!ind] <- 1+eg[!ind]/2 + eg[!ind]^2/12
    g
  }
  ind <- g < log(.Machine$double.eps)/3
  ghi <- log(log(.Machine$double.xmax)) + 1 
  ## ... above ghi alpha(g) is simply exp(g) 
  ii <- g>ghi 
  a <- alpha(g)
  eg <- exp(g)
  l2 <- a*(a-eg-1)
  egi <- eg[ind]
  ## in the lower tail alpha = 1 + b, where b = eg/2 + eg^2/12
  ## so l'' = alpha*(b-eg)...
  b <- egi*(1+egi/6)/2
  l2[ind] <- a[ind]*(b-egi)
  l2[ii] <- -exp(g[ii])
  l3 <- l4 <- NULL
  ## in a similar vein l3 can be robustified...  
  if (deriv>1) {
    l3 <- a*(a*(-2*a + 3*(eg+1)) - 3*eg - eg^2 - 1) 
    l3[ind] <- a[ind]*(-b-2*b^2+3*b*egi-egi^2)
    l3[ii] <- -exp(g[ii])
  }
  ## finally l4, which requires a similar approach...
  if (deriv>2) {
    l4 <- a*(6*a^3 - 12*(eg+1)*a^2+4*eg*a+7*(eg+1)^2*a-(4+3*eg)*eg -(eg+1)^3)
    l4[ind] <- a[ind]*(6*b*(3+3*b+b^2) - 12*egi*(1+2*b+b^2) - 12*b*(2-b) + 4*egi*(1+b)+
                     7*(egi^2+2*egi+b*egi^2+2*b*egi+b)-(4+3*egi)*egi-egi*(3+3*egi+egi^2))
   
    l4[ii] <- -exp(g[ii])
  }
  l1=-a
  ghi <- log(.Machine$double.xmax)/5
  ii <- g > ghi
  if (sum(ii)) {
    l1[ii] <- l2[ii] <- l3[ii] <- l4[ii] <- -exp(ghi)
  }
  list(l1=l1,l2=l2,l3=l3,l4=l4)
} ## ldg

lde <- function(eta,deriv=4) {
  ## llog lik derivs w.r.t. eta
  ind <- eta < log(.Machine$double.eps)/3
  ii <- eta > log(.Machine$double.xmax)
  l1 <- et <- exp(eta);eti <- et[ind]
  l1[!ind] <- et[!ind]/(exp(et[!ind])-1)
  b <- -eti*(1+eti/6)/2
  l1[ind] <- 1+b
  l1[ii] <- 0
  ## l2 ...   
  l2 <- l1*((1-et)-l1)
  l2[ind] <- -b*(1+eti+b) - eti
  l2[ii] <- 0
  l3 <- l4 <- NULL
  ## l3 ...
  if (deriv>1) {
    ii <- eta > log(.Machine$double.xmax)/2
    l3 <- l1*((1-et)^2-et - 3*(1-et)*l1 + 2*l1^2)
    l3[ind] <- l1[ind]*(-3*eti+eti^2 -3*(-eti+b-eti*b) + 2*b*(2+b))
    l3[ii] <- 0
  }
  ## l4 ...
  if (deriv>2) {
    ii <- eta > log(.Machine$double.xmax)/3
    l4 <- l1*((3*et-4)*et + 4*et*l1 + (1-et)^3 - 7*(1-et)^2*l1 + 12*(1-et)*l1^2 
            - 6*l1^3)
    l4[ii] <- 0
    l4[ind] <- l1[ind]*(4*l1[ind]*eti - eti^3 - b -7*b*eti^2 - eti^2 - 5*eti -
                 10*b*eti - 12*eti*b^2 - 6*b^2 - 6*b^3)
  }
  list(l1=l1,l2=l2,l3=l3,l4=l4)
} ## lde


zipll <- function(y,g,eta,deriv=0) {
## function to evaluate zero inflated Poisson log likelihood
## and its derivatives w.r.t. g/gamma and eta where 
## 1-p = exp(-exp(eta)) and lambda = exp(gamma), for each datum in vector y.
## p is probability of potential presence. lambda is Poisson mean
## given potential presence. 
## deriv: 0 - eval
##        1 - grad (l,p) and Hess (ll,lp,pp)
##        2 - third derivs lll,llp,lpp,ppp
##        4 - 4th derivs. llll,lllp,llpp,lppp,pppp

   l1 <- El2 <- l2 <- l3 <- l4 <- NULL
   zind <- y == 0 ## the index of the zeroes
   ## yz <- y[zind];
   yp <- y[!zind]
   l <- et <- exp(eta)
   l[zind] <- -et[zind] # -exp(eta[ind])
   l[!zind] <- l1ee(eta[!zind]) + yp*g[!zind] - lee1(g[!zind]) - lgamma(yp+1)
   p <- 1-exp(-et) ## probablity of non-zero

   if (deriv>0) { ## get first and second derivs...
     n <- length(y)
     l1 <- matrix(0,n,2)
     le <- lde(eta,deriv) ## derivs of ll wrt eta     
     lg <- ldg(g,deriv) ## derivs of ll wrt gamma
     l1[!zind,1] <- yp + lg$l1[!zind]  ## l_gamma, y>0
     l1[zind,2] <- l[zind] ## l_eta, y==0
     l1[!zind,2] <- le$l1[!zind]  ## l_eta, y>0    

     El2 <- l2 <- matrix(0,n,3)
     ## order gg, ge, ee... 
     l2[!zind,1] <- lg$l2[!zind]   ## l_gg, y>0
     l2[!zind,3] <- le$l2[!zind]   ## l_ee, y>0
     l2[zind,3]  <- l[zind] ## l_ee, y=0
     El2[,1] <- p*lg$l2             ## E(l_gg)
     El2[,3] <- -(1-p)*et + p*le$l2 ## E(l_ee)
   }
   if (deriv>1) {
      ## the third derivatives
      ## order ggg,gge,gee,eee
      l3 <- matrix(0,n,4) 
      l3[!zind,1] <- lg$l3[!zind]   ## l_ggg, y>0
      l3[!zind,4] <- le$l3[!zind]   ## l_eee, y>0
      l3[zind,4]  <- l[zind]        ## l_eee, y=0
   }
   if (deriv>3) {
      ## the fourth derivatives
      ## order gggg,ggge,ggee,geee,eeee
      l4 <- matrix(0,n,5) 
      l4[!zind,1] <- lg$l4[!zind]   ## l_gggg, y>0
      l4[!zind,5] <- le$l4[!zind]   ## l_eeee, y>0
      l4[zind,5]  <- l[zind]        ## l_eeee, y=0
   }
   list(l=l,l1=l1,l2=l2,l3=l3,l4=l4,El2=El2)
} ## zipll


ziplss <-  function(link=list("identity","identity")) {
## Extended family for Zero Inflated Poisson fitted as gamlss 
## type model.
## mu1 is Poisson mean, while mu2 is zero inflation parameter.
  ## first deal with links and their derivatives...
  if (length(link)!=2) stop("ziplss requires 2 links specified as character strings")
  okLinks <- list(c("identity"),c("identity"))
  stats <- list()
  param.names <- c("Poisson mean","binary probability")
  for (i in 1:2) {
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
    stop(link[[i]]," link not available for ",param.names[i]," parameter of ziplss")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
           mu.eta=stats[[i]]$mu.eta),
           class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  } 

  residuals <- function(object,type=c("deviance","response")) {
      ls <- function(y) {
        ## compute saturated likelihood for ziplss model 
        l <- y;l[y<2] <- 0
        ind <- y > 1 & y < 18
        ## lambda maximizing likelihood for y = 2 to 17 
        glo <- c(1.593624,2.821439,3.920690,4.965114,5.984901,6.993576,
                 7.997309,8.998888,9.999546,10.999816,11.999926,12.999971,
                 13.999988,14.999995,15.999998,16.999999)
        g <- y ## maximizing lambda essentially y above this
        g[ind] <- glo[y[ind]-1]
        ind <- y > 1
        l[ind] <- zipll(y[ind],log(g[ind]),g[ind]*0+1e10,deriv=0)$l
        l
      } ## ls

      type <- match.arg(type)
      p <- exp(-exp(object$fitted[,2]));
      lam <- exp(object$fitted[,1])
      ind <- lam > .Machine$double.eps^.5
      ## compute E(y)
      Ey <- p ## very small lambda causes conditional expectation to be 1
      Ey[ind] <- p[ind]*lam[ind]/(1-exp(-lam[ind])) 
      rsd <- object$y - Ey ## raw residuals
      if (type=="response") return(rsd)
      else { ## compute deviance residuals
        sgn <- sign(rsd)
        ind <- object$y == 0
        rsd <- pmax(0,2*(ls(object$y) - zipll(object$y,object$fitted[,1],object$fitted[,2],deriv=0)$l))
        rsd <- sqrt(rsd)*sgn
      }
      rsd
  } ## residuals

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family$data
  ## contains any family specific extra information. 
  ## if se = FALSE returns one item list containing matrix otherwise 
  ## list of two matrices "fit" and "se.fit"... 

    if (is.null(eta)) { 
      lpi <- attr(X,"lpi") 
      X1 <- X[,lpi[[1]],drop=FALSE]
      X2 <- X[,lpi[[2]],drop=FALSE]
      gamma <- drop(X1%*%beta[lpi[[1]]]) ## linear predictor for poisson parameter 
      eta <- drop(X2%*%beta[lpi[[2]]]) ## linear predictor for presence parameter 
      if (se) {
        v.g <- drop(pmax(0,rowSums((X1%*%Vb[lpi[[1]],lpi[[1]]])*X1))) ## var of gamma
        v.e <- drop(pmax(0,rowSums((X1%*%Vb[lpi[[1]],lpi[[1]]])*X1))) ## var of eta
        v.eg <- drop(pmax(0,rowSums((X1%*%Vb[lpi[[1]],lpi[[2]]])*X2))) ## cov of eta, gamma
      }
    } else { 
      se <- FALSE
      gamma <- eta[,1]
      eta <- eta[,2]
    }
    et <- exp(eta)
    mu <- p <- 1 - exp(-et)
    fv <- lambda <- exp(gamma)  
    ind <- gamma < log(.Machine$double.eps)/2
    mu[!ind] <- lambda[!ind]/(1-exp(-lambda[!ind]))
    mu[ind] <- 1
    fv <- list(p*mu)    ## E(y)    
    if (!se) return(fv) else {
      df.de <- p  
      ind <- eta < log(.Machine$double.xmax)/2
      df.de[!ind] <- 0
      df.de[ind] <- exp(-et[ind])*et[ind]
      df.de <- df.de * mu
      df.dg <- ((lambda + 1)*mu - mu^2)*p
      fv[[2]] <- sqrt(df.dg^2*v.g + df.de^2*v.e + 2 * df.de * df.dg * v.eg)   
      names(fv) <- c("fit","se.fit")
      return(fv)
    }
  } ## predict


  rd <- function(mu,wt,scale) {
    ## simulate data given fitted latent variable in mu 
      rzip <- function(gamma,eta) { ## generate ziP deviates according to model and lp gamma
        y <- gamma; n <- length(y)
        lambda <- exp(gamma)
        p <- 1- exp(-exp(eta))
        ind <- p > runif(n)
        y[!ind] <- 0
        np <- sum(ind)
        ## generate from zero truncated Poisson, given presence...
        y[ind] <- qpois(runif(np,dpois(0,lambda[ind]),1),lambda[ind])
        y
      } 
      rzip(mu[,1],mu[,2])
  } ## rd

  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    ## null model really has two parameters... probably need to newton iterate
    ls <- function(y) {
        ## compute saturated likelihood for ziplss model 
        l <- y;l[y<2] <- 0
        ind <- y > 1 & y < 18
        ## lambda maximizing likelihood for y = 2 to 17 
        glo <- c(1.593624,2.821439,3.920690,4.965114,5.984901,6.993576,
                 7.997309,8.998888,9.999546,10.999816,11.999926,12.999971,
                 13.999988,14.999995,15.999998,16.999999)
        g <- y ## maximizing lambda essentially y above this
        g[ind] <- glo[y[ind]-1]
        ind <- y > 1
        l[ind] <- zipll(y[ind],log(g[ind]),g[ind]*0+1e10,deriv=0)$l
        l
    } ## ls

    fp <- function(p,y) {
    ## compute zero related part of log likelihood
       eps <- .Machine$double.eps^.5 
       l1p <- if (p>eps) log(1-p) else -p - p^2/2 
       l1p*sum(y==0) + log(p)*sum(y>0)
    } ## fp

    flam <- function(lam,y) {
      ## compute >0 part of log likelihood
      y <- y[y>0]
      sum(y*log(lam) - log(exp(lam)-1) - lgamma(y+1))
    } ## flam

    ## optimize zero repated part of likelihood w.r.t. p...
    lnull <- optimize(fp,interval=c(1e-60,1-1e-10),y=object$y,maximum=TRUE)$objective
    ## optimize >0 part for lambda...
    my <- mean(object$y[object$y>0])
    lnull <- lnull + optimize(flam,interval=c(my/2,my*2),y=object$y,maximum=TRUE)$objective
    object$null.deviance <- 2*(sum(ls(object$y)) - lnull)
   
  }) ## postproc


  ll <- function(y,X,coef,wt,family,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
  ## function defining the gamlss ZIP model log lik. 
  ## First l.p. defines Poisson mean, given presence (lambda)
  ## Second l.p. defines probability of presence (p)
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.
    jj <- attr(X,"lpi") ## extract linear predictor index
    eta <- X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]]
    lambda <- family$linfo[[1]]$linkinv(eta)
    eta1 <- X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]] 
    p <-  family$linfo[[2]]$linkinv(eta1) 
    
    ##n <- length(y)
    ## l1 <- matrix(0,n,2)
    zl <- zipll(y,lambda,p,deriv)

    if (deriv>0) {
      ## need some link derivatives for derivative transform
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(lambda),family$linfo[[2]]$d2link(p))
    }
 
    ## l3 <- l4 <- 
    g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
      ## the third derivatives
      ## order lll,llp,lpp,ppp
      g3 <- cbind(family$linfo[[1]]$d3link(lambda),family$linfo[[2]]$d3link(p))
    }

    if (deriv>3) {
      ## the fourth derivatives
      ## order llll,lllp,llpp,lppp,pppp
      g4 <- cbind(family$linfo[[1]]$d4link(lambda),family$linfo[[2]]$d4link(p))
    }
    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3
      i4 <- family$tri$i4
   
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- gamlss.etamu(zl$l1,zl$l2,zl$l3,zl$l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)

      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                      d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
    } else ret <- list()
    ret$l <- sum(zl$l); ret
  } ## end ll for ZIP

  initialize <- expression({ ## for ZIP
  ## Idea is to regress binarized y on model matrix for p. 
  ## Then downweight any y=0 with p<0.5 and regress g(y) on 
  ## the model matrix for lambda - don't drop as this may
  ## induce rank deficiency in model matrix! 
  ## May be called in both gam.fit5 and initial.spg... 
  ## note that appropriate E scaling
  ## for full calculation may be inappropriate for initialization 
  ## which is basically penalizing something different here.
  ## best we can do here is to use E only as a regularizer.
      n <- rep(1, nobs)
      if (all.equal(y,round(y))!=TRUE) {
          stop("Non-integer response variables are not allowed with ziplss ")
      }
      if ((min(y)==0&&max(y)==1)) stop("Using ziplss for binary data makes no sense")
      ## should E be used unscaled or not?..
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      if (is.null(start)) {
        jj <- attr(x,"lpi")
        start <- rep(0,ncol(x))
        x1 <- x[,jj[[2]],drop=FALSE]
        e1 <- E[,jj[[2]],drop=FALSE] ## square root of total penalty
        yt1 <- as.numeric(as.logical(y)) ## binarized response
        if (use.unscaled) {
          qrx <- qr(rbind(x1,e1))
          x1 <- rbind(x1,e1)
          startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
          startji[!is.finite(startji)] <- 0       
        } else startji <- pen.reg(x1,e1,yt1)
        start[jj[[2]]] <- startji
        p <- drop(x1[1:nobs,,drop=FALSE] %*% startji) ## probability of presence
        ind <- y==0 & p < 0.5 ## downweight these for estimating lambda
        w <- rep(1,nobs); w[ind] <- .1

        ## note assumption that working scale is log...
        yt1 <-  family$linfo[[1]]$linkfun(log(abs(y)+(y==0)*.2))

        yt1 <- yt1*w        

        x1 <-  w*x[,jj[[1]],drop=FALSE];e1 <- E[,jj[[1]],drop=FALSE]
        if (use.unscaled) {
          x1 <- rbind(x1,e1)
          startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))   
          startji[!is.finite(startji)] <- 0
        } else startji <- pen.reg(x1,e1,yt1)
        start[jj[[1]]] <- startji
      }
  }) ## initialize ziplss

  structure(list(family="ziplss",ll=ll,link=paste(link),nlp=2,
    tri = trind.generator(2), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,rd=rd,predict=predict,
    linfo = stats, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    available.derivs = 2 ## can use full Newton here
    ),class = c("general.family","extended.family","family"))
} ## ziplss





