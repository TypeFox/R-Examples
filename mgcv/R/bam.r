## routines for very large dataset generalized additive modelling.
## (c) Simon N. Wood 2009-2015


ls.size <- function(x) {
## If `x' is a list, return the size of its elements, in bytes, in a named array
## otherwise return the size of the object
 if (is.list(x)==FALSE) return(object.size(x))

 xn <- names(x)
 n <- length(x)
 sz <- rep(-1,n)
 for (i in 1:n) sz[i] <- object.size(x[[i]])
 names(sz) <- xn
 sz
} ## ls.size

rwMatrix <- function(stop,row,weight,X,trans=FALSE) {
## Routine to recombine the rows of a matrix X according to info in 
## stop, row and weight. Consider the ith row of the output matrix 
## ind <- 1:stop[i] if i==1 and ind <- (stop[i-1]+1):stop[i]
## otherwise. The ith output row is then X[row[ind],]*weight[ind]
  if (is.matrix(X)) { n <- nrow(X);p<-ncol(X);ok <- TRUE} else { n<- length(X);p<-1;ok<-FALSE}
  stop <- stop - 1;row <- row - 1 ## R indices -> C indices
  oo <-.C(C_rwMatrix,as.integer(stop),as.integer(row),as.double(weight),X=as.double(X),
          as.integer(n),as.integer(p),trans=as.integer(trans),work=as.double(rep(0,n*p)))
  if (ok) return(matrix(oo$X,n,p)) else
  return(oo$X) 
} ## rwMatrix

chol2qr <- function(XX,Xy,nt=1) {
## takes X'X and X'y and returns R and f
## equivalent to qr update.  
  op <- options(warn = -1) ## otherwise warns if +ve semidef
  R <- if (nt) pchol(XX,nt=nt) else chol(XX,pivot=TRUE)
  options(op)
  p <- length(Xy)
  ipiv <- piv <- attr(R,"pivot");ipiv[piv] <- 1:p
  rank <- attr(R,"rank");ind <- 1:rank
  if (rank<p) R[(rank+1):p,] <- 0 ## chol is buggy (R 3.1.0) with pivot=TRUE
  f <- c(forwardsolve(t(R[ind,ind]),Xy[piv[ind]]),rep(0,p-rank))[ipiv]
  R <- R[ipiv,ipiv]
  list(R=R,f=f)
} ## chol2qr

qr.update <- function(Xn,yn,R=NULL,f=rep(0,0),y.norm2=0,use.chol=FALSE,nt=1)
## Let X = QR and f = Q'y. This routine updates f and R
## when Xn is appended to X and yn appended to y. If R is NULL
## then initial QR of Xn is performed. ||y||^2 is accumulated as well.
## if use.chol==TRUE then quicker but less stable accumulation of X'X and
## X'y are used. Results then need post processing, to get R =chol(X'X)
## and f= R^{-1} X'y.
## if nt>1 and use.chol=FALSE then parallel QR is used 
{ p <- ncol(Xn)  
  y.norm2 <- y.norm2+sum(yn*yn)
  if (use.chol) { 
    if (is.null(R)) { 
      R <- crossprod(Xn)
      fn <- as.numeric(t(Xn)%*%yn) 
    } else {
      R <- R + crossprod(Xn)
      fn <- f + as.numeric(t(Xn)%*%yn)
    } 
    return(list(R=R,f=fn,y.norm2=y.norm2))
  } else { ## QR update
    if (!is.null(R)) {
      Xn <- rbind(R,Xn)
      yn <- c(f,yn)
    }
    qrx <- if (nt==1) qr(Xn,tol=0,LAPACK=TRUE) else pqr2(Xn,nt)
    fn <- qr.qty(qrx,yn)[1:p]
    rp <- qrx$pivot;rp[rp] <- 1:p # reverse pivot
    return(list(R = qr.R(qrx)[,rp],f=fn,y.norm2=y.norm2))
  }
} ## qr.update


qr.up <- function(arg) {
## routine for parallel computation of the QR factorization of 
## a large gam model matrix, suitable for calling with parLapply.
  wt <- rep(0,0) 
  dev <- 0    
  for (b in 1:arg$n.block) {
    ind <- arg$start[b]:arg$stop[b]
    X <- predict(arg$G,newdata=arg$mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,block.size=length(ind))
    rownames(X) <- NULL
    if (is.null(arg$coef)) eta1 <- arg$eta[ind] else eta1 <- drop(X%*%arg$coef) + arg$offset[ind]
    mu <- arg$linkinv(eta1) 
    y <- arg$G$y[ind] ## arg$G$model[[arg$response]] 
    weights <- arg$G$w[ind]
    mu.eta.val <- arg$mu.eta(eta1)
    good <- (weights > 0) & (mu.eta.val != 0)
    z <- (eta1 - arg$offset[ind])[good] + (y - mu)[good]/mu.eta.val[good]
    w <- (weights[good] * mu.eta.val[good]^2)/arg$variance(mu)[good]
    dev <- dev + sum(arg$dev.resids(y,mu,weights))
    wt <- c(wt,w)
    w <- sqrt(w)
    ## note assumption that nt=1 in following qr.update - i.e. each cluster node is strictly serial
    if (b == 1) qrx <- qr.update(w*X[good,],w*z,use.chol=arg$use.chol) 
    else qrx <- qr.update(w*X[good,],w*z,qrx$R,qrx$f,qrx$y.norm2,use.chol=arg$use.chol)
    rm(X);if(arg$gc.level>1) gc() ## X can be large: remove and reclaim
  }
  qrx$dev <- dev;qrx$wt <- wt
  if (arg$gc.level>1) { rm(arg,ind,mu,y,weights,mu.eta.val,good,z,w,wt,w);gc()}
  qrx
} ## qr.up

compress.df <- function(dat,m=NULL) {
## Takes dataframe in dat and compresses it by rounding and duplicate 
## removal. For metric variables we first find the unique cases.
## If there are <= m of these then these are employed, otherwise 
## rounding is used. Factors are always reduced to the number of 
## levels present in the data. Idea is that this function is called 
## with columns of dataframes corresponding to single smooths or marginals. 
  d <- ncol(dat) ## number of variables to deal with
  n <- nrow(dat) ## number of data/cases
  if (is.null(m)) m <- if (d==1) 1000 else if (d==2) 100 else 25 else
  if (d>1) m <- round(m^{1/d}) + 1
  
  mf <- mm <- 1 ## total grid points for factor and metric
  for (i in 1:d) if (is.factor(dat[,i])) {  
    mf <- mf * length(unique(as.vector(dat[,i]))) 
  } else {
    mm <- mm * m 
  } 
  if (is.matrix(dat[[1]])) { ## must replace matrix terms with vec(dat[[i]])
    dat0 <- data.frame(as.vector(dat[[1]]))
    if (d>1) for (i in 2:d) dat0[[i]] <- as.vector(dat[[i]])
    names(dat0) <- names(dat)
    dat <- dat0;rm(dat0)
  }
  xu <- uniquecombs(dat)
  if (nrow(xu)>mm*mf) { ## too many unique rows to use only unique
    for (i in 1:d) if (!is.factor(dat[,i])) { ## round the metric variables
      xl <- range(dat[,i])
      xu <- seq(xl[1],xl[2],length=m)
      dx <- xu[2]-xu[1]
      kx <- round((dat[,i]-xl[1])/dx)+1
      dat[,i] <- xu[kx] ## rounding the metric variables
    }
    xu <- uniquecombs(dat)
  }  
  k <- attr(xu,"index")
  ## shuffle rows in order to avoid induced dependencies between discretized
  ## covariates (which can mess up gam.side)...
  seed <- try(get(".Random.seed",envir=.GlobalEnv),silent=TRUE) ## store RNG seed
  if (inherits(seed,"try-error")) {
     runif(1)
     seed <- get(".Random.seed",envir=.GlobalEnv)
  }
  kind <- RNGkind(NULL)
  RNGkind("default","default")
  set.seed(1) ## ensure repeatability
  
  ii <- sample(1:nrow(xu),nrow(xu),replace=FALSE) ## shuffling index
  
  RNGkind(kind[1],kind[2])
  assign(".Random.seed",seed,envir=.GlobalEnv) ## RNG behaves as if it had not been used
  
  xu[ii,] <- xu  ## shuffle rows of xu
  k <- ii[k]     ## correct k index accordingly
  ## ... finished shuffle
  ## if arguments were matrices, then return matrix index
  if (length(k)>n) k <- matrix(k,nrow=n) 
  k -> attr(xu,"index")
  xu
} ## compress.df

check.term <- function(term,rec) {
## utility function for discrete.mf. Checks whether variables in "term"
## have already been discretized, and if so whether this discretization 
## can be re-used for the current "term". Stops if term already discretized
## but we can't re-use discretization. Otherwise returns index of k index
## or 0 if the term is not in the existing list.
  ii <- which(rec$vnames%in%term)
  if (length(ii)) { ## at least one variable already discretized
    if (length(term)==rec$d[min(ii)]) { ## dimensions match previous discretization
      if (sum(!(term%in%rec$vnames[ii]))) ("bam can not discretize with this nesting structure")
      else return(rec$ki[min(ii)]) ## all names match previous - retun index of previous
    } else stop("bam can not discretize with this nesting structure")
  } else return(0) ## no match
} ## check.term

discrete.mf <- function(gp,mf,pmf,m=NULL,full=TRUE) {
## discretize the covariates for the terms specified in smooth.spec
## id and factor by not allowed. pmf is a model frame for just the 
## parametric terms --- mini.mf is applied to this.
## if full is FALSE then parametric and response terms are ignored
## and what is returned is a list where columns can be of 
## different lengths.
## On exit... 
## * mf is a model frame containing the unique discretized covariate
##   values, in randomized order, padded to all be same length
## * nr records the number of unique discretized covariate values
##   i.e. the number of rows before the padding starts
## * k.start contains the starting column in index vector k, for
##   each variable.
## * k is the index matrix. The ith record of the 1st column of the 
##   jth variable is in row k[i,k.start[j]] of the corresponding 
##   column of mf.
## ... there is an element of nr and k.start for each variable of 
## each smooth, but varaibles are onlt discretized and stored in mf
## once. If there are no matrix variables then k.start = 1:(ncol(k)+1) 

  mf0 <- list()
  nk <- 0 ## count number of index vectors to avoid too much use of cbind
  for (i in 1:length(gp$smooth.spec)) nk <- nk + as.numeric(gp$smooth.spec[[i]]$by!="NA") +
    if (inherits(gp$smooth.spec[[i]],"tensor.smooth.spec")) length(gp$smooth.spec[[i]]$margin) else 1
  k <- matrix(0,nrow(mf),nk) ## each column is an index vector
  k.start <- 1:(nk+1) ## record last column for each term
  ik <- 0 ## index counter
  nr <- rep(0,nk) ## number of rows for term
  ## structure to record terms already processed...
  rec <- list(vnames = rep("",0), ## variable names
              ki = rep(0,0),      ## index of original index vector var relates to  
              d = rep(0,0))       ## dimension of terms involving this var
  ## loop through the terms discretizing the covariates...
  for (i in 1:length(gp$smooth.spec)) {
    nmarg <- if (inherits(gp$smooth.spec[[i]],"tensor.smooth.spec")) length(gp$smooth.spec[[i]]$margin) else 1
    maxj <- if (gp$smooth.spec[[i]]$by=="NA") nmarg else nmarg + 1 
    mi <- if (is.null(m)||length(m)==1) m else m[i]
    if (!is.null(gp$smooth.spec[[i]]$id)) stop("discretization can not handle smooth ids")
    j <- 0
    for (jj in 1:maxj) { ## loop through marginals
      if (jj==1&&maxj!=nmarg) termi <- gp$smooth.spec[[i]]$by else {
        j <- j + 1
        termi <- if (inherits(gp$smooth.spec[[i]],"tensor.smooth.spec")) gp$smooth.spec[[i]]$margin[[j]]$term else 
                 gp$smooth.spec[[i]]$term          
      } 
      ik.prev <- check.term(termi,rec) ## term already discretized?
      ik <- ik + 1 ## increment index counter
      if (ik.prev==0) { ## new discretization required
         mfd <- compress.df(mf[termi],m=mi)
         ki <- attr(mfd,"index")
         if (is.matrix(ki)) {
           ind <- (ik+1):length(k.start) 
           k.start[ind] <- k.start[ind] + ncol(ki)-1    ## adjust start indices
           k <- cbind(k,matrix(0,nrow(k),ncol(ki)-1)) ## extend index matrix
           ind <- k.start[ik]:(k.start[ik+1]-1) 
           k[,ind] <- ki 
         } else {
           k[,k.start[ik]] <- ki
         }
         nr[ik] <- nrow(mfd)
         mf0 <- c(mf0,mfd) 
         ## record variable discretization info...
         d <- length(termi)
         rec$vnames <- c(rec$vnames,termi)
         rec$ki <- c(rec$ki,rep(ik,d))
         rec$d <- c(rec$d,rep(d,d))
       } else { ## re-use an earlier discretization...
         ind.prev <- k.start[ik.prev]:(k.start[ik.prev+1]-1)
         ind <- (ik+1):length(k.start)
         k.start[ind] <- k.start[ind] + length(ind.prev)-1
         ind <- k.start[ik]:(k.start[ik+1]-1)
         k[,ind] <- k[,ind.prev]
         #k[,ik] <- k[,ik.prev]
         nr[ik] <- nr[ik.prev]
       }
    } ## end marginal jj loop
  } ## term loop (i)

## old long winded code...
## deal with any by variable (should always go first as basically a 1D marginal)... 
#    if (gp$smooth.spec[[i]]$by!="NA") {
#      termi <- gp$smooth.spec[[i]]$by ## var name
#      ik.prev <- check.term(termi,rec) ## term already discretized?
#      ik <- ik + 1 ## increment index counter
#      if (ik.prev==0) { ## new discretization required
#         mfd <- compress.df(mf[termi],m=mi)
#         k[,ik] <- attr(mfd,"index")
#         nr[ik] <- nrow(mfd)
#         mf0 <- c(mf0,mfd) 
#         ## record variable discretization info... 
#         rec$vnames <- c(rec$vnames,termi)
#         rec$ki <- c(rec$ki,ik)
#         rec$d <- c(rec$d,1)
#        
#       } else { ## re-use an earlier discretization...
#         k[,ik] <- k[,ik.prev]
#         nr[ik] <- nr[ik.prev]
#       }
#    }  
#    if (inherits(gp$smooth.spec[[i]],"tensor.smooth.spec")) { ## tensor branch
#      for (j in 1:length(gp$smooth.spec[[i]]$margin)) { ## loop through margins
#        termj <- gp$smooth.spec[[i]]$margin[[j]]$term
#        ik.prev <- check.term(termj,rec) 
#        ik <- ik + 1
#        if (ik.prev==0) { ## new discretization required
#          mfd <- compress.df(mf[termj],m=mi)
#          k[,ik] <- attr(mfd,"index")
#          nr[ik] <- nrow(mfd)
#          mf0 <- c(mf0,mfd)
#          ## record details...
#          d <- length(termj)
#          rec$vnames <- c(rec$vnames,termj)
#          rec$ki <- c(rec$ki,rep(ik,d))
#          rec$d <- c(rec$d,rep(d,d))
#        } else { ## re-use an earlier discretization... 
#          k[,ik] <- k[,ik.prev]
#          nr[ik] <- nr[ik.prev]
#        }
#      }
#    } else { ## not te or ti...
#      termi <- gp$smooth.spec[[i]]$term 
#      ik.prev <- check.term(termi,rec) 
#      ik <- ik + 1 ## index counter
#      if (ik.prev==0) { ## new discretization required
#        mfd <- compress.df(mf[termi],m=mi)
 #       k[,ik] <- attr(mfd,"index")
#        nr[ik] <- nrow(mfd)
#        mf0 <- c(mf0,mfd) 
#        d <- length(termi)
#        rec$vnames <- c(rec$vnames,termi)
#        rec$ki <- c(rec$ki,rep(ik,d))
#        rec$d <- c(rec$d,rep(d,d))
#      } else { ## re-use an earlier discretization...
#        k[,ik] <- k[,ik.prev]
#        nr[ik] <- nr[ik.prev]
#      }
#    }
#   
#  } ## main term loop

  ## obtain parametric terms and..
  ## pad mf0 so that all rows are the same length
  ## padding is necessary if gam.setup is to be used for setup

  if (full) {
    maxr <- max(nr) 
    pmf0 <- mini.mf(pmf,maxr) ## deal with parametric components
    if (nrow(pmf0)>maxr) maxr <- nrow(pmf0)
    mf0 <- c(mf0,pmf0) ## add parametric terms to end of mf0
    seed <- try(get(".Random.seed",envir=.GlobalEnv),silent=TRUE) ## store RNG seed
    if (inherits(seed,"try-error")) {
       runif(1)
       seed <- get(".Random.seed",envir=.GlobalEnv)
    }
    kind <- RNGkind(NULL)
    RNGkind("default", "default")
    set.seed(9)  
    for (i in 1:length(mf0)) {
      me <- length(mf0[[i]]) 
      if (me < maxr) mf0[[i]][(me+1):maxr] <- sample(mf0[[i]],maxr-me,replace=TRUE)
    }
    ## add response so that gam.setup can do its thing... 
  
    mf0[[gp$response]] <- sample(mf[[gp$response]],maxr,replace=TRUE)
    
    ## mf0 is the discretized model frame (actually a list), padded to have equal length rows
    ## k is the index vector for each sub-matrix, only the first nr rows of which are
    ## to be retained... Use of check.names=FALSE ensures, e.g. 'offset(x)' not changed...

    ## now copy back into mf so terms unchanged
    #mf <- mf[1:maxr,]
    mf <- mf[sample(1:nrow(mf),maxr,replace=TRUE),]
    for (na in names(mf0)) mf[[na]] <- mf0[[na]] 
    RNGkind(kind[1], kind[2])
    assign(".Random.seed", seed, envir = .GlobalEnv)
  } else mf <- mf0

  ## finally one more pass through, expanding k, k.start and nr to deal with replication that
  ## will occur with factor by variables...
  ik <- ncol(k)+1 ## starting index col for this term in k.start
  for (i in length(gp$smooth.spec):1) { ## work down through terms so insertion painless
    if (inherits(gp$smooth.spec[[i]],"tensor.smooth.spec")) nd <-  
         length(gp$smooth.spec[[i]]$margin) else nd <- 1 ## number of indices
    ik <- ik - nd ## starting index if no by  
    if (gp$smooth.spec[[i]]$by!="NA") {
      ik <- ik - 1 ## first index
      nd <- nd + 1 ## number of indices
      byvar <- mf[[gp$smooth.spec[[i]]$by]]
      if (is.factor(byvar)) { ## then need to expand nr and index matrix
        nex <- length(levels(byvar))  ## number of copies of term indices
        if (is.ordered(byvar)) nex <- nex - 1 ## first level dropped
        if (nex>0) { ## insert index copies
          ii0 <- if (ik>1) 1:(ik-1) else rep(0,0) ## earlier
          ii1 <- if (ik+nd-1 < length(nr)) (ik+nd):length(nr) else rep(0,0) ## later
          ii <- ik:(ik+nd-1) ## cols for this term    
          ## indices for columns of k... 
          kk0 <- if (ik>1) 1:(k.start[ik]-1) else rep(0,0) ## earlier
          kk1 <- if (ik+nd-1 < length(nr)) k.start[ik+nd]:ncol(k) else rep(0,0) ## later
          kk <- k.start[ik]:(k.start[ik+nd]-1) ## cols for this term
          k <- cbind(k[,kk0,drop=FALSE],k[,rep(kk,nex),drop=FALSE],k[,kk1,drop=FALSE])
          nr <- c(nr[ii0],rep(nr[ii],nex),nr[ii1])
          ## expand k.start...
          nkk <- length(kk) ## number of k columns in term to be repeated
          k.start <- c(k.start[ii0],rep(k.start[ii],nex)+rep(0:(nex-1),each=nkk)*nkk,
                       (nex-1)*nkk+c(k.start[ii1],k.start[length(k.start)]))
        }
      } ## factor by 
    } ## existing by
  } ## smooth.spec loop
  list(mf=mf,k=k,nr=nr,k.start=k.start)
} ## discrete.mf

mini.mf <-function(mf,chunk.size) {
## takes a model frame and produces a representative subset of it, suitable for 
## basis setup.
  ## first count the minimum number of rows required for representiveness
  mn <- 0
  for (j in 1:length(mf)) mn <- mn + if (is.factor(mf[[j]])) length(levels(mf[[j]])) else 2
  if (chunk.size < mn) chunk.size <- mn
  n <- nrow(mf)
  if (n <= chunk.size) return(mf)

  seed <- try(get(".Random.seed",envir=.GlobalEnv),silent=TRUE) ## store RNG seed
  if (inherits(seed,"try-error")) {
     runif(1)
     seed <- get(".Random.seed",envir=.GlobalEnv)
  }
  kind <- RNGkind(NULL)
  RNGkind("default", "default")
  set.seed(66)  
  ## randomly sample from original frame...
  ind <- sample(1:n,chunk.size)
  mf0 <- mf[ind,,drop=FALSE]
  ## ... now need to ensure certain sorts of representativeness

  ## work through elements collecting the rows containing 
  ## max and min for each variable, and a random row for each 
  ## factor level....

  ind <- sample(1:n,n,replace=FALSE) ## randomized index for stratified sampling w.r.t. factor levels
  fun <- function(X,fac,ind) ind[fac[ind]==X][1] ## stratified sampler
  k <- 0 
  for (j in 1:length(mf)) if (is.numeric(mf0[[j]])) {
    if (is.matrix(mf0[[j]])) { ## find row containing minimum
      j.min <- min((1:n)[as.logical(rowSums(mf[[j]]==min(mf[[j]])))])
      j.max <- min((1:n)[as.logical(rowSums(mf[[j]]==max(mf[[j]])))])
    } else { ## vector
      j.min <- min(which(mf[[j]]==min(mf[[j]])))
      j.max <- min(which(mf[[j]]==max(mf[[j]])))
    }
    k <- k + 1; mf0[k,] <- mf[j.min,]
    k <- k + 1; mf0[k,] <- mf[j.max,] 
  } else if (is.factor(mf[[j]])) { ## factor variable...
    ## randomly sample one row from each factor level...
    find <- apply(X=as.matrix(levels(mf[[j]])),MARGIN=1,FUN=fun,fac=mf[[j]],ind=ind)
    find <- find[is.finite(find)] ## in case drop.unused.levels==FALSE, so that there ar levels without rows
    nf <- length(find)
    mf0[(k+1):(k+nf),] <- mf[find,]
    k <- k + nf
  }

  RNGkind(kind[1], kind[2])
  assign(".Random.seed", seed, envir = .GlobalEnv)

  mf0
} ## mini.mf


bgam.fitd <- function (G, mf, gp ,scale , coef=NULL,etastart = NULL,
    mustart = NULL, offset = rep(0, nobs),rho=0, control = gam.control(), intercept = TRUE, 
    gc.level=0,nobs.extra=0,npt=1) {
## This is a version of bgam.fit1 designed for use with discretized covariates. 
## Difference to bgam.fit1 is that XWX, XWy and Xbeta are computed in C
## code using compressed versions of X. Parallelization of XWX formation
## is performed at the C level using openMP.
## Alternative fitting iteration using Choleski only, including for REML.
## Basic idea is to take only one Newton step for parameters per iteration
## and to control the step length to ensure that at the end of the step we
## are not going uphill w.r.t. the REML criterion...
    
    y <- mf[[gp$response]]
    weights <- G$w 
    conv <- FALSE
    nobs <- nrow(mf)
    offset <- G$offset 

    if (rho!=0) { ## AR1 error model
      
      ld <- 1/sqrt(1-rho^2) ## leading diagonal of root inverse correlation
      sd <- -rho*ld         ## sub diagonal
      N <- nobs    
      ## see rwMatrix() for how following are used...
      ar.row <- c(1,rep(1:N,rep(2,N))[-c(1,2*N)]) ## index of rows to reweight
      ar.weight <- c(1,rep(c(sd,ld),N-1))     ## row weights
      ar.stop <- c(1,1:(N-1)*2+1)    ## (stop[i-1]+1):stop[i] are the rows to reweight to get ith row
      if (!is.null(mf$"(AR.start)")) { ## need to correct the start of new AR sections...
        ii <- which(mf$"(AR.start)"==TRUE)
        if (length(ii)>0) {
          if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
          ar.weight[ii*2-2] <- 0 ## zero sub diagonal
          ar.weight[ii*2-1] <- 1 ## set leading diagonal to 1
        }
      }
    } else {## AR setup complete
      ar.row <- ar.weight <- ar.stop <- -1 ## signal no re-weighting
    }

    family <- G$family
    additive <- if (family$family=="gaussian"&&family$link=="identity") TRUE else FALSE
    variance <- family$variance
    dev.resids <- family$dev.resids
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object")
    valideta <- family$valideta
    if (is.null(valideta))
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu))
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
 
    eta <- if (!is.null(etastart))
         etastart
    else family$linkfun(mustart)
    
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
       stop("cannot find valid starting values: please specify some")
    dev <- sum(dev.resids(y, mu, weights))*2 ## just to avoid converging at iter 1

    conv <- FALSE
   
    G$coefficients <- rep(0,ncol(G$X))
    class(G) <- "gam"  
    
    ## need to reset response and weights to post initialization values
    ## in particular to deal with binomial properly...
    G$y <- y
    G$w <- weights

    Sl <- Sl.setup(G) ## setup block diagonal penalty object
    rank <- 0
    for (b in 1:length(Sl)) rank <- rank + Sl[[b]]$rank
    Mp <- ncol(G$X) - rank ## null space dimension
    Nstep <- 0
    for (iter in 1L:control$maxit) { ## main fitting loop 
      devold <- dev
      dev <- 0
      ## accumulate the QR decomposition of the weighted model matrix
      if (iter==1||!additive) {
        qrx <- list() 
        if (iter>1) {
          ## form eta = X%*%beta
          eta <- Xbd(G$Xd,coef,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop)
        }
        mu <- linkinv(eta)
        mu.eta.val <- mu.eta(eta)
        good <- mu.eta.val != 0
        mu.eta.val[!good] <- .1 ## irrelvant as weight is zero
        z <- (eta - offset) + (G$y - mu)/mu.eta.val
        w <- (G$w * mu.eta.val^2)/variance(mu)
        dev <- sum(dev.resids(G$y,mu,G$w))
      
        qrx$y.norm2 <- if (rho==0) sum(w*z^2) else   ## AR mod needed
          sum(rwMatrix(ar.stop,ar.row,ar.weight,sqrt(w)*z,trans=FALSE)^2) 
       
        ## form X'WX efficiently...
        qrx$R <- XWXd(G$Xd,w,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,npt,G$drop,ar.stop,ar.row,ar.weight)
        ## form X'Wz efficiently...
        qrx$f <- XWyd(G$Xd,w,z,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop,ar.stop,ar.row,ar.weight)
        if(gc.level>1) gc()
     
        ## following reparameterizes X'X and f=X'y, according to initial reparameterizarion...
        qrx$XX <- Sl.initial.repara(Sl,qrx$R,inverse=FALSE,both.sides=TRUE,cov=FALSE,nt=npt)
        qrx$Xy <- Sl.initial.repara(Sl,qrx$f,inverse=FALSE,both.sides=FALSE,cov=FALSE,nt=npt)  
        
        G$n <- nobs
      } else {  ## end of if (iter==1||!additive)
        dev <- qrx$y.norm2 - sum(coef*qrx$f) ## actually penalized deviance
      }
  
      if (control$trace)
         message(gettextf("Deviance = %s Iterations - %d", dev, iter, domain = "R-mgcv"))

      if (!is.finite(dev)) stop("Non-finite deviance")

      ## preparation for working model fit is ready, but need to test for convergence first
      if (iter>2 && abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
          conv <- TRUE
          #coef <- start
          break
      }

      ## use fast REML code
      ## block diagonal penalty object, Sl, set up before loop

      if (iter==1) { ## need to get initial smoothing parameters 
        lambda.0 <- initial.sp(qrx$R,G$S,G$off,XX=TRUE) ## note that this uses the untrasformed X'X in qrx$R
        ## convert intial s.p.s to account for L 
        lsp0 <- log(lambda.0) ## initial s.p.
        if (!is.null(G$L)) lsp0 <- 
          if (ncol(G$L)>0) as.numeric(coef(lm(lsp0 ~ G$L-1+offset(G$lsp0)))) else rep(0,0)
        n.sp <- length(lsp0) 
      }
     
      ## carry forward scale estimate if possible...
      if (scale>0) log.phi <- log(scale) else {
        if (iter==1) {
            if (is.null(coef)||qrx$y.norm2==0) lsp0[n.sp+1] <- log(var(as.numeric(G$y))*.05) else
               lsp0[n.sp+1] <- log(qrx$y.norm2/(nobs+nobs.extra))
        }
      }

      ## get beta, grad and proposed Newton step... 
      repeat { ## Take a Newton step to update log sp and phi
        lsp <- lsp0 + Nstep
        if (scale<=0) log.phi <- lsp[n.sp+1] 
        prop <- Sl.fitChol(Sl,qrx$XX,qrx$Xy,rho=lsp[1:n.sp],yy=qrx$y.norm2,L=G$L,rho0=G$lsp0,log.phi=log.phi,
                 phi.fixed=scale>0,nobs=nobs,Mp=Mp,nt=npt,tol=dev*.Machine$double.eps^.7)
        if (max(Nstep)==0) { 
          Nstep <- prop$step;lsp0 <- lsp;
          break 
        } else {
          if (sum(prop$grad*Nstep)>dev*1e-7) Nstep <- Nstep/2 else {
            Nstep <- prop$step;lsp0 <- lsp;break;
          }
        }
      } ## end of sp update

      coef <- Sl.initial.repara(Sl,prop$beta,inverse=TRUE,both.sides=FALSE,cov=FALSE)

      if (any(!is.finite(coef))) {
          conv <- FALSE
          warning(gettextf("non-finite coefficients at iteration %d",
                  iter))
          break
      }
    } ## end fitting iteration

    if (!conv)
       warning("algorithm did not converge")
   
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
         if (any(mu > 1 - eps) || any(mu < eps))
                warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (family$family == "poisson") {
            if (any(mu < eps))
                warning("fitted rates numerically 0 occurred")
    }
  Mp <- G$nsdf
  if (length(G$smooth)>1) for (i in 1:length(G$smooth)) Mp <- Mp + G$smooth[[i]]$null.space.dim
  scale <- exp(log.phi)
  reml <- (dev/scale - prop$ldetS + prop$ldetXXS + (length(y)-Mp)*log(2*pi*scale))/2
  if (rho!=0) { ## correct REML score for AR1 transform
    df <- if (is.null(mf$"(AR.start)")) 1 else sum(mf$"(AR.start)")
    reml <- reml - (nobs-df)*log(ld)
  }

  for (i in 1:ncol(prop$db)) prop$db[,i] <- ## d beta / d rho matrix
        Sl.initial.repara(Sl,as.numeric(prop$db[,i]),inverse=TRUE,both.sides=TRUE,cov=TRUE,nt=npt) 

  object <- list(db.drho=prop$db,
                 gcv.ubre=reml,mgcv.conv=conv,rank=prop$r,
                 scale.estimated = scale<=0,outer.info=NULL,
                 optimizer=c("perf","chol"))
  object$coefficients <- coef
  ## form linear predictor efficiently...
  object$linear.predictors <- Xbd(G$Xd,coef,G$kd,G$ks,G$ts,G$dt,G$v,G$qc,G$drop) + G$offset
  PP <- Sl.initial.repara(Sl,prop$PP,inverse=TRUE,both.sides=TRUE,cov=TRUE,nt=npt)
  F <- pmmult(PP,qrx$R,FALSE,FALSE,nt=npt)  ##crossprod(PP,qrx$R) - qrx$R contains X'WX in this case
  object$edf <- diag(F)
  object$edf1 <- 2*object$edf - rowSums(t(F)*F) 
  object$sp <- exp(lsp[1:n.sp]) 
  object$sig2 <- object$scale <- scale
  object$Vp <- PP * scale
  object$Ve <- pmmult(F,object$Vp,FALSE,FALSE,nt=npt) ## F%*%object$Vp
  ## sp uncertainty correction... 
  if (!is.null(G$L)) prop$db <- prop$db%*%G$L
  M <- ncol(prop$db) 
  if (M>0) {
    ev <- eigen(prop$hess,symmetric=TRUE)
    ind <- ev$values <= 0
    ev$values[ind] <- 0;ev$values[!ind] <- 1/sqrt(ev$values[!ind])
    rV <- (ev$values*t(ev$vectors))[,1:M]
    Vc <- pcrossprod(rV%*%t(prop$db),nt=npt)
  } else Vc <- 0
  Vc <- object$Vp + Vc  ## Bayesian cov matrix with sp uncertainty
  object$edf2 <- rowSums(Vc*qrx$R)/scale
  object$Vc <- Vc
  object$outer.info <- list(grad = prop$grad,hess=prop$hess)  
  object$AR1.rho <- rho
  object$R <- pchol(qrx$R,npt)
  piv <- attr(object$R,"pivot") 
  object$R[,piv] <- object$R   
  object$iter <- iter 
  object$wt <- w
  object$y <- G$y
  object$prior.weights <- G$w
  rm(G);if (gc.level>0) gc()
  object
} ## end bgam.fitd



bgam.fit <- function (G, mf, chunk.size, gp ,scale ,gamma,method, coef=NULL,etastart = NULL,
    mustart = NULL, offset = rep(0, nobs), control = gam.control(), intercept = TRUE, 
    cl = NULL,gc.level=0,use.chol=FALSE,nobs.extra=0,samfrac=1,npt=1)
{   y <- mf[[gp$response]]
    weights <- G$w
    conv <- FALSE
    nobs <- nrow(mf)
    ##nvars <- ncol(G$X)
    offset <- G$offset
    family <- G$family
    G$family <- gaussian() ## needed if REML/ML used
    variance <- family$variance
    dev.resids <- family$dev.resids
    ## aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object")
    valideta <- family$valideta
    if (is.null(valideta))
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu))
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
 
    ##coefold <- NULL
    eta <- if (!is.null(etastart))
         etastart
    else family$linkfun(mustart)
    
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
       stop("cannot find valid starting values: please specify some")
    dev <- sum(dev.resids(y, mu, weights))*2 ## just to avoid converging at iter 1
    ##boundary <- 
    conv <- FALSE
   
    G$coefficients <- rep(0,ncol(G$X))
    class(G) <- "gam"  
    
    ## need to reset response and weights to post initialization values
    ## in particular to deal with binomial properly...
    G$y <- y
    G$w <- weights

    ## set up cluster for parallel computation...

    if (!is.null(cl)&&inherits(cl,"cluster")) {
      n.threads <- length(cl)
      while(nobs/n.threads < ncol(G$X)) n.threads <- n.threads - 1
      if (n.threads < length(cl)) { 
        warning("Too many cluster nodes to use all efficiently")
      }
    } else n.threads <- 1

    if (n.threads>1) { ## set up thread argument lists
      ## number of obs per thread
      nt <- rep(ceiling(nobs/n.threads),n.threads)
      nt[n.threads] <- nobs - sum(nt[-n.threads])
      arg <- list()
      n1 <- 0
      for (i in 1:n.threads) {
        n0 <- n1+1;n1 <- n1+nt[i]
        ind <- n0:n1 ## this threads data block from mf
        n.block <- nt[i]%/%chunk.size ## number of full sized blocks
        stub <- nt[i]%%chunk.size ## size of end block
        if (n.block>0) {
          start <- (0:(n.block-1))*chunk.size+1
          stop <- (1:n.block)*chunk.size
          if (stub>0) {
            start[n.block+1] <- stop[n.block]+1
            stop[n.block+1] <- nt[i]
            n.block <- n.block+1
          } 
        } else {
          n.block <- 1
          start <- 1
          stop <- nt[i]
        }
        arg[[i]] <- list(nobs= nt[i],start=start,stop=stop,n.block=n.block,
                         linkinv=linkinv,dev.resids=dev.resids,gc.level=gc.level,
                         mu.eta=mu.eta,variance=variance,mf = mf[ind,],
                         eta = eta[ind],offset = offset[ind],G = G,use.chol=use.chol)
        arg[[i]]$G$w <- G$w[ind];arg[[i]]$G$model <- NULL
        arg[[i]]$G$y <- G$y[ind]
      }
    } else { ## single thread, requires single indices
      ## construct indices for splitting up model matrix construction... 
      n.block <- nobs%/%chunk.size ## number of full sized blocks
      stub <- nobs%%chunk.size ## size of end block
      if (n.block>0) {
        start <- (0:(n.block-1))*chunk.size+1
        stop <- (1:n.block)*chunk.size
        if (stub>0) {
          start[n.block+1] <- stop[n.block]+1
          stop[n.block+1] <- nobs
          n.block <- n.block+1
        } 
      } else {
        n.block <- 1
        start <- 1
        stop <- nobs
      }
   } ## single thread indices complete
 
    conv <- FALSE

    if (method=="fREML") Sl <- Sl.setup(G) ## setup block diagonal penalty object
   
    for (iter in 1L:control$maxit) { ## main fitting loop
       ## accumulate the QR decomposition of the weighted model matrix
       wt <- rep(0,0) 
       devold <- dev
       dev <- 0
       if (n.threads == 1) { ## use original serial update code     
         for (b in 1:n.block) {
           ind <- start[b]:stop[b]
           X <- predict(G,newdata=mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,block.size=length(ind))
           rownames(X) <- NULL
           if (is.null(coef)) eta1 <- eta[ind] else eta1 <- drop(X%*%coef) + offset[ind]
           mu <- linkinv(eta1) 
           y <- G$y[ind] ## G$model[[gp$response]] ## - G$offset[ind]
           weights <- G$w[ind]
           mu.eta.val <- mu.eta(eta1)
           good <- (weights > 0) & (mu.eta.val != 0)
           z <- (eta1 - offset[ind])[good] + (y - mu)[good]/mu.eta.val[good]
           w <- (weights[good] * mu.eta.val[good]^2)/variance(mu)[good]
           dev <- dev + sum(dev.resids(y,mu,weights))
           wt <- c(wt,w)
           w <- sqrt(w)
           ## note that QR may be parallel using npt>1, even under serial accumulation...
           if (b == 1) qrx <- qr.update(w*X[good,],w*z,use.chol=use.chol,nt=npt) 
           else qrx <- qr.update(w*X[good,],w*z,qrx$R,qrx$f,qrx$y.norm2,use.chol=use.chol,nt=npt)
           rm(X);if(gc.level>1) gc() ## X can be large: remove and reclaim
        }
        if (use.chol) { ## post proc to get R and f...
          y.norm2 <- qrx$y.norm2 
          qrx <- chol2qr(qrx$R,qrx$f,nt=npt)
          qrx$y.norm2 <- y.norm2
        }
      } else { ## use parallel accumulation 
         for (i in 1:length(arg)) arg[[i]]$coef <- coef
         res <- parallel::parLapply(cl,arg,qr.up) 
         ## single thread debugging version 
         #res <- list()
         #for (i in 1:length(arg)) {
         #  res[[i]] <- qr.up(arg[[i]])
         #}
        ## now consolidate the results from the parallel threads...
        if (use.chol) {
          R <- res[[1]]$R;f <- res[[1]]$f;dev <- res[[1]]$dev
          wt <- res[[1]]$wt;y.norm2 <- res[[1]]$y.norm2
          for (i in 2:n.threads) {
            R <- R + res[[i]]$R; f <- f + res[[i]]$f
            wt <- c(wt,res[[i]]$wt); dev <- dev + res[[i]]$dev
            y.norm2 <- y.norm2 + res[[i]]$y.norm2
          }         
          qrx <- chol2qr(R,f,nt=npt)
          qrx$y.norm2 <- y.norm2
        } else { ## proper QR
          R <- res[[1]]$R;f <- res[[1]]$f;dev <- res[[1]]$dev
          wt <- res[[1]]$wt;y.norm2 <- res[[1]]$y.norm2
          for (i in 2:n.threads) {
            R <- rbind(R,res[[i]]$R); f <- c(f,res[[i]]$f)
            wt <- c(wt,res[[i]]$wt); dev <- dev + res[[i]]$dev
            y.norm2 <- y.norm2 + res[[i]]$y.norm2
          }         
          ## use parallel QR here if npt>1...
          qrx <- if (npt>1) pqr2(R,npt) else qr(R,tol=0,LAPACK=TRUE) 
          f <- qr.qty(qrx,f)[1:ncol(R)]
          rp <- qrx$pivot;rp[rp] <- 1:ncol(R) # reverse pivot
          qrx <- list(R=qr.R(qrx)[,rp],f=f,y.norm2=y.norm2)
        }
      } 

      ## if the routine has been called with only a random sample of the data, then 
      ## R, f and ||y||^2 can be corrected to estimate the full versions...
 
      qrx$R <- qrx$R/sqrt(samfrac)
      qrx$f <- qrx$f/sqrt(samfrac)
      qrx$y.norm2 <- qrx$y.norm2/samfrac

      G$n <- nobs
      
      rss.extra <- qrx$y.norm2 - sum(qrx$f^2)
      
      if (control$trace)
         message(gettextf("Deviance = %s Iterations - %d", dev, iter, domain = "R-mgcv"))

      if (!is.finite(dev)) stop("Non-finite deviance")

      ## preparation for working model fit is ready, but need to test for convergence first
      if (iter>2 && abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
          conv <- TRUE
          coef <- start
          break
      }

      if (method=="GCV.Cp") {
         fit <- magic(qrx$f,qrx$R,G$sp,G$S,G$off,L=G$L,lsp0=G$lsp0,rank=G$rank,
                      H=G$H,C=matrix(0,0,ncol(qrx$R)),     ##C=G$C,
                      gamma=gamma,scale=scale,gcv=(scale<=0),
                      extra.rss=rss.extra,n.score=nobs+nobs.extra)
 
         post <- magic.post.proc(qrx$R,fit,qrx$f*0+1) 
      } else if (method=="fREML") { ## use fast REML code
        ## block diagonal penalty object, Sl, set up before loop
        um <- Sl.Xprep(Sl,qrx$R,nt=npt)
        lambda.0 <- initial.sp(qrx$R,G$S,G$off)
        lsp0 <- log(lambda.0) ## initial s.p.
        ## carry forward scale estimate if possible...
        if (scale>0) log.phi <- log(scale) else {
          if (iter>1) log.phi <- log(object$scale) else {
            if (is.null(coef)||qrx$y.norm2==0) log.phi <- log(var(as.numeric(G$y))*.05) else
               log.phi <- log(qrx$y.norm2/(nobs+nobs.extra))
          }
        }
        fit <- fast.REML.fit(um$Sl,um$X,qrx$f,rho=lsp0,L=G$L,rho.0=G$lsp0,
                             log.phi=log.phi,phi.fixed=scale>0,rss.extra=rss.extra,
                             nobs =nobs+nobs.extra,Mp=um$Mp,nt=npt)
        res <- Sl.postproc(Sl,fit,um$undrop,qrx$R,cov=FALSE,L=G$L,nt=npt)
        object <- list(coefficients=res$beta,db.drho=fit$d1b,
                       gcv.ubre=fit$reml,mgcv.conv=list(iter=fit$iter,
                       message=fit$conv),rank=ncol(um$X),
                       Ve=NULL,scale.estimated = scale<=0,outer.info=fit$outer.info,
                        optimizer=c("perf","newton"))
 
        if (scale<=0) { ## get sp's and scale estimate
          nsp <- length(fit$rho)
          object$sig2 <- object$scale <- exp(fit$rho[nsp])
          object$sp <- exp(fit$rho[-nsp]) 
          nsp <- length(fit$rho.full)
          object$full.sp <- exp(fit$rho.full[-nsp])
        } else { ## get sp's
          object$sig2 <- object$scale <- scale  
          object$sp <- exp(fit$rho)
          object$full.sp <- exp(fit$rho.full)
        }
        class(object)<-c("gam")               
      } else { ## method is one of "ML", "P-REML" etc...
        y <- G$y; w <- G$w; n <- G$n;offset <- G$offset
        G$y <- qrx$f
        G$w <- G$y*0+1
        G$X <- qrx$R
        G$n <- length(G$y)
        G$offset <- G$y*0
        G$dev.extra <- rss.extra
        G$pearson.extra <- rss.extra
        G$n.true <- nobs+nobs.extra
        object <- gam(G=G,method=method,gamma=gamma,scale=scale,control=gam.control(nthreads=npt))
        y -> G$y; w -> G$w; n -> G$n;offset -> G$offset
      }
     
      if (method=="GCV.Cp") { 
        object <- list()
        object$coefficients <- fit$b
        object$edf <- post$edf 
        object$edf1 <- post$edf1
        ##object$F <- post$F
        object$full.sp <- fit$sp.full
        object$gcv.ubre <- fit$score
        object$hat <- post$hat
        object$mgcv.conv <- fit$gcv.info 
        object$optimizer="magic"
        object$rank <- fit$gcv.info$rank
        object$Ve <- post$Ve
        object$Vp <- post$Vb
        object$sig2 <- object$scale <- fit$scale
        object$sp <- fit$sp
        names(object$sp) <- names(G$sp)
        class(object)<-c("gam")
      }

      coef <- object$coefficients
        
      if (any(!is.finite(coef))) {
          conv <- FALSE
          warning(gettextf("non-finite coefficients at iteration %d",
                  iter))
          break
      }
    } ## end fitting iteration

    if (method=="fREML") { ## do expensive cov matrix cal only at end
      res <- Sl.postproc(Sl,fit,um$undrop,qrx$R,cov=TRUE,scale=scale,L=G$L,nt=npt)
      object$edf <- res$edf
      object$edf1 <- res$edf1
      object$edf2 <- res$edf2
      ##object$F <- res$F
      object$hat <- res$hat
      object$Vp <- res$Vp
      object$Ve <- res$Ve
      object$Vc <- res$Vc
    }

    if (!conv)
       warning("algorithm did not converge")
   
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
         if (any(mu > 1 - eps) || any(mu < eps))
                warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (family$family == "poisson") {
            if (any(mu < eps))
                warning("fitted rates numerically 0 occurred")
    }
  object$R <- qrx$R    
  object$iter <- iter 
  object$wt <- wt
  object$y <- G$y
  object$prior.weights <- G$w
  rm(G);if (gc.level>0) gc()
  object
} ## end bgam.fit



bgam.fit2 <- function (G, mf, chunk.size, gp ,scale ,gamma,method, etastart = NULL,
    mustart = NULL, offset = rep(0, nobs), control = gam.control(), intercept = TRUE,npt=1)
## version using sparse full model matrix in place of QR update...
## not multi-threaded, due to anyway disappointing performance
{   G$y <- y <- mf[[gp$response]]
    weights <- G$w
    conv <- FALSE
    nobs <- nrow(mf)
    ##nvars <- ncol(G$X)
    offset <- G$offset
    family <- G$family
    G$family <- gaussian() ## needed if REML/ML used
    variance <- family$variance
    dev.resids <- family$dev.resids
    ##aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object")
    valideta <- family$valideta
    if (is.null(valideta))
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu))
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
 
    ##coefold <- NULL
    eta <- if (!is.null(etastart))
         etastart
    else family$linkfun(mustart)
    
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
       stop("cannot find valid starting values: please specify some")
    dev <- sum(dev.resids(y, mu, weights))*2 ## just to avoid converging at iter 1
    ##boundary <- 
    conv <- FALSE

    G$n <- nobs
    X <- G$X 
    ## need to reset response and weights to post initialization values
    ## in particular to deal with binomial properly...
    G$y <- y
    G$w <- weights

    conv <- FALSE
    for (iter in 1L:control$maxit) { ## main fitting loop
      devold <- dev
      if (iter>1) eta <- as.numeric(X%*%coef) + offset
      mu <- linkinv(eta)
      mu.eta.val <- mu.eta(eta)
      good <- (G$w > 0) & (mu.eta.val != 0)
      z <- (eta - offset)[good] + (y - mu)/mu.eta.val
      w <- (G$w[good] * mu.eta.val[good]^2)/variance(mu)[good]
      dev <- sum(dev.resids(y,mu,G$w))
      W <- Diagonal(length(w),sqrt(w))
      if (sum(good)<nobs) {
        XWX <- as(Matrix::crossprod(W%*%X[good,]),"matrix")
      } else {
        XWX <- as(Matrix::crossprod(W%*%X),"matrix")
      }      
      qrx <- list(R = chol(XWX))
      Wz <- W%*%z

      ## in following note that Q = WXR^{-1} 
      if (sum(good)<nobs) {
        qrx$f <- forwardsolve(t(qrx$R),as.numeric(t(X[good,])%*%(W%*%Wz)))
      } else {
        qrx$f <- forwardsolve(t(qrx$R),as.numeric(t(X)%*%(W%*%Wz)))
      }
      qrx$y.norm2 <- sum(Wz^2)
   
      rss.extra <- qrx$y.norm2 - sum(qrx$f^2)
      
      if (control$trace)
         message(gettextf("Deviance = %s Iterations - %d\n", dev, iter, domain = "R-mgcv"))

      if (!is.finite(dev)) stop("Non-finite deviance")

      ## preparation for working model fit is ready, but need to test for convergence first
      if (iter>2 && abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
          conv <- TRUE
         # coef <- start
          break
      }

      if (method=="GCV.Cp") {
         fit <- magic(qrx$f,qrx$R,G$sp,G$S,G$off,L=G$L,lsp0=G$lsp0,rank=G$rank,
                      H=G$H,C= matrix(0,0,ncol(qrx$R)), ##C=G$C,
                      gamma=gamma,scale=scale,gcv=(scale<=0),
                      extra.rss=rss.extra,n.score=G$n)
 
         post <- magic.post.proc(qrx$R,fit,qrx$f*0+1) 
      } else { ## method is "REML" or "ML"
        y <- G$y; w <- G$w; n <- G$n;offset <- G$offset
        G$y <- qrx$f
        G$w <- G$y*0+1
        G$X <- qrx$R
        G$n <- length(G$y)
        G$offset <- G$y*0
        G$dev.extra <- rss.extra
        G$pearson.extra <- rss.extra
        G$n.true <- n
        object <- gam(G=G,method=method,gamma=gamma,scale=scale,control=gam.control(nthreads=npt))
        y -> G$y; w -> G$w; n -> G$n;offset -> G$offset
      }
      gc()

      if (method=="GCV.Cp") { 
        object <- list()
        object$coefficients <- fit$b
        object$edf <- post$edf
        object$edf1 <- post$edf1
        #object$F <- post$F
        object$full.sp <- fit$sp.full
        object$gcv.ubre <- fit$score
        object$hat <- post$hat
        object$mgcv.conv <- fit$gcv.info 
        object$optimizer="magic"
        object$rank <- fit$gcv.info$rank
        object$Ve <- post$Ve
        object$Vp <- post$Vb
        object$sig2 <- object$scale <- fit$scale
        object$sp <- fit$sp
        names(object$sp) <- names(G$sp)
        class(object)<-c("gam")
      }

      coef <- object$coefficients
        
      if (any(!is.finite(coef))) {
          conv <- FALSE
          warning("non-finite coefficients at iteration ",
                  iter)
          break
      }
    } ## fitting iteration

    if (!conv)
       warning("algorithm did not converge")
   
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
         if (any(mu > 1 - eps) || any(mu < eps))
                warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (family$family == "poisson") {
            if (any(mu < eps))
                warning("fitted rates numerically 0 occurred")
    }
      
  object$iter <- iter  
  object$wt <- w
  object$R <- qrx$R
  object$y <- G$y
  object$prior.weights <- G$w
  rm(G);gc()
  object
} ## end bgam.fit2

ar.qr.up <- function(arg) {
## function to perform QR updating with AR residuals, on one execution thread
  if (arg$rho!=0) { ## AR1 error model
     ld <- 1/sqrt(1 - arg$rho^2) ## leading diagonal of root inverse correlation
     sd <- -arg$rho * ld         ## sub diagonal
  } 
  yX.last <- NULL
  qrx <- list(R=NULL,f=array(0,0),y.norm2=0) ## initial empty qr object
  for (i in 1:arg$n.block) {
    ind <- arg$start[i]:arg$end[i] 
    if (arg$rho!=0) { ## have to find AR1 transform...
       N <- arg$end[i]-arg$start[i]+1
       ## note first row implied by this transform
       ## is always dropped, unless really at beginning of data.
       row <- c(1,rep(1:N,rep(2,N))[-c(1,2*N)])
       weight <- c(1,rep(c(sd,ld),N-1))
       stop <- c(1,1:(N-1)*2+1)
       if (!is.null(arg$mf$"(AR.start)")) { ## need to correct the start of new AR sections...
           ii <- which(arg$mf$"(AR.start)"[ind]==TRUE)
           if (length(ii)>0) {
             if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
             weight[ii*2-2] <- 0 ## zero sub diagonal
             weight[ii*2-1] <- 1 ## set leading diagonal to 1
           }
       }
     } 
     ## arg$G$model <- arg$mf[ind,]
     w <- sqrt(arg$G$w[ind])
     X <- w*predict(arg$G,newdata=arg$mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,block.size=length(ind))
     y <- w*(arg$mf[ind,arg$response] - arg$offset[ind]) ## w*(arg$G$model[[arg$response]] - arg$offset[ind])
     if (arg$rho!=0) {
       ## Apply transform...
       if (arg$last&&arg$end[i]==arg$nobs) yX.last <- 
           c(y[nrow(X)],X[nrow(X),]) ## store final row, in case of update
       if (arg$first&&i==1) {
          X <- rwMatrix(stop,row,weight,X)
          y <- rwMatrix(stop,row,weight,y)
       } else {
          X <- rwMatrix(stop,row,weight,X)[-1,]
          y <- rwMatrix(stop,row,weight,y)[-1]
       } 
     } ## dealt with AR1      
     qrx <- qr.update(X,y,qrx$R,qrx$f,qrx$y.norm2,use.chol=arg$use.chol)
     rm(X);if (arg$gc.level>1) {gc()} ## X can be large: remove and reclaim
  } ## all blocks dealt with
  qrx$yX.last <- yX.last
  if (arg$gc.level>1) {rm(arg,w,y,ind);gc()}
  qrx
} ## ar.qr.up

pabapr <- function(arg) {
## function for parallel calling of predict.gam
## QUERY: ... handling?
  predict.gam(arg$object,newdata=arg$newdata,type=arg$type,se.fit=arg$se.fit,terms=arg$terms,
                        block.size=arg$block.size,newdata.guaranteed=arg$newdata.guaranteed,
                        na.action=arg$na.action)
}

predict.bam <- function(object,newdata,type="link",se.fit=FALSE,terms=NULL,exclude=NULL,
                        block.size=50000,newdata.guaranteed=FALSE,na.action=na.pass,
                        cluster=NULL,discrete=TRUE,n.threads=1,...) {
## function for prediction from a bam object, possibly in parallel
  ## remove some un-needed stuff from object
  if (discrete && !is.null(object$dinfo)) {
    return(predict.bamd(object,newdata,type,se.fit,terms,exclude,
                        block.size,newdata.guaranteed,na.action,n.threads,...))
  }
  object$Sl <- object$qrx <- object$R <- object$F <- object$Ve <-
  object$Vc <- object$G <- object$residuals <- object$fitted.values <-
  object$linear.predictors <- NULL
  gc()

  if (!is.null(cluster)&&inherits(cluster,"cluster")) { 
     ## require(parallel)
     n.threads <- length(cluster)
  } else n.threads <- 1
  if (missing(newdata)) n <- nrow(object$model) else n <- nrow(newdata)
  if (n < 100*n.threads) n.threads <- 1 ## not worth the overheads
  if (n.threads==1) { ## single threaded call
    if (missing(newdata)) return(
      predict.gam(object,newdata=object$model,type=type,se.fit=se.fit,terms=terms,exclude=exclude,
                        block.size=block.size,newdata.guaranteed=newdata.guaranteed,
                        na.action=na.action,...)
    ) else return(
      predict.gam(object,newdata=newdata,type=type,se.fit=se.fit,terms=terms,exclude=exclude,
                        block.size=block.size,newdata.guaranteed=newdata.guaranteed,
                        na.action=na.action,...))
  } else { ## parallel call...
    nt <- rep(floor(n/n.threads),n.threads)
    nt[1] <- n - sum(nt[-1])
    arg <- list()
    n1 <- 0
    for (i in 1:n.threads) { 
      n0 <- n1+1;n1 <- n1+nt[i]
      ind <- n0:n1 ## this thread's data block from mf
      arg[[i]] <- list(object=object,type=type,se.fit=se.fit,terms=terms,exclude=exclude,
                        block.size=block.size,newdata.guaranteed=newdata.guaranteed,
                        na.action=na.action)
      arg[[i]]$object$model <- object$model[1:2,] ## save space
      if (missing(newdata)) {
        arg[[i]]$newdata <- object$model[ind,]
      } else {
        arg[[i]]$newdata <- newdata[ind,]
      }
    } ## finished setting up arguments
    ## newdata and object no longer needed - all info in thread lists...
    if (!missing(newdata)) rm(newdata)
    rm(object)
    gc()
    res <- parallel::parLapply(cluster,arg,pabapr) ## perform parallel prediction
    gc()
    ## and splice results back together...
    if (type=="lpmatrix") {
      X <- res[[1]]
      for (i in 2:length(res)) X <- rbind(X,res[[i]])
      return(X)
    } else if (se.fit==TRUE) {
      rt <- list(fit = res[[1]]$fit,se.fit = res[[1]]$se.fit)
      if (type=="terms") {
        for (i in 2:length(res)) { 
          rt$fit <- rbind(rt$fit,res[[i]]$fit)
          rt$se.fit <- rbind(rt$se.fit,res[[i]]$se.fit)
        }
      } else {
        for (i in 2:length(res)) { 
          rt$fit <- c(rt$fit,res[[i]]$fit)
          rt$se.fit <- c(rt$se.fit,res[[i]]$se.fit)
        }
      }
      return(rt)
    } else { ## no se's returned
      rt <- res[[1]]
       if (type=="terms") {
        for (i in 2:length(res)) rt <- rbind(rt,res[[i]])
      } else {
        for (i in 2:length(res)) rt <- c(rt,res[[i]])
      }
      return(rt)
    } 
  }
} ## end predict.bam 


bam.fit <- function(G,mf,chunk.size,gp,scale,gamma,method,rho=0,
                    cl=NULL,gc.level=0,use.chol=FALSE,npt=1) 
## function that does big additive model fit in strictly additive case
{  ## first perform the QR decomposition, blockwise....
   n <- nrow(mf)
   if (rho!=0) { ## AR1 error model
     ld <- 1/sqrt(1-rho^2) ## leading diagonal of root inverse correlation
     sd <- -rho*ld         ## sub diagonal
   }

   if (n>chunk.size) { ## then use QR accumulation approach
     if (!is.null(cl)&&inherits(cl,"cluster")) { 
       n.threads <- length(cl)
       while(n/n.threads < ncol(G$X)) n.threads <- n.threads - 1
       if (n.threads < length(cl)) { 
         warning("Too many cluster nodes to use all efficiently")
       }
     } else n.threads <- 1

     G$coefficients <- rep(0,ncol(G$X))
     class(G) <- "gam"

     if (n.threads>1) { ## set up thread argument lists
       ## number of obs per thread
       nt <- rep(ceiling(n/n.threads),n.threads)
       nt[n.threads] <- n - sum(nt[-n.threads])
       arg <- list()
       n1 <- 0
       for (i in 1:n.threads) { 
         n0 <- n1+1;n1 <- n1+nt[i]
         if (i>1&&rho!=0) { ## need to start from end of last block if rho!=0
           n0 <- n0-1;nt[i] <- nt[i]+1 
         }   
         ind <- n0:n1 ## this thread's data block from mf
         n.block <- nt[i]%/%chunk.size ## number of full sized blocks
         stub <- nt[i]%%chunk.size ## size of end block
         if (n.block>0) { 
           ## each block is of size 
           start <- (0:(n.block-1))*chunk.size+1
           end <- start + chunk.size - 1
           if (stub>0) {
             start[n.block+1] <- end[n.block]+1
             end[n.block+1] <- nt[i]
             n.block <- n.block+1
           } 
           if (rho!=0) { ## then blocks must overlap
             ns <- length(start)
             if (ns>1) start[2:ns] <- start[2:ns]-1 
           }
         } else {
           n.block <- 1
           start <- 1
           end <- nt[i]
         }
         arg[[i]] <- list(nobs= nt[i],start=start,end=end,n.block=n.block,
                         rho=rho,mf = mf[ind,],gc.level=gc.level,
                         offset = G$offset[ind],G = G,response=gp$response,
                         first=FALSE,last=FALSE,use.chol=use.chol)
         if (i==1) arg[[1]]$first <- TRUE
         if (i==n.threads) arg[[i]]$last <- TRUE 
         arg[[i]]$G$w <- G$w[ind];arg[[i]]$G$model <- NULL
       }
     } else { ## single thread, requires single indices 
       n.block <- n%/%chunk.size ## number of full sized blocks
       stub <- n%%chunk.size ## size of end block
       if (stub>0) n.block <- n.block + 1
       start <- 0:(n.block-1)*chunk.size    ## block starts
       end <- start + chunk.size;           ## block ends
       end[n.block] <- n
       if (rho==0) start <- start + 1  ## otherwise most blocks go to 1 before block start
       start[1] <- 1  
     } 
    
     if (n.threads==1) { ## use original single thread method...
       qrx <- list(R=NULL,f=array(0,0),y.norm2=0) ## initial empty qr object
       for (i in 1:n.block) {
         ind <- start[i]:end[i] 
         if (rho!=0) {
           N <- end[i]-start[i]+1

           row <- c(1,rep(1:N,rep(2,N))[-c(1,2*N)])
           weight <- c(1,rep(c(sd,ld),N-1))
           stop <- c(1,1:(N-1)*2+1) 
           if (!is.null(mf$"(AR.start)")) { ## need to correct the start of new AR sections...
             ii <- which(mf$"(AR.start)"[ind]==TRUE)
             if (length(ii)>0) {
               if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
               weight[ii*2-2] <- 0 ## zero sub diagonal
               weight[ii*2-1] <- 1 ## set leading diagonal to 1
             }
           }
         } 
         #G$model <- mf[ind,]
         w <- sqrt(G$w[ind])
         X <- w*predict(G,newdata=mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,block.size=length(ind))
         y <- w*(mf[ind,gp$response]-G$offset[ind])  ## w*(G$model[[gp$response]] - G$offset[ind])
         if (rho!=0) {
           ## Apply transform...
           if (end[i]==n) yX.last <- c(y[nrow(X)],X[nrow(X),]) ## store final row, in case of update
           if (i==1) {
             X <- rwMatrix(stop,row,weight,X)
             y <- rwMatrix(stop,row,weight,y)
           } else {
             X <- rwMatrix(stop,row,weight,X)[-1,]
             y <- rwMatrix(stop,row,weight,y)[-1]
           } 
         }      

         qrx <- qr.update(X,y,qrx$R,qrx$f,qrx$y.norm2,use.chol=use.chol,nt=npt)
         rm(X)
         if (gc.level>1) {gc()} ## X can be large: remove and reclaim
       } ## end of single thread block loop
       if (use.chol) { ## post proc to get R and f...
          y.norm2 <- qrx$y.norm2 
          qrx <- chol2qr(qrx$R,qrx$f,nt=npt)
          qrx$y.norm2 <- y.norm2
        }
     } else { ## use parallel accumulation
     
       res <- parallel::parLapply(cl,arg,ar.qr.up)
       ## Single thread de-bugging...
       # res <- list()
       # for (i in 1:length(arg)) {
       #   res[[i]] <- ar.qr.up(arg[[i]])
       # }

       ## now consolidate the results from the parallel threads...
       R <- res[[1]]$R;f <- res[[1]]$f; ## dev <- res[[1]]$dev
       y.norm2 <- res[[1]]$y.norm2
       for (i in 2:n.threads) {
         if (use.chol) {
           R <- R + res[[i]]$R; f <- f + res[[i]]$f
         } else {
           R <- rbind(R,res[[i]]$R); f <- c(f,res[[i]]$f)
         }
         y.norm2 <- y.norm2 + res[[i]]$y.norm2
       } 
       if (use.chol) {
         qrx <- chol2qr(R,f,nt=npt)
         qrx$y.norm2 <- y.norm2
       } else { ## proper QR         
         ## use parallel QR if npt>1...
         qrx <- if (npt>1) pqr2(R,npt) else qr(R,tol=0,LAPACK=TRUE) 
         f <- qr.qty(qrx,f)[1:ncol(R)]
         rp <- qrx$pivot;rp[rp] <- 1:ncol(R) # reverse pivot
         qrx <- list(R=qr.R(qrx)[,rp],f=f,y.norm2=y.norm2)
       }
       yX.last <- res[[n.threads]]$yX.last
     } 
     G$n <- n
     G$y <- mf[[gp$response]]
   
   } else { ## n <= chunk.size
     if (rho==0) qrx <- qr.update(sqrt(G$w)*G$X,sqrt(G$w)*(G$y-G$offset),use.chol=use.chol,nt=npt) else {
       row <- c(1,rep(1:n,rep(2,n))[-c(1,2*n)])
       weight <- c(1,rep(c(sd,ld),n-1))
       stop <- c(1,1:(n-1)*2+1)
       if (!is.null(mf$"(AR.start)")) { ## need to correct the start of new AR sections...
         ii <- which(mf$"(AR.start)"==TRUE)
         if (length(ii)>0) {
           if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
           weight[ii*2-2] <- 0 ## zero sub diagonal
           weight[ii*2-1] <- 1 ## set leading diagonal to 1
         }
       }
       yX.last <- c(G$y[n],G$X[n,])  ## store final row, in case of update
       X <- rwMatrix(stop,row,weight,sqrt(G$w)*G$X)
       y <- rwMatrix(stop,row,weight,sqrt(G$w)*G$y)
       qrx <- qr.update(X,y,use.chol=use.chol,nt=npt)
   
       rm(X); if (gc.level>1) gc() ## X can be large: remove and reclaim
     } 
     if (use.chol) { ## post proc to get R and f...
        y.norm2 <- qrx$y.norm2 
        qrx <- chol2qr(qrx$R,qrx$f,nt=npt)
        qrx$y.norm2 <- y.norm2
     }
   }

   rss.extra <- qrx$y.norm2 - sum(qrx$f^2)
 
   if (method=="GCV.Cp") {
     fit <- magic(qrx$f,qrx$R,G$sp,G$S,G$off,L=G$L,lsp0=G$lsp0,rank=G$rank,
                H=G$H,C=matrix(0,0,ncol(qrx$R)),  ##C=G$C,
                gamma=gamma,scale=scale,gcv=(scale<=0),
                extra.rss=rss.extra,n.score=n)
 
     post <- magic.post.proc(qrx$R,fit,qrx$f*0+1) 
   } else if (method=="fREML"){ ## use fast REML code
     Sl <- Sl.setup(G) ## setup block diagonal penalty object
     um <- Sl.Xprep(Sl,qrx$R,nt=npt)
     lambda.0 <- initial.sp(qrx$R,G$S,G$off)
     lsp0 <- log(lambda.0) ## initial s.p.
     if (scale<=0) log.phi <- log(var(as.numeric(G$y))*.05) else ## initial phi guess
                   log.phi <- log(scale)
     fit <- fast.REML.fit(um$Sl,um$X,qrx$f,rho=lsp0,L=G$L,rho.0=G$lsp0,
            log.phi=log.phi,phi.fixed=scale>0,rss.extra=rss.extra,
            nobs =n,Mp=um$Mp,nt=npt)
     res <- Sl.postproc(Sl,fit,um$undrop,qrx$R,cov=TRUE,scale=scale,L=G$L,nt=npt)
     object <- list(coefficients=res$beta,edf=res$edf,edf1=res$edf1,edf2=res$edf2,##F=res$F,
                    db.drho=fit$d1b,
                    gcv.ubre=fit$reml,hat=res$hat,mgcv.conv=list(iter=fit$iter,
                    message=fit$conv),rank=ncol(um$X),
                    Ve=res$Ve,Vp=res$Vp,Vc=res$Vc,
                    scale.estimated = scale<=0,outer.info=fit$outer.info,
                    optimizer=c("perf","newton"))
     if (scale<=0) { ## get sp's and scale estimate
       nsp <- length(fit$rho)
       object$sig2 <- object$scale <- exp(fit$rho[nsp])
       object$sp <- exp(fit$rho[-nsp])
       nsp <- length(fit$rho.full)
       object$full.sp <- exp(fit$rho.full[-nsp])
     } else { ## get sp's
       object$sig2 <- object$scale <- scale  
       object$sp <- exp(fit$rho)
       object$full.sp <- exp(fit$rho.full)
     }
     
     if (rho!=0) { ## correct RE/ML score for AR1 transform
       df <- if (is.null(mf$"(AR.start)")) 1 else sum(mf$"(AR.start)")
       object$gcv.ubre <- object$gcv.ubre - (n-df)*log(ld)
     }

     G$X <- qrx$R;G$dev.extra <- rss.extra
     G$pearson.extra <- rss.extra;G$n.true <- n
     object$Sl <- Sl ## to allow for efficient update
     class(object)<-c("gam")
   } else { ## method is "ML", "P-REML" or similar
     y <- G$y; w <- G$w; n <- G$n;offset <- G$offset
     G$y <- qrx$f
     G$w <- G$y*0+1
     G$X <- qrx$R
     G$n <- length(G$y)
     G$offset <- G$y*0
     G$dev.extra <- rss.extra
     G$pearson.extra <- rss.extra
     G$n.true <- n
     object <- gam(G=G,method=method,gamma=gamma,scale=scale,control=gam.control(nthreads=npt))
     y -> G$y; w -> G$w; n -> G$n;offset -> G$offset
     if (rho!=0) { ## correct RE/ML score for AR1 transform 
       df <- if (is.null(mf$"(AR.start)")) 1 else sum(mf$"(AR.start)")
       object$gcv.ubre <- object$gcv.ubre - (n-df)*log(ld)
     }
   }
   if (method=="GCV.Cp") { 
     object <- list()
     object$coefficients <- fit$b
     object$edf <- post$edf
     object$edf1 <- post$edf1
     ##object$F <- post$F
     object$full.sp <- fit$sp.full
     object$gcv.ubre <- fit$score
     object$hat <- post$hat
     object$mgcv.conv <- fit$gcv.info 
     object$optimizer="magic"
     object$rank <- fit$gcv.info$rank
     object$Ve <- post$Ve
     object$Vp <- post$Vb
     object$sig2 <- object$scale <- fit$scale
     object$sp <- fit$sp
     class(object)<-c("gam")
   } else {
    
   }
   G$smooth <- G$X <- NULL
   object$prior.weights <- G$w
   object$AR1.rho <- rho
   if (rho!=0) { ## need to store last model matrix row, to allow update
     object$yX.last <- yX.last
   }
  
   object$R <- qrx$R
   object$gamma <- gamma;object$G <- G;object$qrx <- qrx ## to allow updating of the model
   object$y <- mf[[gp$response]]
   object$iter <- 1
   object
} # end of bam.fit

predict.bamd <- function(object,newdata,type="link",se.fit=FALSE,terms=NULL,exclude=NULL,
                        block.size=50000,newdata.guaranteed=FALSE,na.action=na.pass,n.threads=1,...) {
## function for prediction from a bam object, by discrete methods
  ## remove some un-needed stuff from object
  object$Sl <- object$qrx <- object$R <- object$F <- object$Ve <-
  object$Vc <- object$G <- object$residuals <- object$fitted.values <-
  object$linear.predictors <- NULL
  gc()
  if (missing(newdata)) newdata <- object$model

  if (type=="iterms") {
    type <- "terms"
    warning("iterms reset to terms")
  }
  
  if (!is.null(exclude)) warning("exclude ignored by discrete prediction at present")

  newdata <- predict.gam(object,newdata=newdata,type="newdata",se.fit=se.fit,terms=terms,exclude=exclude,
            block.size=block.size,newdata.guaranteed=newdata.guaranteed,
            na.action=na.action,...) 

  ## Parametric terms have to be dealt with safely, but without forming all terms 
  ## or a full model matrix. Strategy here is to use predict.gam, having removed
  ## key smooth related components from model object, so that it appears to be
  ## a parametric model... 
  offset <- 0
  if (object$nsdf) { ## deal with parametric terms...
    ## save copies of smooth info...
    smooth <- object$smooth; coef <- object$coefficients; Vp <- object$Vp
    ## remove key smooth info from object 
    object$coefficients <-  object$coefficients[1:object$nsdf]
    object$Vp <- object$V[1:object$nsdf,1:object$nsdf]
    object$smooth <- NULL
    ## get prediction for parametric component. Always "lpmatrix", unless terms required.
    ptype <- if (type %in% c("terms","iterms")) type else "lpmatrix"
    pp <- predict.gam(object,newdata=newdata,type=ptype,se.fit=se.fit,terms=terms,exclude=exclude,
            block.size=block.size,newdata.guaranteed=TRUE,
            na.action=na.action,...)  
    ## restore smooths to 'object'
    object$coefficients <- coef
    object$Vp <- Vp
    object$smooth <- smooth
    if (ptype=="lpmatrix") {
      offset <- attr(pp,"model.offset")
      if (is.null(offset)) offset <- 0
    }
  } ## parametric component dealt with

  ## now discretize covariates... 
  dk <- discrete.mf(object$dinfo$gp,mf=newdata,pmf=NULL,full=FALSE)
    
  Xd <- list() ### list of discrete model matrices...
  if (object$nsdf>0) {
     Xd[[1]] <- if (type%in%c("term","iterms")) matrix(0,0,0) else pp 
     kd <- cbind(1:nrow(newdata),dk$k) ## add index for parametric part to index list
     kb <- k <- 2; 
     dk$k.start <- c(1,dk$k.start+1) ## and adjust k.start accordingly
     dk$nr <- c(NA,dk$nr) ## need array index to match elements of Xd
  } else {
    kb <- k <- 1;  
    kd <- dk$k
  }
  ## k[,ks[j,1]:ks[j,2]] gives index columns for term j, thereby allowing 
  ## summation over matrix covariates....
  ks <- cbind(dk$k.start[-length(dk$k.start)],dk$k.start[-1])

  ts <- object$dinfo$ts
  dt <- object$dinfo$dt   
  for (i in 1:length(object$smooth)) { ## work through the smooth list
    ## first deal with any by variable (as first marginal of tensor)...
    if (object$smooth[[i]]$by!="NA") {
      by.var <- dk$mf[[object$smooth[[i]]$by]][1:dk$nr[k]]
      if (is.factor(by.var)) { 
         ## create dummy by variable...
         by.var <- as.numeric(by.var==object$smooth[[i]]$by.level)  
      }
      Xd[[k]] <- matrix(by.var,dk$nr[k],1)
      k <- k + 1
    }
    ## ... by done
    if (inherits(object$smooth[[i]],"tensor.smooth")) { 
      nmar <- length(object$smooth[[i]]$margin) 
      if (!is.null(object$smooth[[i]]$rind)) {
         ## terms re-ordered for efficiency, so the same has to be done on indices...
         rind <- k:(k+dt[kb]-1) ## could use object$dinfo$dt[kb]   
         dk$nr[rind] <- dk$nr[k+object$smooth[[i]]$rind-1] 
         ks[rind,] <- ks[k+object$smooth[[i]]$rind-1,] # either this line or next not both
         ##kd[,rind] <- kd[,k+object$smooth[[i]]$rind-1]
      } 
      XP <- object$smooth[[i]]$XP         
      for (j in 1:nmar) {
        Xd[[k]] <- PredictMat(smooth[[i]]$margin[[j]],dk$mf,n=dk$nr[k])
        if (!is.null(XP)&&(j<=length(XP))&&!is.null(XP[[j]])) Xd[[k]] <- Xd[[k]]%*%XP[[j]]
        k <- k + 1 
      }
    } else { ## not a tensor smooth
      object$smooth[[i]]$by <- "NA" ## have to ensure by not applied here!
      Xd[[k]] <- PredictMat(object$smooth[[i]],dk$mf,n=dk$nr[k])
      k <- k + 1
    }
    kb <- kb + 1
  }
   
  ## end of discrete set up
  se <- se.fit
  if (type=="terms") {
    if (object$nsdf>0) {
      if (se) {
        fit <- cbind(pp$fit,matrix(0,nrow(kd),length(object$smooth)))
        se.fit <- cbind(pp$se.fit,matrix(0,nrow(kd),length(object$smooth))) 
      } else fit <- cbind(pp,matrix(0,nrow(kd),length(object$smooth)))
      k <- 2; ## starting Xd
      kk <- ncol(fit) - length(object$smooth) + 1 ## starting col of fit for smooth terms
    } else {
      if (se) {
        fit <- matrix(0,nrow(kd),length(object$smooth))
        se.fit <- matrix(0,nrow(kd),length(object$smooth)) 
      } else fit <- matrix(0,nrow(kd),length(object$smooth))
      k <- 1; ## starting Xd
      kk <- 1 ## starting col of fit for smooth terms
    }
    for (i in 1:length(object$smooth)) {
      ii <- ts[k]:(ts[k]+dt[k]-1) ## index components for this term
      ind <- object$smooth[[i]]$first.para:object$smooth[[i]]$last.para ## index coefs for this term
      if (!is.null(object$dinfo$drop)) { 
        drop <- object$dinfo$drop-object$smooth[[i]]$first.para+1
        drop <- drop[drop<=length(ii)]
      } else drop <- NULL
      fit[,kk] <- Xbd(Xd[ii],object$coefficients[ind],kd,ks[ii,],                           ##kd[,ii,drop=FALSE]
                      1,dt[k],object$dinfo$v[k],object$dinfo$qc[k],drop=drop)
      if (se) se.fit[,kk] <- diagXVXd(Xd[ii],object$Vp[ind,ind],kd,ks[ii,],                 #kd[,ii,drop=FALSE],
                       1,dt[k],object$dinfo$v[k],object$dinfo$qc[k],drop=drop,n.threads=n.threads)^.5
      k <-  k + 1; kk <- kk + 1
    } 
    fit.names <- c(if (se) colnames(pp$fit) else colnames(pp),unlist(lapply(object$smooth,function(x) x$label)))
    colnames(fit) <- fit.names
    if (se) { 
      colnames(se.fit) <- fit.names
      fit <- list(fit=fit,se.fit=se.fit)
    }
  } else if (type=="lpmatrix") {
    fit <- Xbd(Xd,diag(length(object$coefficients)),kd,ks,ts,dt,object$dinfo$v,object$dinfo$qc,drop=object$dinfo$drop)
  } else { ## link or response
    fit <- Xbd(Xd,object$coefficients,kd,ks,ts,dt,object$dinfo$v,object$dinfo$qc,drop=object$dinfo$drop) + offset
    if (type=="response") {
      linkinv <- object$family$linkinv
      dmu.deta <- object$family$mu.eta
    } else linkinv <- dmu.deta <- NULL
    if (se==TRUE) {
      se.fit <- diagXVXd(Xd,object$Vp,kd,ks,ts,dt,object$dinfo$v,object$dinfo$qc,drop=object$dinfo$drop,n.threads=n.threads)^.5
      if (type=="response") {
        se.fit <- se.fit * abs(dmu.deta(fit))
        fit <- linkinv(fit)
      }
      fit <- list(fit=fit,se.fit=se.fit)
    } else if (type=="response") fit <- linkinv(fit)
  }
  fit
} ## end predict.bamd 




sparse.model.matrix <- function(G,mf,chunk.size) {
## create a whole sparse model matrix
  nobs = nrow(mf)
  n.block <- nobs%/%chunk.size ## number of full sized blocks
  stub <- nobs%%chunk.size ## size of end block
  if (n.block>0) {
    start <- (0:(n.block-1))*chunk.size+1
      stop <- (1:n.block)*chunk.size
      if (stub>0) {
        start[n.block+1] <- stop[n.block]+1
        stop[n.block+1] <- nobs
        n.block <- n.block+1
      } 
  } else {
    n.block <- 1
    start <- 1
    stop <- nobs
  }
  G$coefficients <- rep(0,ncol(G$X))
  class(G) <- "gam"

  X <- Matrix(0,nobs,ncol(G$X))   
  for (b in 1:n.block) {    
    ind <- start[b]:stop[b]
    #G$model <- mf[ind,]
    X[ind,] <- as(predict(G,newdata=mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,blocksize=length(ind)),"dgCMatrix")
    gc()
  }
  X
} # sparse.model.matrix

tero <- function(sm) {
## te smooth spec re-order so that largest marginal is last.
  maxd <- 0
  ns <- length(sm$margin)
  for (i in 1:ns) if (sm$margin[[i]]$bs.dim>=maxd) {
    maxi <- i;maxd <- sm$margin[[i]]$bs.dim
  }
  if (maxi<ns) { ## re-ordering required
    ind <- 1:ns;ind[maxi] <- ns;ind[ns] <- maxi
    sm$margin <- sm$margin[ind]
    sm$fix <- sm$fix[ind]
    sm$term <- rep("",0)
    for (i in 1:ns) sm$term <- c(sm$term,sm$margin[[i]]$term)
    sm$label <- paste0(substr(sm$label,1,3),paste0(sm$term,collapse=","),")",collapse="")
  }
  sm
} ## tero

AR.resid <- function(rsd,rho=0,AR.start=NULL) {
## standardised residuals for AR1 model
  if (rho==0) return(rsd)
  ld <- 1/sqrt(1-rho^2) ## leading diagonal of root inverse correlation
  sd <- -rho*ld         ## sub diagonal
  N <- length(rsd)    
  ## see rwMatrix() for how following are used...
  ar.row <- c(1,rep(1:N,rep(2,N))[-c(1,2*N)]) ## index of rows to reweight
  ar.weight <- c(1,rep(c(sd,ld),N-1))     ## row weights
  ar.stop <- c(1,1:(N-1)*2+1)    ## (stop[i-1]+1):stop[i] are the rows to reweight to get ith row
  if (!is.null(AR.start)) { ## need to correct the start of new AR sections...
    ii <- which(AR.start==TRUE)
    if (length(ii)>0) {
          if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
          ar.weight[ii*2-2] <- 0 ## zero sub diagonal
          ar.weight[ii*2-1] <- 1 ## set leading diagonal to 1
    }
  }
  rwMatrix(ar.stop,ar.row,ar.weight,rsd)
} ## AR.resid

bam <- function(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action=na.omit,
                offset=NULL,method="fREML",control=list(),select=FALSE,scale=0,gamma=1,knots=NULL,sp=NULL,
                min.sp=NULL,paraPen=NULL,chunk.size=10000,rho=0,AR.start=NULL,discrete=FALSE,
                sparse=FALSE,cluster=NULL,nthreads=NA,gc.level=1,use.chol=FALSE,samfrac=1,
                drop.unused.levels=TRUE,G=NULL,fit=TRUE,...)

## Routine to fit an additive model to a large dataset. The model is stated in the formula, 
## which is then interpreted to figure out which bits relate to smooth terms and which to 
## parametric terms.
## This is a modification of `gam' designed to build the QR decompostion of the model matrix 
## up in chunks, to keep memory costs down.
## If cluster is a parallel package cluster uses parallel QR build on cluster. 
## 'n.threads' is number of threads to use for non-cluster computation (e.g. combining 
## results from cluster nodes). If 'NA' then is set to max(1,length(cluster)).
{ control <- do.call("gam.control",control)
  if (control$trace) t3 <- t2 <- t1 <- t0 <- proc.time()
  if (is.null(G)) { ## need to set up model!
    if (is.character(family))
            family <- eval(parse(text = family))
    if (is.function(family))
            family <- family()
    if (is.null(family$family))
            stop("family not recognized")
    ##family = gaussian() ## no choise here
    if (family$family=="gaussian"&&family$link=="identity") am <- TRUE else am <- FALSE
    if (scale==0) { if (family$family%in%c("poisson","binomial")) scale <- 1 else scale <- -1} 
    if (!method%in%c("fREML","GCV.Cp","REML",
                    "ML","P-REML","P-ML")) stop("un-supported smoothness selection method")
    if (is.logical(discrete)) {
      discretize <- discrete
      discrete <- NULL ## use default discretization, if any
    } else {
      discretize <- if (is.numeric(discrete)) TRUE else FALSE
    }
    if (discretize) { 
      if (method!="fREML") { 
        discretize <- FALSE
        warning("discretization only available with fREML")
      } else {
        if (!is.null(cluster)) warning("discrete method does not use parallel cluster - use nthreads instead")
      }
    }
    if (method%in%c("fREML")&&!is.null(min.sp)) {
      min.sp <- NULL
      warning("min.sp not supported with fast REML computation, and ignored.")
    }
    if (sparse&&method%in%c("fREML")) {
      method <- "REML"
      warning("sparse=TRUE not supported with fast REML, reset to REML.")
    }
    gp <- interpret.gam(formula) # interpret the formula 
    if (discretize) { 
      ## re-order the tensor terms for maximum efficiency, and 
      ## signal that "re"/"fs" terms should be constructed with marginals
      ## also for efficiency
      if (length(gp$smooth.spec)>0) for (i in 1:length(gp$smooth.spec)) { 
        if (inherits(gp$smooth.spec[[i]],"tensor.smooth.spec")) 
        gp$smooth.spec[[i]] <- tero(gp$smooth.spec[[i]])
        if (inherits(gp$smooth.spec[[i]],c("re.smooth.spec","fs.smooth.spec"))&&gp$smooth.spec[[i]]$dim>1) {
          #gp$smooth.spec[[i]]$xt <- "tensor"
          class(gp$smooth.spec[[i]]) <- c(class(gp$smooth.spec[[i]]),"tensor.smooth.spec")
                                        ##c("re.smooth.spec","tensor.smooth.spec")
          gp$smooth.spec[[i]]$margin <- list()
          ## only ok for 'fs' with univariate metric variable (caught in 'fs' construcor)...
          for (j in 1:gp$smooth.spec[[i]]$dim) gp$smooth.spec[[i]]$margin[[j]] <- list(term=gp$smooth.spec[[i]]$term[j])
        }
      }
    }
    cl <- match.call() # call needed in gam object for update to work
    mf <- match.call(expand.dots=FALSE)
    mf$formula <- gp$fake.formula 
    mf$method <-  mf$family<-mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp <- mf$gc.level <-
    mf$gamma <- mf$paraPen<- mf$chunk.size <- mf$rho <- mf$sparse <- mf$cluster <- mf$discrete <-
    mf$use.chol <- mf$samfrac <- mf$nthreads <- mf$G <- mf$fit <- mf$select <- mf$...<-NULL
    mf$drop.unused.levels <- drop.unused.levels
    mf[[1]]<-as.name("model.frame")
    pmf <- mf
 
    pmf$formula <- gp$pf
    pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part
    pterms <- attr(pmf,"terms") ## pmf only used for this and discretization, if selected.
   
    if (gc.level>0) gc()

    mf <- eval(mf, parent.frame()) # the model frame now contains all the data
   # if ("matrix"%in%unlist(lapply(mf,class))) {
   #   mfattr <- attributes(mf)
   #   mf <- lapply(mf,drop) # avoid single column matrices
   #   mfattr -> attributes(mf)
   # }
    if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
    terms <- attr(mf,"terms")
    if (gc.level>0) gc()  
    if (rho!=0&&!is.null(mf$"(AR.start)")) if (!is.logical(mf$"(AR.start)")) stop("AR.start must be logical")
    
    ## summarize the *raw* input variables
    ## note can't use get_all_vars here -- buggy with matrices
    vars <- all.vars(gp$fake.formula[-2]) ## drop response here
    inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))

    ## allow a bit of extra flexibility in what `data' is allowed to be (as model.frame actually does)
    if (!is.list(data)&&!is.data.frame(data)) data <- as.data.frame(data) 

    dl <- eval(inp, data, parent.frame())
    if (!control$keepData) { rm(data);gc()} ## save space
    names(dl) <- vars ## list of all variables needed
    var.summary <- variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data
    rm(dl); if (gc.level>0) gc() ## save space    

    ## need mini.mf for basis setup, then accumulate full X, y, w and offset
    if (discretize) {
      ## discretize the data, creating list mf0 with discrete values
      ## and indices giving the discretized value for each element of model frame.
      ## 'discrete' can be null, or contain a discretization size, or
      ## a discretization size per smooth term.   
      dk <- discrete.mf(gp,mf,pmf,m=discrete)
      mf0 <- dk$mf ## padded discretized model frame
      sparse.cons <- 0 ## default constraints required for tensor terms

    } else { 
      mf0 <- mini.mf(mf,chunk.size)
      if (sparse) sparse.cons <- 2 else sparse.cons <- -1
    }
    rm(pmf); ## no further use
    if (control$trace) t1 <- proc.time()
    reset <- TRUE
    while (reset) {
      G <- gam.setup(gp,pterms=pterms,
                 data=mf0,knots=knots,sp=sp,min.sp=min.sp,
                 H=NULL,absorb.cons=TRUE,sparse.cons=sparse.cons,select=select,
                 idLinksBases=TRUE,scale.penalty=control$scalePenalty,
                 paraPen=paraPen,apply.by=!discretize)
      if (!discretize&&ncol(G$X)>=chunk.size) { ## no point having chunk.size < p
        chunk.size <- 4*ncol(G$X)
        warning(gettextf("chunk.size < number of coefficients. Reset to %d",chunk.size))
        if (chunk.size>=nrow(mf)) { ## no sense splitting up computation
          mf0 <- mf ## just use full dataset
        } else reset <- FALSE
      } else reset <- FALSE
    }
    if (control$trace) t2 <- proc.time()
    if (discretize) {
    
      v <- G$Xd <- list()
      ## have to extract full parametric model matrix from pterms and mf
      G$Xd[[1]] <- model.matrix(G$pterms,mf) 
      G$kd <- cbind(1:nrow(mf),dk$k) ## add index for parametric part to index list
      dk$k.start <- c(1,dk$k.start+1) ## and adjust k.start accordingly
      ## k[,ks[j,1]:ks[j,2]] gives index columns for term j, thereby allowing 
      ## summation over matrix covariates....
      G$ks <- cbind(dk$k.start[-length(dk$k.start)],dk$k.start[-1])
      ## create data object suitable for discrete data methods, from marginal model 
      ## matrices in G$smooth and G$X (stripping out padding, of course)
      if (ncol(G$Xd[[1]])) {
        kb <- k <- 2; qc <- dt <- ts <- rep(0,length(G$smooth)+1)
        dt[1] <- ts[1] <- 1;
        dk$nr <- c(NA,dk$nr) ## need array index to match elements of Xd
      } else {
        kb <- k <- 1; qc <- dt <- ts <- rep(0,length(G$smooth))
      }
      drop <- rep(0,0) ## index of te related columns to drop
      for (i in 1:length(G$smooth)) {
        ts[kb] <- k
        ## first deal with any by variable (as first marginal of tensor)...
        if (G$smooth[[i]]$by!="NA") {
          dt[kb] <- 1
          by.var <- dk$mf[[G$smooth[[i]]$by]][1:dk$nr[k]]
          if (is.factor(by.var)) { 
            ## create dummy by variable...
            by.var <- as.numeric(by.var==G$smooth[[i]]$by.level)  
          }
          G$Xd[[k]] <- matrix(by.var,dk$nr[k],1)
          k <- k + 1
        } else dt[kb] <- 0
        ## ... by done
        if (inherits(G$smooth[[i]],"tensor.smooth")) { 
          nmar <- length(G$smooth[[i]]$margin) 
          dt[kb] <- dt[kb] + nmar
          if (inherits(G$smooth[[i]],"fs.interaction")&&which(G$smooth[[i]]$fterm==G$smooth[[i]]$term)!=1) {
            ## have to reverse the terms because tensor representation assumes factor is first
            G$smooth[[i]]$rind <- 2:1 ## (k+1):k
          }          
          if (!is.null(G$smooth[[i]]$rind)) {
            ## terms re-ordered for efficiency, so the same has to be done on indices...
            rind <- k:(k+dt[kb]-1)    
            dk$nr[rind] <- dk$nr[k+G$smooth[[i]]$rind-1]
            G$ks[rind,] <- G$ks[k+G$smooth[[i]]$rind-1,] # either this line or next not both
            #G$kd[,rind] <- G$kd[,k+G$smooth[[i]]$rind-1]
          }       
          for (j in 1:nmar) {
            G$Xd[[k]] <- G$smooth[[i]]$margin[[j]]$X[1:dk$nr[k],,drop=FALSE]
            k <- k + 1 
          }
          ## deal with any side constraints on tensor terms  
          di <- attr(G$smooth[[i]],"del.index")
          if (!is.null(di)&&length(di>0)) {
            di <- di + G$smooth[[i]]$first.para + length(drop)  - 1
            drop <- c(drop,di)
          }
          ## deal with tensor smooth constraint
          qrc <- attr(G$smooth[[i]],"qrc")
          ## compute v such that Q = I-vv' and Q[,-1] is constraint null space basis
          if (inherits(qrc,"qr")) {
            v[[kb]] <- qrc$qr/sqrt(qrc$qraux);v[[kb]][1] <- sqrt(qrc$qraux)
            qc[kb] <- 1 ## indicate a constraint
          } else { 
            v[[kb]] <- rep(0,0) ##
            if (!inherits(qrc,"character")||qrc!="no constraints") warning("unknown tensor constraint type")
          } 
        } else { ## not a tensor smooth
          v[[kb]] <- rep(0,0)
          dt[kb] <- dt[kb] + 1
          G$Xd[[k]] <- G$X[1:dk$nr[k],G$smooth[[i]]$first.para:G$smooth[[i]]$last.para]
          k <- k + 1
        }
        kb <- kb + 1
      }
      if (length(drop>0)) G$drop <- drop ## index of terms to drop as a result of side cons on tensor terms

      ## ... Xd is the list of discretized model matrices, or marginal model matrices
      ## kd contains indexing vectors, so the ith model matrix or margin is Xd[[i]][kd[i,],]
      ## ts[i] is the starting matrix in Xd for the ith model matrix, while dt[i] is the number 
      ## of elements of Xd that make it up (1 for a singleton, more for a tensor). 
      ## v is list of Householder vectors encoding constraints and qc the constraint indicator.
      G$v <- v;G$ts <- ts;G$dt <- dt;G$qc <- qc
    } ## if (discretize)

    if (control$trace) t3 <- proc.time()
    G$sparse <- sparse

    ## no advantage to "fREML" with no free smooths...
    if (((!is.null(G$L)&&ncol(G$L) < 1)||(length(G$sp)==0))&&method=="fREML") method <- "REML"

    G$var.summary <- var.summary 
    G$family <- family
    G$terms<-terms;
    G$pred.formula <- gp$pred.formula

    n <- nrow(mf)
  
    if (is.null(mf$"(weights)")) G$w<-rep(1,n)
    else G$w<-mf$"(weights)"    
  
    G$offset <- model.offset(mf)  
    if (is.null(G$offset)) G$offset <- rep(0,n)

    if (ncol(G$X)>nrow(mf)) stop("Model has more coefficients than data") 
  
    if (ncol(G$X) > chunk.size && !discretize) { ## no sense having chunk.size < p
      chunk.size <- 4*ncol(G$X)
      warning(gettextf("chunk.size < number of coefficients. Reset to %d",chunk.size))    
    }

    G$cl<-cl;
    G$am <- am
     
    G$min.edf<-G$nsdf #-dim(G$C)[1]
    if (G$m) for (i in 1:G$m) G$min.edf<-G$min.edf+G$smooth[[i]]$null.space.dim
    G$discretize <- discretize
    G$formula<-formula
    ## environment(G$formula)<-environment(formula)
    environment(G$pterms) <- environment(G$terms) <- environment(G$pred.formula) <- 
    environment(G$formula) <- .BaseNamespaceEnv

  } else { ## G supplied
    scale <- G$scale
    mf <- G$mf; G$mf <- NULL
    gp <- G$gp; G$gp <- NULL
    na.action <- G$na.action; G$na.action <- NULL
  } ## end of G setup 

  if (!fit) {
    G$scale <- scale
    G$mf <- mf;G$na.action <- na.action;G$gp <- gp
    return(G)
  }


  ## number of threads to use for non-cluster node computation
  if (!is.finite(nthreads)||nthreads<1) nthreads <- max(1,length(cluster))

  G$conv.tol<-control$mgcv.tol      # tolerence for mgcv
  G$max.half<-control$mgcv.half     # max step halving in bfgs optimization


  ## now build up proper model matrix, and deal with y, w, and offset...

  if (control$trace) cat("Setup complete. Calling fit\n")
  
  colnamesX <- colnames(G$X)  

  if (G$sparse) { ## Form a sparse model matrix...
    warning("sparse=TRUE is deprecated")
    if (sum(G$X==0)/prod(dim(G$X))<.5) warning("model matrix too dense for any possible benefit from sparse")
    if (nrow(mf)<=chunk.size) G$X <- as(G$X,"dgCMatrix") else 
      G$X <- sparse.model.matrix(G,mf,chunk.size)
    if (rho!=0) warning("AR1 parameter rho unused with sparse fitting")
    object <- bgam.fit2(G, mf, chunk.size, gp ,scale ,gamma,method=method,
                       control = control,npt=nthreads,...)
  } else if (G$am&&!G$discretize) {
    if (nrow(mf)>chunk.size) G$X <- matrix(0,0,ncol(G$X)); if (gc.level>1) gc() 
    object <- bam.fit(G,mf,chunk.size,gp,scale,gamma,method,rho=rho,cl=cluster,
                      gc.level=gc.level,use.chol=use.chol,npt=nthreads)
  } else if (G$discretize) {
    object <- bgam.fitd(G, mf, gp ,scale ,nobs.extra=0,rho=rho,
                       control = control,npt=nthreads,gc.level=gc.level,...)
                       
  } else {
    G$X  <- matrix(0,0,ncol(G$X)); if (gc.level>1) gc()
    if (rho!=0) warning("AR1 parameter rho unused with generalized model")
    coef <- NULL
    if (samfrac<1 && samfrac>0) { ## sub-sample first to get close to right answer...
      ind <- sample(1:nrow(mf),ceiling(nrow(mf)*samfrac))
      if (length(ind)<2*ncol(G$X)) warning("samfrac too small - ignored") else {
        Gw <- G$w;Goffset <- G$offset
        G$w <- G$w[ind];G$offset <- G$offset[ind]
        control1 <- control
        control1$epsilon <- 1e-2
        object <- bgam.fit(G, mf[ind,], chunk.size, gp ,scale ,gamma,method=method,nobs.extra=0,
                       control = control1,cl=cluster,npt=nthreads,gc.level=gc.level,
                       use.chol=use.chol,samfrac=1,...)
        G$w <- Gw;G$offset <- Goffset
        coef <- object$coefficients
      }
    }
    ## fit full dataset
    object <- bgam.fit(G, mf, chunk.size, gp ,scale ,gamma,method=method,coef=coef,
                       control = control,cl=cluster,npt=nthreads,gc.level=gc.level,
                       use.chol=use.chol,...)
  }

  if (gc.level>0) gc()

  if (control$trace) t4 <- proc.time()

  if (control$trace) cat("Fit complete. Finishing gam object.\n")

  if (scale < 0) { object$scale.estimated <- TRUE;object$scale <- object$scale.est} else {
    object$scale.estimated <- FALSE; object$scale <- scale
  }

  object$assign <- G$assign # applies only to pterms  
  object$boundary <- FALSE  # always FALSE for this case
  object$call<-G$cl # needed for update() to work 
  object$cmX <- G$cmX ## column means of model matrix --- useful for CIs
 
  object$contrasts <- G$contrasts
  object$control <- control
  object$converged <- TRUE ## no iteration
  object$data <- NA ## not saving it in this case
  object$df.null <- nrow(mf)
  object$df.residual <- object$df.null - sum(object$edf) 
 
  object$family <- family
  object$formula<-G$formula 
 
  if (method=="GCV.Cp") {
    if (scale<=0) object$method <- "GCV" else object$method <- "UBRE"
  } else {
    object$method <- method
  }
  object$min.edf<-G$min.edf
  object$model <- mf;rm(mf);if (gc.level>0) gc()
  object$na.action <- attr(object$model,"na.action") # how to deal with NA's
  object$nsdf <- G$nsdf
  if (G$nsdf>0) names(object$coefficients)[1:G$nsdf] <- colnamesX[1:G$nsdf]
  object$offset <- G$offset
  ##object$prior.weights <- G$w
  object$pterms <- G$pterms
  object$pred.formula <- G$pred.formula 
  object$smooth <- G$smooth

  object$terms <- G$terms
  object$var.summary <- G$var.summary 
  if (is.null(object$wt)) object$weights <- object$prior.weights else
  object$weights <- object$wt
  object$xlevels <- G$xlevels
  #object$y <- object$model[[gp$response]]
  object$NA.action <- na.action ## version to use in bam.update
  names(object$sp) <- names(G$sp)
  if (!is.null(object$full.sp)) names(object$full.sp) <- names(G$lsp0)

  names(object$coefficients) <- G$term.names
  names(object$edf) <- G$term.names

  ## note that predict.gam assumes that it must be ok not to split the 
  ## model frame, if no new data supplied, so need to supply explicitly
  class(object) <- c("bam","gam","glm","lm")
  if (!G$discretize) { object$linear.predictors <- 
          as.numeric(predict.bam(object,newdata=object$model,block.size=chunk.size,cluster=cluster))
  } else { ## store discretization specific information to help with discrete prediction
    object$dinfo <- list(gp=gp, v = G$v, ts = G$ts, dt = G$dt, qc = G$qc, drop = G$drop)
  } 
  rm(G);if (gc.level>0) gc()

  object$fitted.values <- family$linkinv(object$linear.predictors)
   
  object$residuals <- sqrt(family$dev.resids(object$y,object$fitted.values,object$prior.weights)) * 
                      sign(object$y-object$fitted.values)
  if (rho!=0) object$std.rsd <- AR.resid(object$residuals,rho,object$model$"(AR.start)")

  dev <- object$deviance <- sum(object$residuals^2)
  if (rho!=0&&family$family=="gaussian") dev <- sum(object$std.rsd^2)
  object$aic <- family$aic(object$y,1,object$fitted.values,object$prior.weights,dev) -
                2 * (length(object$y) - sum(sum(object$model[["(AR.start)"]])))*log(1/sqrt(1-rho^2)) + ## correction for AR
                2*sum(object$edf)
  if (!is.null(object$edf2)&&sum(object$edf2)>sum(object$edf1)) object$edf2 <- object$edf1
  object$null.deviance <- sum(family$dev.resids(object$y,weighted.mean(object$y,object$prior.weights),object$prior.weights))
  if (!is.null(object$full.sp)) {
    if (length(object$full.sp)==length(object$sp)&&
        all.equal(object$sp,object$full.sp)==TRUE) object$full.sp <- NULL
  }
  environment(object$formula) <- environment(object$pred.formula) <-
  environment(object$terms) <- environment(object$pterms) <- 
  environment(attr(object$model,"terms"))  <- .GlobalEnv  
  if (control$trace) { 
    t5 <- proc.time()
    t5 <- rbind(t1-t0,t2-t1,t3-t2,t4-t3,t5-t4)[,1:3]
    row.names(t5) <- c("initial","gam.setup","pre-fit","fit","finalise")
    print(t5)
  }
  names(object$gcv.ubre) <- method
  object
} ## end of bam


bam.update <- function(b,data,chunk.size=10000) {
## update the strictly additive model `b' in the light of new data in `data'
## Need to update modelframe (b$model) 
  if (is.null(b$qrx)) { 
    stop("Model can not be updated")
  }
  gp<-interpret.gam(b$formula) # interpret the formula 
  
  X <- predict(b,newdata=data,type="lpmatrix",na.action=b$NA.action) ## extra part of model matrix
  rownames(X) <- NULL
  cnames <- names(b$coefficients)

  AR.start <- NULL ## keep R checks happy

  ## now get the new data in model frame form...
  getw <- "(weights)"%in%names(b$model)
  getARs <- "(AR.start)"%in%names(b$model)
  if (getw&&getARs) {
    mf <- model.frame(gp$fake.formula,data,weights=weights,AR.start=AR.start,
                      xlev=b$xlev,na.action=b$NA.action)
    w <- mf[["(weights)"]]
  } else if (getw) { 
    mf <- model.frame(gp$fake.formula,data,weights=weights,xlev=b$xlev,na.action=b$NA.action)
    w <- mf[["(weights)"]]
  } else if (getARs) {
    mf <- model.frame(gp$fake.formula,data,AR.start=AR.start,xlev=b$xlev,na.action=b$NA.action)
    w <- rep(1,nrow(mf))
  } else {
    mf <- model.frame(gp$fake.formula,data,xlev=b$xlev,na.action=b$NA.action)
    w <- rep(1,nrow(mf))
  }
  
  b$model <- rbind(b$model,mf) ## complete model frame --- old + new

  ## get response and offset...

  off.col <- attr(attr(b$model,"terms"),"offset")
  if (is.null(off.col)) offset <- rep(0,nrow(mf)) else offset <-  mf[,off.col]
  y <-  mf[,attr(attr(b$model,"terms"),"response")] - offset
  
  ## update G
  b$G$y <- c(b$G$y,y)
  b$G$offset <- c(b$G$offset,offset)
  b$G$w <- c(b$G$w,w)
  b$G$n <- nrow(b$model)
  n <- b$G$n;
  ## update the qr decomposition...

  w <- sqrt(w)

  if (b$AR1.rho!=0) { ## original model had AR1 error structure...
    rho <- b$AR1.rho
    ld <- 1/sqrt(1-rho^2) ## leading diagonal of root inverse correlation
    sd <- -rho*ld         ## sub diagonal
    ## append the final row of weighted X and y from original fit, first
    wy <- c(b$yX.last[1],w*y)
    wX <- rbind(b$yX.last[-1],w*X)
    m <- nrow(wX)
    b$yX.last <- c(wy[m],wX[m,])

    row <- c(1,rep(1:m,rep(2,m))[-c(1,2*m)])
    weight <- c(1,rep(c(sd,ld),m-1))
    stop <- c(1,1:(m-1)*2+1)
    if (!is.null(mf$"(AR.start)")) { ## need to correct the start of new AR sections...
         ii <- which(mf$"(AR.start)"==TRUE)
         if (length(ii)>0) {
           if (ii[1]==1) ii <- ii[-1] ## first observation does not need any correction
           weight[ii*2-2] <- 0 ## zero sub diagonal
           weight[ii*2-1] <- 1 ## set leading diagonal to 1
         }
    }
   
    ## re-weight to independence....
    wX <- rwMatrix(stop,row,weight,wX)[-1,]
    wy <- rwMatrix(stop,row,weight,wy)[-1]    

    ## update
    b$qrx <- qr.update(wX,wy,b$qrx$R,b$qrx$f,b$qrx$y.norm2)
  } else {
    b$qrx <- qr.update(w*X,w*y,b$qrx$R,b$qrx$f,b$qrx$y.norm2)
  }

  ## now do the refit...
  rss.extra <- b$qrx$y.norm2 - sum(b$qrx$f^2)

  if (b$method=="GCV"||b$method=="UBRE") method <- "GCV.Cp" else method <- b$method

 
  if (method=="GCV.Cp") {
    if (b$method=="GCV") scale <- -1 else scale = b$sig2
   
    fit <- magic(b$qrx$f,b$qrx$R,b$sp,b$G$S,b$G$off,L=b$G$L,lsp0=b$G$lsp0,rank=b$G$rank,
               H=b$G$H,C= matrix(0,0,ncol(b$qrx$R)),##C=b$G$C,
               gamma=b$gamma,scale=scale,gcv=(scale<=0),
               extra.rss=rss.extra,n.score=n)
 
    post <- magic.post.proc(b$qrx$R,fit,b$qrx$f*0+1) 
    b$y <- b$G$y;b$offset <- b$G$offset; b$G$w -> b$weights -> b$prior.weights;
    
  } else if (method=="fREML") { ## fast REML

     um <- Sl.Xprep(b$Sl,b$qrx$R)
     lsp0 <- log(b$sp) ## initial s.p.
     log.phi <- log(b$sig2) ## initial or fixed scale
     fit <- fast.REML.fit(um$Sl,um$X,b$qrx$f,rho=lsp0,L=b$G$L,rho.0=b$G$lsp0,
            log.phi=log.phi,phi.fixed = !b$scale.estimated,rss.extra=rss.extra,
            nobs =n,Mp=um$Mp,nt=1)
     if (b$scale.estimated) scale <- -1 else scale=b$sig2
     res <- Sl.postproc(b$Sl,fit,um$undrop,b$qrx$R,cov=TRUE,scale=scale,L=b$g$L)


     object <- list(coefficients=res$beta,edf=res$edf,edf1=res$edf1,edf2=res$edf2,##F=res$F,
                    gcv.ubre=fit$reml,hat=res$hat,outer.info=list(iter=fit$iter,
                    message=fit$conv),optimizer="fast-REML",rank=ncol(um$X),
                    Ve=NULL,Vp=res$V,Vc=res$Vc,db.drho=fit$d1b,scale.estimated = scale<=0)
     if (scale<=0) { ## get sp's and scale estimate
       nsp <- length(fit$rho)
       object$sig2 <- object$scale <- exp(fit$rho[nsp])
       object$sp <- exp(fit$rho[-nsp]) 
       nsp <- length(fit$rho.full)
       object$full.sp <- exp(fit$rho.full[-nsp])
     } else { ## get sp's
       object$sig2 <- object$scale <- scale  
       object$sp <- exp(fit$rho)
       object$full.sp <- exp(fit$rho.full)
     }
     
     if (b$AR1.rho!=0) { ## correct RE/ML score for AR1 transform
       df <- if (getARs) sum(b$model$"(AR.start)") else 1
       object$gcv.ubre <- object$gcv.ubre - (n-df)*log(ld)
     }

     b$G$X <- b$qrx$R;b$G$dev.extra <- rss.extra
     b$G$pearson.extra <- rss.extra;b$G$n.true <- n
     b$y <- b$G$y;b$offset <- b$G$offset; b$G$w -> b$weights -> b$prior.weights;

  } else { ## method is "REML" or "ML"
    y <- b$G$y; w <- b$G$w;offset <- b$G$offset
    b$G$y <- b$qrx$f
    b$G$w <- b$G$y*0+1
    b$G$X <- b$qrx$R
    b$G$n <- length(b$G$y)
    b$G$offset <- b$G$y*0
    b$G$dev.extra <- rss.extra
    b$G$pearson.extra <- rss.extra
    b$G$n.true <- n
    if (b$scale.estimated) scale <- -1 else scale = b$sig2
    in.out <- list(sp=b$sp,scale=b$reml.scale)
    object <- gam(G=b$G,method=method,gamma=b$gamma,scale=scale,in.out=in.out) 
    if (b$AR1.rho!=0) { ## correct RE/ML score for AR1 transform
       df <- if (getARs) sum(b$model$"(AR.start)") else 1
       object$gcv.ubre <- object$gcv.ubre - (n-df)*log(ld)
    }
    offset -> b$G$offset -> b$offset
    w -> b$G$w -> b$weights -> b$prior.weights; n -> b$G$n
    y -> b$G$y -> b$y;
  }
 
  if (method=="GCV.Cp") { 

    b$coefficients <- fit$b
    b$edf <- post$edf
    b$edf1 <- post$edf1
    ##b$F <- post$F
    b$full.sp <- fit$sp.full
    b$gcv.ubre <- fit$score
    b$hat <- post$hat
    b$mgcv.conv <- fit$gcv.info 
    b$optimizer="magic"
    b$rank <- fit$gcv.info$rank
    b$Ve <- post$Ve
    b$Vp <- post$Vb
    b$sig2 <- b$scale <- fit$scale
    b$sp <- fit$sp

  } else { ## REML or ML
    b$coefficients <- object$coefficients
    b$edf <- object$edf
    b$edf1 <- object$edf1
    ##b$F <- object$F
    b$full.sp <- object$sp.full
    b$gcv.ubre <- object$gcv.ubre
    b$hat <- object$hat
    b$outer.info <- object$outer.info 
    b$rank <- object$rank
    b$Ve <- object$Ve
    b$Vp <- object$Vp
    b$sig2 <- b$scale <- object$sig2
    b$sp <- object$sp
    if (b$AR1.rho!=0) { ## correct RE/ML score for AR1 transform
      b$gcv.ubre <- b$gcv.ubre - (n-1)*log(ld)
    }
  }

  b$R <- b$qrx$R
  b$G$X <- NULL
  b$linear.predictors <- as.numeric(predict.gam(b,newdata=b$model,block.size=chunk.size))
  b$fitted.values <- b$linear.predictor ## strictly additive only!
  
  b$residuals <- sqrt(b$family$dev.resids(b$y,b$fitted.values,b$prior.weights)) * 
                      sign(b$y-b$fitted.values)
  b$deviance <- sum(b$residuals^2)
  b$aic <- b$family$aic(b$y,1,b$fitted.values,b$prior.weights,b$deviance) +
           2 * sum(b$edf) 
  if (b$AR1.rho!=0) { ## correct aic for AR1 transform
    df <- if (getARs) sum(b$model$"(AR.start)") else 1
    b$aic <- b$aic + 2*(n-df)*log(ld)
  }
  b$null.deviance <- sum(b$family$dev.resids(b$y,mean(b$y),b$prior.weights))
  names(b$coefficients) <- names(b$edf) <- cnames
  b
} ## end of bam.update


#### ISSUES:   
## ? negative binomial support --- docs say it's there...
## offset unused in bam/bgam.fit, also gp only needed for "response",
## so could efficiently be replaced
