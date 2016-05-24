#******************************************************************************* 
#
# Estimation for Multivariate Normal Data with Monotone Missingness
# Copyright (C) 2007, University of Cambridge
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## bmonomvn:
##
## Under a MVN model, sample from the posterior distribution of the
## mean and variance matrix from a data matrix y that is potentially
## subject to monotone missingness.  The rows need not be sorted to
## give a monotone appearance, but the pattern must be monotone.

'bmonomvn' <-
function(y, pre=TRUE, p=0.9, B=100, T=200, thin=1, economy=FALSE,
         method=c("lasso", "ridge", "lsr", "factor", "hs", "ng"),
         RJ=c("p", "bpsn", "none"), capm=TRUE, #capm=method!="lasso",
         start=NULL, mprior= 0, rd=NULL, theta=0, rao.s2=TRUE, QP=NULL,
         verb=1, trace=FALSE)
  {
    ## (quitely) double-check that blasso is clean before-hand
    bmonomvn.cleanup(1)
    
    ## what to do if fatally interrupted?
    on.exit(bmonomvn.cleanup(verb))
    
    ## save column names in a data frame, and then work with a matrix
    nam <- colnames(y)
    y <- as.matrix(y)
  
    ## dimensions of the inputs
    M <- ncol(y)
    N <- nrow(y)
    
    ## check pre
    if(length(pre) != 1 || !is.logical(pre))
      stop("pre must be a scalar logical")
    
    ## check B
    if(length(B) != 1 || B < 0)
      stop("B must be a scalar integer >= 0")
    
    ## check T
    if(length(T) != 1 || (T <= 0 && B > 0))
      stop("if B>0 when T must be a scalar integer >= 1")

    ## check thin
    if(length(thin) != 1 || thin <= 0)
      stop("thin must be a scalar > 0")

    ## check rao.s2
    if(length(rao.s2) != 1 || !is.logical(rao.s2))
      stop("rao.s2 must be a scalar logical")
    
    ## check economy
    if(length(economy) != 1 || !is.logical(economy))
      stop("economy must be a scalar logical")

    ## check trace
    if(length(trace) != 1 || !is.logical(trace))
      stop("trace must be a scalar logical")

    ## check method, and change to an integer
    method <- match.arg(method)
    mi <- 0  ## lasso
    if(method == "ridge") mi <- 1
    else if(method == "lsr") mi <- 2
    else if(method == "factor") mi <- 3
    else if(method == "hs") mi <- 4
    else if(method == "ng") mi <- 5

    ## check p argument
    if(method == "factor" && (length(p) != 1 || p < 1)) {
      warning("should have scalar p >= 1 for factor method, using default p=1")
      p <- 1
    } else if(length(p) != 1 || p > 1 || p < 0) {
      warning("should have scalar 0 <= p <= 1, using default p=1")
      p <- 1
    }
    
    ## check RJ, and change to an integer
    RJ <- match.arg(RJ)
    RJi <- 0  ## bpsn: "big-p small-n"
    if(RJ == "p") RJi <- 1 ## only do RJ when parsimonious is activated
    else if(RJ == "none") RJi <- 2 ## no RJ

    ## disallow bad RJ combination
    if(method == "lsr" && RJ == "none")
      stop("bad method (", method, ") and RJ (", RJ, ") combination", sep="")
    
    ## check capm
    if(length(capm) != 1 || !is.logical(capm))
      stop("capm must be a scalar logical or \"p\"\n")
    if(capm == TRUE && RJ == "none")
      stop("capm = TRUE is invalid when RJ = \"none\"")

    ## check mprior
    if(any(mprior < 0)) stop("must have all(0 <= mprior)");
    if(length(mprior) == 1) {
      if(mprior != 0 && RJ == FALSE)
        warning(paste("setting mprior=", mprior,
                      " ignored since RJ=FALSE", sep=""))
      if(mprior > 1) stop("must have scalar 0 <= mprior < 1")
      mprior <- c(mprior, 0)
    } else if(length(mprior) != 2)
      stop("mprior should be a scalar or 2-vector in [0,1]")

    ## check r and delta (rd), or default
    if(is.null(rd)) {
      if(method == "lasso" || method == "ng") rd <- c(2,0.1)
      else if(method == "ridge") rd <- c(5, 10)
      else rd <- c(0,0)
    }
    if(length(rd) != 2 || (method=="lasso" && any(rd <= 0)))
      stop("rd must be a positive 2-vector")
    if(method == "ng" && rd[1] != 2)
      stop("must have rd[1] = 2 for NG prior")

    ## check theta
    if(length(theta) != 1)# || theta < 0)
      stop("theta must be a non-negative scalar")
    
    ## save the call
    cl <- match.call()

    ## reorder the rows and columsn of y
    dap <- da.perm(y, pre)
    nas <- dap$nas; nao <- dap$nao;
    y <- dap$y; R <- dap$R
    n <- N - apply(R, 2, function(x){ sum(x == 1) }) ## number of non-ones
    ## if(sum(R == 2) != 0) stop("missingness pattern in y is not monotone")

    ## create the reverse ordering
    if(!is.null(nao)) oo <- order(nao)
    else oo <- NULL

    ## get the indices of the factors in the case of method = "factor"
    if(method == "factor") {
      if(is.null(oo)) facts <- 1:p
      else facts <- oo[1:p]
    } else facts <- double(0)
    
    ## check the start argument
    start <- check.start(start, nao, M)

    ## check the QP argument
    QPin <- check.QP(QP, M, nao, oo)
    if(is.logical(QP) && QP == FALSE) QP <- NULL

    ## save old y and then replace NA with zeros
    Y <- y
    Y[is.na(Y)] <- 0

    ## possibly add a new long if verbose printing
    if(verb >=1) cat("\n")
    
    ## call the C routine
    r <- .C("bmonomvn_R",

            ## begin estimation inputs
            B = as.integer(B),
            T = as.integer(T),
            thin = as.double(thin),
            M = as.integer(M),
            N = as.integer(N),
            Y = as.double(t(Y)),
            n = as.integer(n),
            R = as.integer(t(R)),
            p = as.double(p),
            mi = as.integer(mi),
            facts.len = as.integer(length(facts)),
            facts = as.integer(facts-1),
            RJi = as.integer(RJi),
            capm = as.integer(capm),
            smu.len = as.integer(length(start$mu)),
            smu = as.double(start$mu),
            sS.len = as.integer(length(start$S)),
            sS = as.double(start$S),
            sncomp.len = as.integer(length(start$ncomp)),
            sncomp = as.integer(start$ncomp),
            slambda.len = as.integer(length(start$lambda)),
            slambda = as.double(start$lambda),
            mprior = as.double(mprior),
            rd = as.double(rd),
            theta = as.double(theta),
            rao.s2 = as.integer(rao.s2),
            economy = as.integer(economy),
            verb = as.integer(verb),
            trace = as.integer(trace),

            ## begin Quadratic Progamming inputs
            QPcols = as.integer(c(length(QPin$cols), QPin$cols-1)),
            QPd = as.double(QPin$dvec),
            QPdmu = as.integer(QPin$dmu),
            QPA = as.double(QPin$Amat),
            QPb0 = as.double(QPin$b0),
            QPmc = as.integer(QPin$mu.constr),
            QPq = as.integer(QPin$q),
            QPmeq = as.integer(QPin$meq),
            
            ## begin estimation outputs
            mu = double(M),
            mu.var = double(M),
            mu.cov = double(M*M),
            S = double(M*M),
            S.var = double(M*M),
            mu.map = double(M),
            S.map = double(M*M),
            S.nz = double(M*M),
            Si.nz = double(M*M),
            lpost.map = double(1),
            which.map = integer(1),
            llik = double(T),
            llik.norm.len = as.integer(T * (theta != 0)),
            llik.norm = double(T * (theta != 0)),
            methods = integer(M),
            thin.act = integer(M),
            nu.len = as.integer(T*(theta < 0)),
            nu = double(T*(theta < 0)),
            lambda2 = double(M),
            ncomp = double(M),

            ## begin Quadratic Programming outputs
            W.len = as.integer(T*length(QPin$cols)*(!is.null(QPin$Amat))),
            W = double(T*length(QPin$cols)*(!is.null(QPin$Amat))),
            
            PACKAGE = "monomvn")

    ## copy the inputs back into the returned R-object
    r$Y <- NULL; r$y <- y;
    if(sum(R == 2) == 0) r$R <- NULL
    else r$R <- R

    ## remove lengths
    r$facts.len <- r$smu.len <- r$sS.len <- r$sncomp.len <- r$slambda.len <- NULL
    r$llik.norm.len <- r$nu.len <- r$W.len <- NULL 

    ## make S and other covars into matrices
    r$mu.cov <- matrix(r$mu.cov, ncol=M)
    r$S <- matrix(r$S, ncol=M)
    r$S.var <- matrix(r$S.var, ncol=M)
    r$S.map <- matrix(r$S.map, ncol=M)
    r$S.nz <- matrix(r$S.nz, ncol=M)
    r$Si.nz <- matrix(r$Si.nz, ncol=M)

    ## possibly add column permutation info from pre-processing
    if(pre) {
      r$na <- nas
      r$o <- nao
    }
    
    ## extract the methods
    mnames <- c("bcomplete", "brjlasso", "brjng", "brjhs", "brjridge", "brjlsr", 
                "blasso", "bng", "bhs", "bridge", "blsr")
    r$methods <- mnames[r$methods]
    
    ## put the original ordering back
    if(pre) {
      r$mu <- r$mu[oo]
      r$mu.var <- r$mu.var[oo]
      r$mu.map <- r$mu.map[oo]
      r$mu.cov <- r$mu.cov[oo,oo]
      r$S <- r$S[oo,oo]
      r$S.var <- r$S.var[oo,oo]
      r$S.map <- r$S.map[oo,oo]
      r$Si.nz <- r$Si.nz[oo,oo]
      r$ncomp <- r$ncomp[oo]
      r$lambda2 <- r$lambda2[oo]
      r$methods <- r$methods[oo]
      r$thin <- r$thin.act[oo]
    }

    ## deal with names
    if(! is.null(nam)) {
      r$mu <- matrix(r$mu, nrow=length(r$mu))
      rownames(r$mu) <- nam
      r$mu.var <- matrix(r$mu.var, nrow=length(r$mu.var))
      rownames(r$mu.var) <- nam
      r$mu.map <- matrix(r$mu.map, nrow=length(r$mu.map))
      rownames(r$mu.map) <- nam
      colnames(r$mu.cov) <- rownames(r$mu.cov) <- nam
      colnames(r$S) <- rownames(r$S) <- nam
      colnames(r$S.var) <- rownames(r$S.var) <- nam
      colnames(r$S.map) <- rownames(r$S.map) <- nam
      colnames(r$Si.nz) <- rownames(r$Si.nz) <- nam
      r$ncomp <- matrix(r$ncomp, nrow=length(r$ncomp))
      rownames(r$ncomp) <- nam
      r$lambda2 <- matrix(r$lambda2, nrow=length(r$lambda2))
      rownames(r$lambda2) <- nam
    }

    ## read the trace in the output files, and then delete them
    if(trace)
      r$trace <- bmonomvn.read.traces(r$N, r$n, r$M, nao, oo, nam,
                                      capm, mprior, R, cl,
                                      method == "hs", r$thin, r$verb)
    else r$trace <- NULL
    
    ## final line
    if(verb >= 1) cat("\n")

    ## null-out redundancies
    r$n <- r$N <- r$M <- r$mi <- r$verb <- NULL
    r$smu <- r$sS <- r$sncomp <- r$slambda <- NULL
    r$thin.act <- NULL
    if(r$theta == 0) { r$theta <- r$llik.norm <- NULL }
    else if(r$theta > 0) r$nu <- NULL

    ## change back to logicals or original inputs
    r$rao.s2 <- as.logical(r$rao.s2)
    r$capm <- as.logical(r$capm)
    r$economy <- as.logical(r$economy)
    r$RJi <- NULL; r$RJ <- RJ
    r$facts <- facts
    if(r$mprior[2] == 0) r$mprior <- r$mprior[-2]

    ## off-by-one
    r$which.map <- r$which.map + 1

    ## record Quadratic Programming info
    r$QPd <- r$QPA <- r$QPb0 <- r$QPq <- r$QPmeq <- NULL;
    r$QPmc <- r$QPcols <- r$QPdmu <- NULL
    if(!is.null(QP)) {
      QPW <- postprocess.QP(QPin, r$W, M, T, nam)
      r$QP <- QPW$QP; r$W <- QPW$W
    } else { r$QP <- r$W <- NULL }
    
    ## assign class, call and methods, and return
    r$call <- cl
    class(r) <- "monomvn"
    return(r)
  }


## check.start:
##
## sanity check the format of the start vector, and then
## re-arrange the components into the monotone order
## specified in nao, which should agree with start$o

check.start <- function(start, nao, M)
{
  s <- list(mu=NULL, S=NULL, ncomp=NULL)
  
  if(!is.null(start)) {
    
    ## make sure orders are the smae
    if(!is.null(start$o) && start$o != nao)
      stop("starting monotone order is not the same as y's order")
    
    ## check and the reorder mu
    if(!is.null(start$mu) && length(start$mu) == M) s$mu <- start$mu[nao]
    else stop("start$mu must be specified and have length ncol(y)")
    
    ## check and then reorder S
    if(!is.null(start$S) && nrow(start$S) == ncol(start$S) &&
       nrow(start$S) == M)
      s$S <- start$S[nao, nao]
    else stop("start$S must be specified, be square, and have dim = ncol(y)")
    
    ## check and then reorder ncomp
    if(!is.null(start$ncomp)) {
      s$ncomp <- start$ncomp
      na <- is.na(s$ncomp)
      s$ncomp <- s$ncomp[nao]
      s$ncomp[na[nao]] <- (0:(length(s$ncomp)-1))[na[nao]]
      if(length(s$ncomp) != M || !is.integer(s$ncomp) || any(s$ncomp < 0) )
        stop("start$ncomp must be non-neg integer vector of length ncol(y)")
    } else s$ncomp <- 0:(M-1)

    ## check and then reorder lambda
    if(!is.null(start$lambda)) {
      s$lambda <- start$lambda
      na <- is.na(s$lambda)
      s$lambda <- s$lambda[nao]
      s$lambda[na[nao]] <- 0
      if(length(s$lambda) != M || any(s$lambda < 0) )
        stop("start$lambda should be a non-neg vector of length ncol(y)")
    } else s$lambda <- rep(0, M)

  }
  return(s)
}


## bmonomvn.cleanup
##
## gets called when the C-side is aborted by the R-side and enables
## the R-side to clean up the memory still allocaed to the C-side,
## as well as whatever files were left open on the C-side

"bmonomvn.cleanup" <-  function(verb)
{
  .C("bmonomvn_cleanup", PACKAGE="monomvn")

  ## should a newline be appended after trace removals
  nl <- FALSE
  
  ## get rid of trace of the mean samples (mu)
  if(file.exists(paste("./", "mu.trace", sep=""))) {
    unlink("mu.trace")
    if(verb >= 1) cat("NOTICE: removed mu.trace\n")
    nl <- TRUE
  }

  ## get rid of trace of the Covar samples (S)
  if(file.exists(paste("./", "S.trace", sep=""))) {
    unlink("S.trace")
    if(verb >= 1) cat("NOTICE: removed S.trace\n")
    nl <- TRUE
  }

  ## get rid of trace of the Covar samples (S)
  if(file.exists(paste("./", "DA.trace", sep=""))) {
    unlink("DA.trace")
    if(verb >= 1) cat("NOTICE: removed DA.trace\n")
    nl <- TRUE
  }

  ## get all of the names of the tree files
  b.files <- list.files(pattern="blasso_M[0-9]+_n[0-9]+.trace")
    
  ## for each tree file
  if(length(b.files > 0)) {
    for(i in 1:length(b.files)) {
      if(verb >= 1) cat(paste("NOTICE: removed ", b.files[i], "\n", sep=""))
      unlink(b.files[i])
      nl <- TRUE
    }    
  }

  ## final newline
  if(verb >= 1 && nl) cat("\n")
}


## da.perm:
## 
## Re-order the columns and rows of y to follow the
## monotone missingness pattern with the least number
## of violations.  Return a re-ordered y, the missingness
## matrix R with entries: 0 (observed), 1 (monotone missing)
## and 2 (requires data augmentation), and vectors
## describing the row and column permutations performed

da.perm <- function(y, pre = TRUE)
  {
    ## dimensions of y
    n <- nrow(y)
    m <- ncol(y)
    
    ## forst re-order by column missingness
    nas <- apply(y, 2, function(x) {sum(is.na(x))} )
    
    ## check for cols with all NAs
    if(sum(nas == n) > 0) {
      cat("cols with no data:\n")
      print((1:m)[nas == n])
      stop("remove these columns and try again")
    }
    
    ## re-order the columns to follow the monotone pattern
    if(pre) {
      nao <- order(nas)
      y <- y[,nao]
    } else nao <- 1:m
    
    ## calculate initial R with ones and zeros only
    R <- matrix(0, ncol=ncol(y), nrow=nrow(y))
    R[is.na(y)] <- 1
    
    ## now see which indicies violate motonicity, and
    ## put a 2 in the places found
    for (j in 2:ncol(y)) {
      
      ## First check each of the next group of columns for monotonicity
      for(i in 1:nrow(y)) {
        if(R[i,j] == 0) { ## entry not missing

          ## get entries which violate monotonicity
          da <- (1:(j-1))[R[i,1:(j-1)] == 1]
          R[i,da] <- 2
        }
      } 
    }

    ## re-order R and y by the numebr of 1s in the rows of R
    r <- apply(R, 1, function(x){ sum(x == 1) })
    or <- order(r)
    y <- y[or,]
    R <- R[or,]
    
    ## return R and y, in addition to the col and row orderings
    return(list(y=y, R=R, nas=nas, nao=nao, or=or))
  }
