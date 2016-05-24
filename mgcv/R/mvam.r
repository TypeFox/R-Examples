## (c) Simon N. Wood (2013-2015) mvn model extended family. 
## Released under GPL2 ...

lpi.expand <- function(X,trailing=TRUE) {
## takes a model matrix X, with "lpi" attribute, and produces 
## full redundant version in which each column block is the full
## model matrix for one linear predictor, which may involve 
## repeating columns between blocks.
## See mvn family (ll) for prototypic application 
  lpi <- attr(X,"lpi")
  if (!attr(lpi,"overlap")) return(X) ## nothing to do
  ip <- unlist(lpi)
  if (trailing&&max(ip)<ncol(X)) { ## include any unindexed trailing blocks
    ii <- (max(ip)+1):ncol(X)
    X <- cbind(X[,ip],X[,ii])
  } else X <- X[,ip]
  k <- 0
  for (i in 1:length(lpi)) {
    lpi[[i]] <- 1:length(lpi[[i]]) + k
    k <- k + length(lpi[[i]])
  } 
  attr(X,"lpi") <- lpi
  X
} ## lpi.expand

lpi.contract <- function(x,lpi,type="rc",trailing=TRUE) {
## takes a vector or matrix x, and applies an lpi contraction to it
## if x is a matrix then type can be "r", "c" or "rc" for row, col
## or row, column contraction.
## See mvn family (ll) for prototypic application 
  ip <- unlist(lpi)
  if (trailing) { ## copy across any un-indexed trailing blocks
    lip <- length(ip) ## last row/col indexed in x
    ## get the length, nt, of any unindexed trailing block
    nt <- 0
    if (is.matrix(x)) { 
      if (type=="r"&&nrow(x)>lip) nt <- nrow(x)-lip else
      if (ncol(x)>lip) nt <- ncol(x) - lip
    } else if (length(x)>lip) nt <- length(x) - lip
    if (nt>0) { ## there is a trailing block - index it in lpi
      lpi[[length(lpi)+1]] <- 1:nt + max(ip)
      ip <- unlist(lpi)  
    }
  }
  p <- max(ip) ## dimension of result
  if (is.matrix(x)) {
    if (type=="c"||type=="rc") { ## column contraction
      k <- 0
      z <- matrix(0,nrow(x),p)
      for (i in 1:length(lpi)) { 
        ii <- 1:length(lpi[[i]]) + k
        k <- k + length(ii) 
        z[,lpi[[i]]] <- z[,lpi[[i]]] + x[,ii]
      } 
      if (type=="rc") x <- z
    }
    if (type=="r"||type=="rc") { ## row contraction
      z <- matrix(0,p,ncol(x))
      k <- 0
      for (i in 1:length(lpi)) { 
        ii <- 1:length(lpi[[i]]) + k
        k <- k + length(ii) 
        z[lpi[[i]],] <- z[lpi[[i]],] + x[ii,]
      }
    } 
  } else { ## vector
    z <- rep(0,p);k <- 0
    for (i in 1:length(lpi)) {
      ii <- 1:length(lpi[[i]]) + k
      k <- k + length(ii) 
      z[lpi[[i]]] <- z[lpi[[i]]] + x[ii]
    }
  }
  z
} ## lpi.contract

mvn <- function(d=2) { 
## Extended family object for multivariate normal additive model.
  if (d<2) stop("mvn requires 2 or more dimensional data")
  stats <- list()
  for (i in 1:d) {
    stats[[i]] <- make.link("identity") 
  }
  
  ##env <- new.env(parent = .GlobalEnv)
  validmu <- function(mu) all(is.finite(mu))

  ## initialization has to add in the extra parameters of 
  ## the cov matrix...
  
    preinitialize <- expression({
    ## code to evaluate in estimate.gam...
    ## extends model matrix with dummy columns and 
    ## finds initial coefficients
      ydim <- ncol(G$y) ## dimension of response
      nbeta <- ncol(G$X)
      ntheta <- ydim*(ydim+1)/2 ## number of cov matrix factor params
      lpi <- attr(G$X,"lpi")
      XX <- crossprod(G$X)
      G$X <- cbind(G$X,matrix(0,nrow(G$X),ntheta)) ## add dummy columns to G$X
      #G$cmX <- c(G$cmX,rep(0,ntheta)) ## and corresponding column means
      G$term.names <- c(G$term.names,paste("R",1:ntheta,sep="."))
      attr(G$X,"lpi") <- lpi
      attr(G$X,"XX") <- XX
      ## pad out sqrt of balanced penalty matrix to account for extra params
      attr(G$Sl,"E") <- cbind(attr(G$Sl,"E"),matrix(0,nbeta,ntheta))
      G$family$data <- list(ydim = ydim,nbeta=nbeta)
      G$family$ibeta = rep(0,ncol(G$X))
      ## now get initial parameters and store in family...
      for (k in 1:ydim) {
        sin <- G$off %in% lpi[[k]]
        #Sk <- G$S[sin]
        um <- magic(G$y[,k],G$X[,lpi[[k]]],rep(-1,sum(sin)),G$S[sin],
                    match(G$off[sin],lpi[[k]]), ## turn G$off global indices into indices for this predictor
                    nt=control$nthreads)
        G$family$ibeta[lpi[[k]]] <- um$b
        G$family$ibeta[nbeta+1] <- -.5*log(um$scale) ## initial log root precision
        nbeta <- nbeta + ydim - k + 1
      }
    })
    
    postproc <- expression({
    ## code to evaluate in estimate.gam, to do with estimated factor of
    ## precision matrix, etc...
      ydim <- G$family$data$ydim
      R <- matrix(0,ydim,ydim)
      ind <- G$family$data$nbeta + 1:(ydim*(ydim+1)/2);
      theta <- object$coefficients[ind]
      k <- 1;for (i in 1:ydim) for (j in i:ydim) {
        if (i==j) R[i,j] <- exp(theta[k]) else R[i,j] <- theta[k]
        k <- k + 1
      }
      object$family$data <- list(R=R) 
      rsd <- R%*%t(object$y-object$fitted.values)
      object$deviance <- sum(rsd^2)
      rsd <- R%*%(t(object$y)-colMeans(object$y))
      object$null.deviance <- sum(rsd^2)
    })
    
    initialize <- expression({
      ## called in gam.fit5 and initial.spg
        n <- rep(1, nobs)
        if (is.null(start)) start <- family$ibeta
        ## need to re-parameterize XX is non-standard
        if (exists("rp",inherits=FALSE)&&length(rp$rp)>0) 
           attr(x,"XX") <- Sl.repara(rp$rp,t(Sl.repara(rp$rp,attr(x,"XX"))))
    })


    residuals <- function(object,type=c("response","deviance")) {
      type <- match.arg(type)
      res <- object$y - object$fitted.values
      if (type=="deviance") res <- res%*%t(object$family$data$R)
      res 
    } ## residuals


    ##rd <- qf <- NULL ## these functions currently undefined for 

    ll <- function(y,X,coef,wt,family,deriv=0,d1b=NULL,d2b=NULL,Hp=NULL,rank=0,fh=NULL,D=NULL) {
    ## function defining the Multivariate Normal model log lik.
    ## Calls C code "mvn_ll"
    ## deriv codes: 0   - evaluate the log likelihood
    ##              1   - evaluate the grad and Hessian, H, of log lik w.r.t. coefs (beta)
    ##              2/3 - evaluate d1H =dH/drho given db/drho in d1b 
    ##                    (2 is diagonal only - not implemented efficiently)
    ##              4 -  given d1b and d2b evaluate trHid2H= tr(Hp^{-1}d2H/drhodrho') (not implemented)
    ## Hp is the preconditioned penalized hessian of the log lik
    ##    which is of rank 'rank'.
    ## fh is a factorization of Hp - either its eigen decomposition 
    ##    or its Choleski factor
    ## D is the diagonal pre-conditioning matrix used to obtain Hp
    ##   if Hr is the raw Hp then Hp = D*t(D*Hr)
      lpi <- attr(X,"lpi") ## lpi[[k]] is index of model matrix columns for dim k 
      overlap <- attr(lpi,"overlap") ## do dimensions share terms?
      drop <- attr(X,"drop") 
      if (!is.null(drop)) { ## the optimizer has dropped some parameters
        ## it will have adjusted lpi automatically, but XX is mvn specific
        attr(X,"XX") <- attr(X,"XX")[-drop,-drop]
      }
      m <- length(lpi)  ## number of dimensions of MVN
      if (overlap) { ## linear predictors share terms - expand to redundant representation 
        ip <- unlist(lpi)
        XX <- attr(X,"XX")[ip,ip]
        X <- lpi.expand(X)
        attr(X,"XX") <- XX;rm(XX)
        lpi0 <- lpi ## need to save this for contraction of results
        lpi <- attr(X,"lpi") ## this indexes the cols of each l.p in the expanded X
        ## need to expand coef beta, leaving m*(m+1)/2 final coefs of R at end
        ind <- (max(ip)+1):length(coef)
        if (length(ind)!=m*(m+1)/2) stop("mvn dimension error")
        coef <- c(coef[ip],coef[ind])
        ## do same for derivatives of coef wrt log smoothing params...
        if (!is.null(d1b)) d1b <- rbind(d1b[ip,],d1b[ind,]) 
      } else ind <- NULL
      lpstart <- rep(0,m)
      for (i in 1:(m-1)) lpstart[i] <- lpi[[i+1]][1]
      lpstart[m] <- lpi[[m]][length(lpi[[m]])]+1 
      nb <- length(coef)  ## total number of parameters
      if (deriv<2) {
        nsp = 0;d1b <- dH <- 0
      } else {
        nsp = ncol(d1b)
        dH = rep(0,nsp*nb*nb)
      }
      #cat("\nderiv=",deriv,"  lpstart=",lpstart," dim(y) = ",dim(y),
      #    "\ndim(XX)=",dim(attr(X,"XX"))," m=",m," nsp=",nsp,"\n")
      oo <- .C("mvn_ll",y=as.double(t(y)),X=as.double(X),XX=as.double(attr(X,"XX")),
               beta=as.double(coef),n=as.integer(nrow(X)),
               lpi=as.integer(lpstart-1),m=as.integer(m),ll=as.double(0),lb=as.double(coef*0),
               lbb=as.double(rep(0,nb*nb)), dbeta = as.double(d1b), dH = as.double(dH), 
               deriv = as.integer(nsp>0),nsp = as.integer(nsp),nt=as.integer(1),PACKAGE="mgcv")
      lb <- oo$lb;lbb <- matrix(oo$lbb,nb,nb)
      if (overlap) { ## need to apply lpi contraction
        lb <- lpi.contract(lb,lpi0)
        ## lpi.contract will automatically carry across the R coef related stuff 
        lbb <- lpi.contract(lbb,lpi0)
      }
      if (nsp==0) d1H <- NULL else if (deriv==2) {
        d1H <- matrix(0,nb,nsp)
        for (i in 1:nsp) { 
          dH <- matrix(oo$dH[ind],nb,nb)
          if (overlap) dH <- lpi.contract(dH,lpi0)
          d1H[,i] <- diag(dH)
          ind <- ind + nb*nb
        }
      } else { ## deriv==3
        d1H <- list();ind <- 1:(nb*nb)
        for (i in 1:nsp) { 
          dH <- matrix(oo$dH[ind],nb,nb)
          if (overlap) dH <- lpi.contract(dH,lpi0)
          d1H[[i]] <- dH
          ind <- ind + nb*nb
        }
      }
      list(l=oo$ll,lb=lb,lbb=lbb,d1H=d1H)
    } ## ll

    # environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
    # environment(rd)<- environment(qf)<- environment(variance) <- environment(putTheta) 
    ##environment(aic) <- 
    ##environment(ll) <- env
    structure(list(family = "Multivariate normal", 
        ## link = linktemp, linkfun = stats$linkfun, linkinv = stats$linkinv, 
        ll=ll,nlp=d,
        initialize = initialize,preinitialize=preinitialize,postproc=postproc,
        residuals=residuals,
        validmu = validmu, ## valideta = stats$valideta, 
        ## rd=rd,qf=qf,  
        linfo = stats, ## link information list
        d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
        ls=1, ## signal ls not needed
        available.derivs = 1 ## signal only first derivatives available...
        ),
        class = c("general.family","extended.family","family"))
} ## mvn

