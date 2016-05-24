## (c) Simon N. Wood (2013, 2014) coxph model general family. 
## Released under GPL2 ...

cox.ph <- function (link = "identity") { 
## Extended family object for Cox PH.
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("identity")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for coxph family; available link is \"identity\" ")
  }
  env <- new.env(parent = .GlobalEnv)
  validmu <- function(mu) all(is.finite(mu))

    

    ## initialization is tough here... need data frame in reverse time order,
    ## and intercept removed from X...
  
    preinitialize <- expression({
    ## code to evaluate in estimate.gam...
      ## sort y (time) into decending order, and
      ## re-order weights and rows of X accordingly
      G$family$data <- list()
      y.order <- order(G$y,decreasing=TRUE)
      G$family$data$y.order <- y.order
      G$y <- G$y[y.order]
      G$X <- G$X[y.order,,drop=FALSE]
      G$w <- G$w[y.order]
    })
    
    postproc <- expression({
    ## code to evaluate in estimate.gam, to do with data ordering and 
    ## baseline hazard estimation...
      ## first get the estimated hazard and prediction information...
      object$family$data <- G$family$hazard(G$y,G$X,object$coefficients,G$w)
      rumblefish <- G$family$hazard(G$y,matrix(0,nrow(G$X),0),object$coefficients,G$w)
      s0.base <- exp(-rumblefish$h[rumblefish$r]) ## no model baseline survival 
      s0.base[s0.base >= 1] <- 1 - 2*.Machine$double.eps ## avoid NA later
      ## now put the survivor function in object$fitted
      object$fitted.values <- exp(-object$family$data$h[object$family$data$r]*exp(object$linear.predictors))
      ## compute the null deviance...
      s.base <- exp(-object$family$data$h[object$family$data$r]) ## baseline survival
      s.base[s.base >= 1] <- 1 - 2*.Machine$double.eps ## avoid NA later
      object$null.deviance <- ## sum of squares of null deviance residuals
      2*sum(abs((object$prior.weights + log(s0.base) + object$prior.weights*(log(-log(s0.base)))))) 
      ## and undo the re-ordering...
      object$linear.predictors[y.order] <- object$linear.predictors
      object$fitted.values[y.order] <- object$fitted.values
      object$y[y.order] <- object$y  
      object$prior.weights[y.order] <- object$prior.weights
    })
    
    initialize <- expression({
        n <- rep(1, nobs)
        if (is.null(start)) start <- rep(0,ncol(x))
    })

    hazard <- function(y, X,beta,wt) {
    ## get the baseline hazard function information, given times in descending order in y
    ## model matrix (same ordering) in X, coefs in beta and censoring in wt (1 = death, 0
    ## = censoring)
      tr <- unique(y);r <- match(y,tr);nt <- length(tr)
      oo <- .C("coxpp",as.double(X%*%beta),A=as.double(X),as.integer(r),d=as.integer(wt),
               h=as.double(rep(0,nt)),q=as.double(rep(0,nt)),km=as.double(rep(0,nt)),
               n=as.integer(nrow(X)),p=as.integer(ncol(X)),
               nt=as.integer(nt),PACKAGE="mgcv")
      p <- ncol(X)
      list(tr=tr,h=oo$h,q=oo$q,a=matrix(oo$A[p*nt],p,nt),nt=nt,r=r,km=oo$km)
    }

    residuals <- function(object,type=c("deviance","martingale")) {
      type <- match.arg(type)
      w <- object$prior.weights;log.s <- log(object$fitted.values)
      res <- w + log.s ## martingale residuals
      if (type=="deviance") { 
        log.s[log.s>-1e-50] <- -1e-50
        res <- sign(res)*sqrt(-2*(res + w * log(-log.s)))
      }
      res 
    }


    predict <- function(family,se=FALSE,eta=NULL,y=NULL,
               X=NULL,beta=NULL,off=NULL,Vb=NULL) {
      ## prediction function.
      if (sum(is.na(y))>0) stop("NA times supplied for cox.ph prediction")
      ii <- order(y,decreasing=TRUE) ## C code expects non-increasing
      n <- nrow(X)
      oo <- .C("coxpred",as.double(X[ii,]),t=as.double(y[ii]),as.double(beta),as.double(Vb),
                a=as.double(family$data$a),h=as.double(family$data$h),q=as.double(family$data$q),
                tr = as.double(family$data$tr),
                n=as.integer(n),p=as.integer(ncol(X)),nt = as.integer(family$data$nt),
                s=as.double(rep(0,n)),se=as.double(rep(0,n)),PACKAGE="mgcv")
      s <- sef <- oo$s
      s[ii] <- oo$s
      sef[ii] <- oo$se    
      if (se) return(list(fit=s,se.fit=sef)) else return(list(fit=s))
    }

    rd <- qf <- NULL ## these functions currently undefined for Cox PH

    ll <- function(y,X,coef,wt,family,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
    ## function defining the cox model log lik.
    ## Calls C code "coxlpl"
    ## deriv codes: 0   - evaluate the log likelihood
    ##              1   - evaluate the grad and Hessian, H, of log lik w.r.t. coefs (beta)
    ##              2/3 - evaluate d1H =dH/drho given db/drho in d1b 
    ##                    (2 is for evaluation of diagonal only)
    ##              4 -  given d1b and d2b evaluate trHid2H= tr(Hp^{-1}d2H/drhodrho')
    ## Hp is the preconditioned penalized hessian of the log lik
    ##    which is of rank 'rank'.
    ## fh is a factorization of Hp - either its eigen decomposition 
    ##    or its Choleski factor
    ## D is the diagonal pre-conditioning matrix used to obtain Hp
    ##   if Hr is the raw Hp then Hp = D*t(D*Hr)

      ##tr <- sort(unique(y),decreasing=TRUE)
      tr <- unique(y)
      r <- match(y,tr)
      p <- ncol(X)
      deriv <- deriv - 1
      mu <- X%*%coef
      g <- rep(0,p);H <- rep(0,p*p)
      if (deriv > 0) {
        M <- ncol(d1b)
        d1H <- if (deriv==1) rep(0,p*M) else rep(0,p*p*M)
      } else M <- d1H <- 0
      if (deriv > 2) {
        d2H <- rep(0,p*M*(M+1)/2)
        if (is.list(fh)) {
          ev <- fh
        } else  { ## need to compute eigen here
          ev <- eigen(Hp,symmetric=TRUE)
          if (rank < p) ev$values[(rank+1):p] <- 0
        } 
        X <- X%*%(ev$vectors*D)
        d1b <- t(ev$vectors)%*%(d1b/D); d2b <- t(ev$vectors)%*%(d2b/D)
      } else trHid2H <- d2H <- 0
      ## note that the following call can not use .C(C_coxlpl,...) since the ll
      ## function is not in the mgcv namespace.
      oo <- .C("coxlpl",as.double(mu),as.double(X),as.integer(r),as.integer(wt),
            as.double(tr),n=as.integer(length(y)),p=as.integer(p),nt=as.integer(length(tr)),
            lp=as.double(0),g=as.double(g),H=as.double(H),
            d1b=as.double(d1b),d1H=as.double(d1H),d2b=as.double(d2b),d2H=as.double(d2H),
            n.sp=as.integer(M),deriv=as.integer(deriv),PACKAGE="mgcv");
      if (deriv==1) d1H <- matrix(oo$d1H,p,M) else
      if (deriv>1) {
        ind <- 1:(p^2)
        d1H <- list()
        for (i in 1:M) { 
          d1H[[i]] <- matrix(oo$d1H[ind],p,p)
          ind <- ind + p^2
        }
      } 
      if (deriv>2) { 
        d2H <- matrix(oo$d2H,p,M*(M+1)/2)
        #trHid2H <- colSums(d2H)
        d <- ev$values
        d[d>0] <- 1/d[d>0];d[d<=0] <- 0
        trHid2H <- colSums(d2H*d)
      }
      assign(".log.partial.likelihood", oo$lp, envir=environment(sys.function()))
      list(l=oo$lp,lb=oo$g,lbb=matrix(oo$H,p,p),d1H=d1H,d2H=d2H,trHid2H=trHid2H)
    }

    # environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
    # environment(rd)<- environment(qf)<- environment(variance) <- environment(putTheta) 
    #environment(aic) <- 
    environment(ll) <- env
    structure(list(family = "Cox PH", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, ll=ll, ## aic = aic, 
        mu.eta = stats$mu.eta, 
        initialize = initialize,preinitialize=preinitialize,postproc=postproc,
        hazard=hazard,predict=predict,residuals=residuals,
        validmu = validmu, valideta = stats$valideta, 
        rd=rd,qf=qf,drop.intercept = TRUE,
        ls=1, ## signal ls not needed
        available.derivs = 2 ## can use full Newton here
        ),
        class = c("general.family","extended.family","family"))
} ## cox.ph

