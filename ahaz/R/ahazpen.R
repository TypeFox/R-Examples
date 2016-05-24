"[.ahazpen"<-function(object,i)
  {
    ## Purpose: Extract ahazpen object at i'th lambda-value
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   object        : ahazpen object
    ##   i             : lambda index
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    if(missing(i))
      return(object)
    if(any(i==0) || any(abs(i)>length(object$lambda)))
      stop("out of range")
    
    small <- object
    small$beta <- object$beta[i,]
    small$lambda <- object$lambda[i]
    small$df <- object$df[i]

    return(small)
  }

"ahazpen" <- function(surv, X, weights,  standardize = TRUE,  penalty = lasso.control(),
                      nlambda = 100, dfmax = nvars, pmax = min(nvars, 2*dfmax),
                      lambda.minf = ifelse(nobs < nvars,0.05, 1e-4), lambda,
                      penalty.wgt = NULL, keep = NULL, control=list()){
  ## Purpose: Penalized semiparametric additive hazards regression
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   surv          : Surv object (right censored/counting process)
  ##   X             : Numeric matrix of covariates
  ##   weights       : Weight vector (nonnegative)
  ##   standardize   : Standardize X?
  ##   penalty       : Control function for initializing penalty
  ##   nlambda       : Size of lambda grid
  ##   dfmax         : Maximal number of covariates to include in path
  ##   pmax          : Maximal number of covariates to ever consider in CCD
  ##   lambda        : Optional lambda vector
  ##   lambda.minf   : Optional minimum lambda value
  ##   penalty.wgt   : Optional penalty factor for each covariate
  ##   keep          : Keep covariates in 'keep' (no penalization) 
  ##   control       : List of control parameters for CCD algorithm
  ## ----------------------------------------------------------------------
  ## Author: Anders Gorst-Rasmussen

  
  this.call<-match.call()
  control <- do.call("ahazpen.fit.control",control)
  
  penalty<-eval(penalty)
  if(is.character(penalty)){
    tmp<-c("lasso.control","sscad.control",penalty)[pmatch(tolower(penalty),c("lasso.control","sscad.control"),nomatch = 3)]
    penalty<-get(tmp,mode="function")
  } 
  if(is.function(penalty)) penalty <- penalty()
  if(is.null(penalty$type)) {
    print(penalty)
    stop("'penalty' not recognized")
  }
  
  
  # Formatting of data, checks etc
  data <- ahaz.readnew(surv=surv, X=X, weights=weights, standardize=standardize)
  nobs <-data$nobs
  nvars <- data$nvars
  
  # Validity checks, more formatting
  if(!missing(lambda) && (!is.numeric(lambda)||any(lambda <= 0))){stop("'lambda' must be numeric and strictly positive")}
  if(nlambda<1){stop("'nlambda' must be strictly positive")}
  if(any(penalty.wgt < 0)){stop("'penalty.wgt' must be nonnegative")}
  if (dfmax < 1) { warning("'dfmax' must be strictly positive")}
  if(lambda.minf <= 0){stop("'lambda.minf' must be strictly positive")}
  if(any(keep > nvars | keep <= 0)){stop("'keep' incorrectly specified")}
  if(is.null(penalty.wgt))
    penalty.wgt <- rep(1, nvars)
  penalize<-rep(1,nvars)
  if(!is.null(keep))
    penalize[keep] <- 0
  if(missing(lambda)) {
    lambda <- numeric(nlambda)
  } else {
    lambda <- lambda
      nlambda <- length(lambda)
    lambda.minf <- 2
  }
  
  # Format penalties, do scalings if needed etc.
  sdw<-function(x,w){ wt<-w/sum(w); sqrt(sum(x^2*wt)-sum(x*wt)^2)}
  if(standardize && (is.function(penalty$init.sol)| is.function(penalty$ada.wgt)))
    sc<-apply(X,2,function(x){sdw(x,data$weights[1:data$nobs])})
  else
    sc<-rep(1,ncol(X))
  # SSCAD
  if(penalty$type=="sscad") { 
    ada.wgt<-rep(1,ncol(X))
    initsol<-penalty$init.sol(surv=surv,X=X,weights=weights)*sc
    if(is.null(penalty$c))
        prefactor<-mean(ahaz(surv=surv,X=X,univariate=TRUE,weights=weights)$D/sc^2)
    else
      prefactor<-penalty$c
   # LASSO
  }    else {
    initsol<-rep(0,ncol(X))
    prefactor <- 1
    if(is.function(penalty$ada.wgt))
      ada.wgt<-penalty$ada.wgt(surv,X,weights)*sc
    else
      ada.wgt<-rep(1,ncol(X))
  }
  
  # Un-penalize variables to be kept, apply ada.factor
  penalty.wgt<-penalty.wgt*penalize*abs(ada.wgt)
  
  # Run CCD
  a <- .C("ahapen",
          "X"=as.double(data$X),
          "time1"=as.double(data$time1),
          "time2"=as.double(data$time2),
          "event"=as.integer(data$event),
          "weights"=as.double(data$weights),
          "n"=as.integer(nobs),
          "p"=as.integer(nvars),
          "standardize"=as.integer(standardize),
          "lambdaminf"=as.double(lambda.minf),
          "nlambda"=as.integer(nlambda),
          "lambda"=as.double(lambda*sum(data$weights[1:nobs])),
          "thresh"=as.double(control$thresh),
          "maxit"=as.integer(control$maxit),
          "estims"=numeric(nlambda*pmax),
          "dfmax"=as.integer(dfmax),
          "pmax"=as.integer(pmax),
          "lambdaflag"=as.integer(nlambda),
          "error"=as.integer(0),
          "activeidx"=as.integer(rep(0,pmax)),
          "rightcens"=as.integer(data$rightcens),
          "penalty"=as.double(ifelse(penalty.wgt==Inf,-1,penalty.wgt)),
          "alpha" = as.double(penalty$alpha),
          "initsol" = as.double(initsol),
          "a" = as.double(penalty$a),
          "passesleft" = as.integer(0),
          "nsteps" = as.integer(penalty$nsteps),
          "prefactor" = as.double(prefactor))
  
  if(a$lambdaflag<=0)
    stop("no nonzero coefficients")
  if(a$error!=0)
    ahazpen.warning(a$error,a$lambdaflag)
  
  # Format dgCMatrix
  active <- (a$activeidx!=0)
  ix <- order(a$activeidx[active])
  z <- rep((1:sum(active))[ix],a$lambdaflag)+rep(a$pmax*(0:(a$lambdaflag-1)),each=sum(active))
  bt <- as.vector(a$estims[z])  
  beta <- new("dgCMatrix",Dim=as.integer(c(a$p,a$lambdaflag)),x=bt,
              p=as.integer(c(0,cumsum(rep(sum(active),a$lambdaflag)))),
              i=as.integer(rep(a$activeidx[active][ix]-1,a$lambdaflag)))
  
  # More output formatting
  df <- apply(beta,2,function(x) sum(x!=0))
  npasses<-control$maxit - a$passesleft
  lambdaout<-a$lambda[1:a$lambdaflag] / sum(data$weights)
  
  return(structure(list("call"=this.call,"beta"=beta,"lambda"=lambdaout, "df"=df, "nvars"=data$nvars,"nobs"=data$nobs,"surv"=data$surv,
                        "npasses"=npasses,"penalty.wgt"=penalty.wgt, "penalty"=penalty,"pmax"=pmax,"dfmax"=dfmax),class="ahazpen"))
}


"ahazpen.warning"<-function(error,lambdaflag)
  {
    if(error==10)
      msg<-paste("Ran out of iterations at ", lambdaflag,"th lambda value - estimates may not have converged",sep="")
    if(error==20)
      msg<-paste("Number of active variables excceeds pmax at ", lambdaflag,"th lambda value - solutions for larger lambda values returned",sep="")
    if(error==30)
      msg<-paste("Number of active variables excceeds nobs-1 at ", lambdaflag,"th lambda value - solutions for larger lambda values returned",sep="")
    warning(msg,call.=FALSE)
}
