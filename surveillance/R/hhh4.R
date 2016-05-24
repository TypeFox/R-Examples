################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### hhh4 is an extended version of algo.hhh for the sts-class
### The function allows the incorporation of random effects and covariates.
###
### Copyright (C) 2010-2012 Michaela Paul, 2012-2015 Sebastian Meyer
### $Revision: 1475 $
### $Date: 2015-09-14 22:38:55 +0200 (Mon, 14. Sep 2015) $
################################################################################

## Error message issued in loglik, score and fisher functions upon NA parameters
ADVICEONERROR <- "\n  Try different starting values, more iterations, or another optimizer.\n"


### Main function to be called by the user

hhh4 <- function (stsObj, control = list(
    ar = list(f = ~ -1,        # a formula "exp(x'lamba)*y_t-lag" (ToDo: matrix)
              offset = 1,      # multiplicative offset
              lag = 1),        # autoregression on y_i,t-lag
    ne = list(f = ~ -1,        # a formula "exp(x'phi) * sum_j w_ji * y_j,t-lag"
              offset = 1,      # multiplicative offset
              lag = 1,         # regression on y_j,t-lag
              weights = neighbourhood(stsObj) == 1,  # weights w_ji
              scale = NULL,    # such that w_ji = scale * weights
              normalize = FALSE), # w_ji -> w_ji / rowSums(w_ji), after scaling
    end = list(f = ~ 1,        # a formula "exp(x'nu) * n_it"
               offset = 1),    # optional multiplicative offset e_it
    family = c("Poisson", "NegBin1", "NegBinM"), # or a factor of length nUnit
    subset = 2:nrow(stsObj),   # epidemic components require Y_{t-lag}
    optimizer = list(stop = list(tol = 1e-5, niter = 100), # control arguments
                     regression = list(method = "nlminb"), # for optimization 
                     variance = list(method = "nlminb")),  # <- or "Nelder-Mead"
    verbose = FALSE,           # level of reporting during optimization
    start = list(fixed = NULL, # list of start values, replacing initial
                 random = NULL,  # values from fe() and ri() in 'f'ormulae
                 sd.corr = NULL), 
    data = list(t = stsObj@epoch - min(stsObj@epoch)), # named list of covariates
    keep.terms = FALSE  # whether to keep interpretControl(control, stsObj)
    ), check.analyticals = FALSE)
{
  ptm <- proc.time()
  
  ## Convert old disProg class to new sts class
  if (inherits(stsObj, "disProg")) {
      stsObj <- disProg2sts(stsObj)
  } else {
      stopifnot(inherits(stsObj, "sts"))
  }
  
  ## check control and set default values (for missing arguments)
  control <- setControl(control, stsObj)
  
  ## get model terms
  model <- interpretControl(control, stsObj)
  dimFixedEffects <- model$nFE + model$nd + model$nOverdisp
  dimRandomEffects <- model$nRE

  ## starting values 
  #* -> better default values possible
  theta.start <- model$initialTheta
  Sigma.start <- model$initialSigma
  
  ## check if initial values are valid
  ## CAVE: there might be NA's in mu if there are missing values in Y
  mu <- meanHHH(theta.start, model, total.only=TRUE)
  if(any(mu==0, na.rm=TRUE) || any(is.infinite(mu)))
    stop("some mean is degenerate (0 or Inf) at initial values")

  ## check score vector and fisher information at starting values
  check.analyticals <- if (isTRUE(check.analyticals)) {
      if (length(theta.start) > 50) "maxLik" else "numDeriv"
  } else if (is.character(check.analyticals)) {
      match.arg(check.analyticals, c("numDeriv", "maxLik"), several.ok=TRUE)
  } else NULL
  if (length(check.analyticals) > 0L) {
      resCheck <- checkAnalyticals(model, theta.start, Sigma.start,
                                   methods=check.analyticals)
      return(resCheck)
  }

  ## maximize loglikelihood (penalized and marginal)
  myoptim <- fitHHH(theta=theta.start,sd.corr=Sigma.start, model=model,
                    cntrl.stop       = control$optimizer$stop,
                    cntrl.regression = control$optimizer$regression,
                    cntrl.variance   = control$optimizer$variance,
                    verbose=control$verbose)

  convergence <- myoptim$convergence == 0
  thetahat <- c(myoptim$fixef, myoptim$ranef)
  loglik <- myoptim$loglik
  cov <- try(solve(myoptim$fisher), silent=TRUE)

  ## check for degenerate fisher info
  if(inherits(cov, "try-error")){ # fisher info is singular
    if (control$verbose)
        cat("WARNING: Final Fisher information matrix is singular!\n")
    convergence <- FALSE
  } else if(any(!is.finite(diag(cov))) || any(diag(cov)<0)){
    if (control$verbose)
        cat("WARNING: non-finite or negative covariance of regression parameters!\n")
    convergence <- FALSE
  }
  if (!convergence) {
      if (control$verbose) {
          cat("Penalized loglikelihood =", loglik, "\n")
          thetastring <- paste(round(thetahat,2), collapse=", ")
          thetastring <- strwrap(thetastring, exdent=10, prefix="\n", initial="")
          cat("theta = (", thetastring, ")\n")
      }
      warning("Results are not reliable!", ADVICEONERROR)
      res <- myoptim
      res$convergence <- convergence
      res$call <- match.call()
      class(res) <- "hhh4"
      return(res)
  }

  ## optimization successful, return a full "hhh4" object
  dimnames(cov) <- list(names(thetahat), names(thetahat))
  if (dimRandomEffects>0) {
      Sigma.orig <- myoptim$sd.corr
      Sigma.cov <- solve(myoptim$fisherVar)
      dimnames(Sigma.cov) <- list(names(Sigma.orig),names(Sigma.orig))
      Sigma.trans <- getSigmai(head(Sigma.orig,model$nVar),
                               tail(Sigma.orig,model$nCorr),
                               model$nVar)
      dimnames(Sigma.trans) <-
          rep.int(list(sub("^sd\\.", "",
                           names(Sigma.orig)[seq_len(model$nVar)])), 2L)
  } else {
      Sigma.orig <- Sigma.cov <- Sigma.trans <- NULL
  }
  
  ## Done
  result <- list(coefficients=thetahat, se=sqrt(diag(cov)), cov=cov, 
                 Sigma=Sigma.trans,     # estimated covariance matrix of ri's
                 Sigma.orig=Sigma.orig, # variance parameters on original scale
                 Sigma.cov=Sigma.cov,   # covariance matrix of Sigma.orig
                 call=match.call(),
                 dim=c(fixed=dimFixedEffects,random=dimRandomEffects),
                 loglikelihood=loglik,
                 margll=marLogLik(Sigma.orig, thetahat, model),
                 convergence=convergence,
                 fitted.values=meanHHH(thetahat, model, total.only=TRUE),
                 control=control,
                 terms=if(control$keep.terms) model else NULL,
                 stsObj=stsObj,
                 lags=sapply(control[c("ar","ne")], function (comp)
                             if (comp$inModel) comp$lag else NA_integer_),
                 nObs=sum(!model$isNA[control$subset,]),
                 nTime=length(model$subset), nUnit=ncol(stsObj),
                 ## CAVE: nTime is not nrow(stsObj) as usual!
                 runtime=proc.time()-ptm)
  class(result) <- "hhh4"
  result
}


## set default values for model specifications in control
setControl <- function (control, stsObj)
{
  stopifnot(is.list(control))
  nTime <- nrow(stsObj)
  nUnit <- ncol(stsObj)
  if(nTime <= 2) stop("too few observations")
  
  ## arguments in 'control' override any corresponding default arguments
  defaultControl <- eval(formals(hhh4)$control)
  environment(defaultControl$ar$f) <- environment(defaultControl$ne$f) <-
      environment(defaultControl$end$f) <- .GlobalEnv
  control <- modifyList(defaultControl, control)

  ## check that component specifications are list objects
  for (comp in c("ar", "ne", "end")) {
      if(!is.list(control[[comp]])) stop("'control$", comp, "' must be a list")
  }

  ## check lags in "ar" and "ne" components
  for (comp in c("ar", "ne")) {
      if (!isScalar(control[[comp]]$lag) || control[[comp]]$lag < (comp=="ar"))
          stop("'control$", comp, "$lag' must be a ",
               if (comp=="ar") "positive" else "non-negative", " integer")
      control[[comp]]$lag <- as.integer(control[[comp]]$lag)
  }

  
  ### check AutoRegressive component

  if (control$ar$isMatrix <- is.matrix(control$ar$f)) {
      ## this form is not implemented -> will stop() in interpretControl()
      if (any(dim(control$ar$f) != nUnit))
          stop("'control$ar$f' must be a square matrix of size ", nUnit)
      if (is.null(control$ar$weights)) { # use identity matrix
          control$ar$weights <- diag(nrow=nUnit)
      } else if (!is.matrix(control$ar$weights) ||
                 any(dim(control$ar$weights) != nUnit)) {
          stop("'control$ar$weights' must be a square matrix of size ", nUnit)
      }
      control$ar$inModel <- TRUE
  } else if (inherits(control$ar$f, "formula")) {
      if (!is.null(control$ar$weights)) {
          warning("argument 'control$ar$weights' is not used")
          control$ar$weights <- NULL
      }
      # check if formula is valid
      control$ar$inModel <- isInModel(control$ar$f)
  } else {
      stop("'control$ar$f' must be either a formula or a matrix")
  }
  
  
  ### check NEighbourhood component
  
  if (!inherits(control$ne$f, "formula"))
      stop("'control$ne$f' must be a formula")
  control$ne$inModel <- isInModel(control$ne$f)
  
  if (control$ne$inModel) {
      if (nUnit == 1)
          stop("\"ne\" component requires a multivariate 'stsObj'")
      ## if ar$f is a matrix it includes neighbouring units => no "ne" component
      if (control$ar$isMatrix)
          stop("there must not be an extra \"ne\" component ",
               "if 'control$ar$f' is a matrix")
      ## check ne$weights specification
      checkWeights(control$ne$weights, nUnit, nTime,
                   neighbourhood(stsObj), control$data,
                   check0diag = control$ar$inModel)
      ## check optional scaling of weights
      if (!is.null(control$ne$scale)) {
          stopifnot(is.numeric(control$ne$scale))
          if (is.vector(control$ne$scale)) {
              stopifnot(length(control$ne$scale) == 1L ||
                            length(control$ne$scale) %% nUnit == 0,
                        !is.na(control$ne$scale))
          } else {
              checkWeightsArray(control$ne$scale, nUnit, nTime)
          }
      }
  } else {
      control$ne[c("weights", "scale", "normalize")] <- list(NULL, NULL, FALSE)
  }

  
  ### check ENDemic component
  
  if (!inherits(control$end$f, "formula"))
      stop("'control$end$f' must be a formula")
  control$end$inModel <- isInModel(control$end$f)


  ### check offsets

  for (comp in c("ar", "ne", "end")) {
      if (is.matrix(control[[comp]]$offset) && is.numeric(control[[comp]]$offset)){
          if (!identical(dim(control[[comp]]$offset), dim(stsObj)))
              stop("'control$",comp,"$offset' must be a numeric matrix of size ",
                   nTime, "x", nUnit)
          if (any(is.na(control[[comp]]$offset)))
              stop("'control$",comp,"$offset' must not contain NA values")
      } else if (!identical(as.numeric(control[[comp]]$offset), 1)) {
          stop("'control$",comp,"$offset' must either be 1 or a numeric ",
               nTime, "x", nUnit, " matrix")
      }
  }


  ### stop if no component is included in the model
  
  if (length(comps <- componentsHHH4(list(control=control))) == 0L)
      stop("none of the components 'ar', 'ne', 'end' is included in the model")
  

  ### check remaining components of the control list

  if (is.factor(control$family)) {
      stopifnot(length(control$family) == nUnit)
      control$family <- droplevels(control$family)
      names(control$family) <- colnames(stsObj)
  } else {
      control$family <- match.arg(control$family, defaultControl$family)
  }

  if (!is.vector(control$subset, mode="numeric") ||
      !all(control$subset %in% seq_len(nTime)))
      stop("'control$subset' must be %in% 1:", nTime)
  lags <- c(ar = control$ar$lag, ne = control$ne$lag)
  maxlag <- suppressWarnings(max(lags[names(lags) %in% comps])) # could be -Inf
  if (control$subset[1L] <= maxlag) {
      warning("'control$subset' should be > ", maxlag, " due to epidemic lags")
  }

  if (!is.list(control$optimizer) ||
      any(! sapply(c("stop", "regression", "variance"),
                   function(x) is.list(control$optimizer[[x]]))))
      stop("'control$optimizer' must be a list of lists")

  control$verbose <- as.integer(control$verbose)
  if (length(control$verbose) != 1L || control$verbose < 0)
      stop("'control$verbose' must be a logical or non-negative numeric value")

  stopifnot(is.list(control$start))
  control$start <- local({
      defaultControl$start[] <- control$start[names(defaultControl$start)]
      defaultControl$start
  })
  if (!all(vapply(X = control$start,
                  FUN = function(x) is.null(x) || is.vector(x, mode="numeric"),
                  FUN.VALUE = TRUE, USE.NAMES = FALSE)))
      stop("'control$start' must be a list of numeric start values")
  
  stopifnot(length(control$keep.terms) == 1L, is.logical(control$keep.terms))

  ## Done
  return(control)
}


# check whether or not one of the three components is included in the model
isInModel <- function(formula, name=deparse(substitute(formula)))
{
  term <- terms.formula(formula)
  if(attr(term,"response") > 0) stop(name, " cannot contain a response")
  attr(term, "intercept") + length(attr(term, "term.labels")) > 0
}

# used to incorporate covariates and unit-specific effects
fe <- function(x,          # covariate
               unitSpecific = FALSE, # TRUE means which = rep.int(TRUE, nUnits)
               which=NULL, # NULL = overall, vector with booleans = unit-specific
               initial=NULL) # vector of inital values for parameters
{
  stsObj <- get("stsObj", envir=parent.frame(1), inherits=TRUE) #checkFormula()
  nTime <- nrow(stsObj)
  nUnits <- ncol(stsObj)
  
  if(!is.numeric(x)){
    stop("Covariate \'",deparse(substitute(x)),"\' is not numeric\n")
  }
  lengthX <- length(x)
  if(lengthX == 1){
    terms <- matrix(x, nTime, nUnits, byrow=FALSE)
    mult <- "*"
  } else if(lengthX == nTime){
    terms <- matrix(x, nTime, nUnits, byrow=FALSE)
    mult <- "*"  
  } else if(lengthX == nTime*nUnits){
    if(!is.matrix(x)){
     stop("Covariate \'",deparse(substitute(x)),"\' is not a matrix\n")    
    }
    # check dimensions of covariate
    if((ncol(x) != nUnits) | (nrow(x) != nTime)){
     stop("Dimension of covariate \'",deparse(substitute(x)),"\' is not suitably specified\n")    
    }
    terms <- x
    mult <- "*"  
  } else {
    stop("Covariate \'",deparse(substitute(x)),"\' is not suitably specified\n")
  }
  
  intercept <- all(terms==1)
  
  # overall or unit-specific effect?
  unitSpecific <- unitSpecific || !is.null(which)
  if (unitSpecific) {
      if (is.null(which)) {
          which <- rep.int(TRUE, nUnits)
      } else {
          stopifnot(is.vector(which, mode="logical"), length(which) == nUnits)
      }
      
      terms[,!which] <- 0
  }
  
  # get dimension of parameter
  dim.fe <- if (unitSpecific) sum(which) else 1
  
  # check length of initial values + set default values
  if (is.null(initial)) {
    initial <- rep.int(0,dim.fe)
  } else if (length(initial) != dim.fe) {
    stop("initial values for '",deparse(substitute(x)),"' must be of length ",dim.fe)
  }
  
  summ <- if (unitSpecific) "colSums" else "sum"
    
  name <- deparse(substitute(x))
  if (unitSpecific)
      name <- paste(name, colnames(stsObj)[which], sep=".")
    
  result <- list(terms=terms,
                name=name,
                Z.intercept=NULL,
                which=which,
                dim.fe=dim.fe,
                initial.fe=initial,
                dim.re=0,
                dim.var=0,
                initial.var=NULL,
                initial.re=NULL,
                intercept=intercept,
                unitSpecific=unitSpecific,
                random=FALSE,
                corr=FALSE,
                summ=summ,
                mult=mult
                )
  return(result)
}

# random intercepts
ri <- function(type=c("iid","car"), 
               corr=c("none","all"),
               initial.fe=0,
               initial.var=-.5,
               initial.re=NULL)
{
  stsObj <- get("stsObj", envir=parent.frame(1), inherits=TRUE) #checkFormula()
  if (ncol(stsObj) == 1)
      stop("random intercepts require a multivariate 'stsObj'")
  type <- match.arg(type)
  corr <- match.arg(corr)
  corr <- switch(corr, 
                  "none"=FALSE,
                  "all"=TRUE)
                  
  if(type=="iid"){
    Z <- 1
    dim.re <- ncol(stsObj)
    mult <- "*"
  } else if(type=="car"){ # construct penalty matrix K
    K <- neighbourhood(stsObj)
    checkNeighbourhood(K)
    K <- K == 1                         # indicate first-order neighbours
    ne <- colSums(K)                    # number of first-order neighbours
    K <- -1*K
    diag(K) <- ne
    
    dimK <- nrow(K)
    
    # check rank of the nhood, only connected neighbourhoods are allowed
    if(qr(K)$rank != dimK-1) stop("neighbourhood matrix contains islands")
    # singular-value decomposition of K
    svdK <- svd(K)
    # just use the positive eigenvalues  of K in descending order 
    # for a the factorisation of the penalty matrix K = LL'
    L <- svdK$u[,-dimK] %*% diag(sqrt(svdK$d[-dimK]))            #* only use non-zero eigenvalues
    
    # Z = L(L'L)^-1, which can't be simplified to Z=(L')^-1 as L is not square
    Z <- L %*% solve(t(L)%*%L)
    
    dim.re <- dimK - 1L
    mult <- "%*%"
  } 
  
  # check length of initial values + set default values
  stopifnot(length(initial.fe) == 1, length(initial.var) == 1)
  if (is.null(initial.re)) {
    initial.re <- rnorm(dim.re,0,sd=sqrt(0.001))
  } else if (length(initial.re) != dim.re) {
    stop("'initial.re' must be of length ", dim.re)
  }

  result <- list(terms=1,
                name=paste("ri(",type,")",sep=""),
                Z.intercept=Z,
                which=NULL,
                dim.fe=1,
                initial.fe=initial.fe,
                dim.re=dim.re,
                dim.var=1,
                initial.var=initial.var,
                initial.re=initial.re,
                intercept=TRUE,
                unitSpecific=FALSE,
                random=TRUE,
                corr=corr,
                summ="colSums",
                mult=mult
                )
  return(result)
}

### check specification of formula
## f: one of the component formulae (ar$f, ne$f, or end$f)
## component: 1, 2, or 3, corresponding to the ar/ne/end component, respectively
## data: the data-argument of hhh4()
## stsObj: the stsObj is not used directly in checkFormula, but in fe() and ri()
checkFormula <- function(f, component, data, stsObj)
{
  term <- terms.formula(f, specials=c("fe","ri"))
  
  # check if there is an overall intercept
  intercept.all <- attr(term, "intercept") == 1
  
  # list of variables in the component
  vars <- as.list(attr(term,"variables"))[-1] # first element is "list"
  nVars <- length(vars)

  # begin with intercept
  res <- if (intercept.all) {
      c(fe(1), list(offsetComp=component))
  } else {
      if (nVars==0)
          stop("formula ", deparse(substitute(f)), " contains no variables")
      NULL
  }
  
  # find out fixed effects without "fe()" specification
  # (only if there are variables in addition to an intercept "1")
  fe.raw <- setdiff(seq_len(nVars), unlist(attr(term, "specials")))
  
  # evaluate covariates
  for(i in fe.raw)
      res <- cbind(res, c(
          eval(substitute(fe(x), list(x=vars[[i]])), envir=data),
          list(offsetComp=component)
          ))
  
  # fixed effects
  for(i in attr(term, "specials")$fe)
      res <- cbind(res, c(
          eval(vars[[i]], envir=data),
          list(offsetComp=component)
          ))
  
  res <- cbind(res, deparse.level=0) # ensure res has matrix dimensions
  
  # random intercepts
  RI <- attr(term, "specials")$ri
  if (sum(unlist(res["intercept",])) + length(RI) > 1)
      stop("There can only be one intercept in the formula ",
           deparse(substitute(f)))
  for(i in RI)
      res <- cbind(res, c(
          eval(vars[[i]], envir=data),
          list(offsetComp=component)
          ))
  
  return(res)
}


## Create function (pars, type = "response") which
## returns the weighted sum of time-lagged counts of neighbours
## (or its derivates, if type = "gradient" or type = "hessian").
## For type="reponse", this is a nTime x nUnits matrix (like Y),
## otherwise a list of such matrices,
## which for the gradient has length length(pars) and
## length(pars)*(length(pars)+1)/2 for the hessian.
## If neweights=NULL (i.e. no NE component in model), the result is always 0.
## offset is a multiplicative offset for \phi_{it}, e.g., the population.
## scale is a nUnit-vector or a nUnit x nUnit matrix scaling neweights.
neOffsetFUN <- function (Y, neweights, scale, normalize,
                         nbmat, data, lag = 1, offset = 1)
{
    if (is.null(neweights)) { # no neighbourhood component
        as.function(alist(...=, 0), envir=.GlobalEnv)
        ## dimY <- dim(Y)
        ## as.function(c(alist(...=),
        ##               substitute(matrix(0, r, c), list(r=dimY[1], c=dimY[2]))),
        ##             envir=.GlobalEnv)
    } else if (is.list(neweights)) { # parametric weights
        wFUN <- scaleNEweights.list(neweights, scale, normalize)
        function (pars, type = "response") {
            name <- switch(type, response="w", gradient="dw", hessian="d2w")
            weights <- wFUN[[name]](pars, nbmat, data)
            ## gradient and hessian are lists if length(pars$d) > 1L
            ## but can be single matrices/arrays if == 1 => _c_onditional lapply
            res <- clapply(weights, function (W)
                           offset * weightedSumNE(Y, W, lag))
            ##<- clapply always returns a list (possibly of length 1)
            if (type=="response") res[[1L]] else res
        }
    } else { # fixed (known) weight structure (0-length pars)
        weights <- scaleNEweights.default(neweights, scale, normalize)
        env <- new.env(hash = FALSE, parent = emptyenv())  # small -> no hash
        env$initoffset <- offset * weightedSumNE(Y, weights, lag)
        as.function(c(alist(...=), quote(initoffset)), envir=env)
    }
}


# interpret and check the specifications of each component
# control must contain all arguments, i.e. setControl was used
interpretControl <- function (control, stsObj)
{
  nTime <- nrow(stsObj)
  nUnits <- ncol(stsObj)

  Y <- observed(stsObj)

  
  ##########################################################################
  ##  get the model specifications for each of the three components
  ##########################################################################
  ar <- control$ar
  ne <- control$ne
  end <- control$end

  ## for backwards compatibility with surveillance < 1.8-0, where the ar and ne
  ## components of the control object did not have an offset
  if (is.null(ar$offset)) ar$offset <- 1
  if (is.null(ne$offset)) ne$offset <- 1
  ## for backward compatibility with surveillance < 1.9-0
  if (is.null(ne$normalize)) ne$normalize <- FALSE
  
  ## create list of offsets of the three components
  Ym1 <- rbind(matrix(NA_integer_, ar$lag, nUnits), head(Y, nTime-ar$lag))
  Ym1.ne <- neOffsetFUN(Y, ne$weights, ne$scale, ne$normalize,
                        neighbourhood(stsObj), control$data, ne$lag, ne$offset)
  offsets <- list(ar=ar$offset*Ym1, ne=Ym1.ne, end=end$offset)
  ## -> offset$ne is a function of the parameter vector 'd', which returns a
  ##    nTime x nUnits matrix -- or 0 (scalar) if there is no NE component
  ## -> offset$end might just be 1 (scalar)
  
  ## Initial parameter vector 'd' of the neighbourhood weight function
  initial.d <- if (is.list(ne$weights)) ne$weights$initial else numeric(0L)
  dim.d <- length(initial.d)
  names.d <- if (dim.d == 0L) character(0L) else {
      paste0("neweights.", if (is.null(names(initial.d))) {
          if (dim.d==1L) "d" else paste0("d", seq_len(dim.d))
      } else names(initial.d))
  }

  ## determine all NA's (FIXME: why do we need this? Why include is.na(Y)?)
  isNA <- is.na(Y)
  if (ar$inModel)
      isNA <- isNA | is.na(offsets[[1L]])
  if (ne$inModel)
      isNA <- isNA | is.na(offsets[[2L]](initial.d))

  ## get terms for all components
  all.term <- NULL
  if(ar$isMatrix) stop("matrix-form of 'control$ar$f' is not implemented")
  if(ar$inModel) # ar$f is a formula
      all.term <- cbind(all.term, checkFormula(ar$f, 1, control$data, stsObj))
  if(ne$inModel)
      all.term <- cbind(all.term, checkFormula(ne$f, 2, control$data, stsObj))
  if(end$inModel)
      all.term <- cbind(all.term, checkFormula(end$f,3, control$data, stsObj))
  
  dim.fe <- sum(unlist(all.term["dim.fe",]))
  dim.re.group <- unlist(all.term["dim.re",], use.names=FALSE)
  dim.re <- sum(dim.re.group)
  dim.var <- sum(unlist(all.term["dim.var",]))
  dim.corr <- sum(unlist(all.term["corr",]))
  
  if(dim.corr>0){
    if(dim.var!=dim.corr) stop("Use corr=\'all\' or corr=\'none\' ")
    dim.corr <- switch(dim.corr,0,1,3) 
  }
  
  # the vector with dims of the random effects must be equal if they are correlated
  if(length(unique(dim.re.group[dim.re.group>0]))!=1 & dim.corr>0){ 
    stop("Correlated effects must have same penalty")
  }
  
  n <- c("ar","ne","end")[unlist(all.term["offsetComp",])]
  names.fe <- names.var <- names.re <- character(0L)
  for(i in seq_along(n)){
    .name <- all.term["name",i][[1]]
    names.fe <- c(names.fe, paste(n[i], .name, sep="."))
    if(all.term["random",i][[1]]) {
        names.var <- c(names.var, paste("sd", n[i], .name, sep="."))
        names.re <- c(names.re, paste(n[i], .name, if (.name == "ri(iid)") {
            colnames(stsObj)
        } else {
            seq_len(all.term["dim.re",i][[1]])
        }, sep = "."))
    }
  }
  index.fe <- rep(1:ncol(all.term), times=unlist(all.term["dim.fe",]))
  index.re <- rep(1:ncol(all.term), times=unlist(all.term["dim.re",])) 
  
  # poisson or negbin model
  if(identical(control$family, "Poisson")){
    ddistr <- function(y,mu,size){
        dpois(y, lambda=mu, log=TRUE)
    }
    dim.overdisp <- 0L
    index.overdisp <- names.overdisp <- NULL
  } else { # NegBin
    ddistr <- function(y,mu,size){
        dnbinom(y, mu=mu, size=size, log=TRUE)
    }
    ## version that can handle size = Inf (i.e. the Poisson special case):
    ## ddistr <- function (y,mu,size) {
    ##     poisidx <- is.infinite(size)
    ##     res <- y
    ##     res[poisidx] <- dpois(y[poisidx], lambda=mu[poisidx], log=TRUE)
    ##     res[!poisidx] <- dnbinom(y[!poisidx], mu=mu[!poisidx],
    ##                              size=size[!poisidx], log=TRUE)
    ##     res
    ## }
    index.overdisp <- if (is.factor(control$family)) {
        control$family
    } else if (control$family == "NegBinM") {
        factor(colnames(stsObj), levels = colnames(stsObj))
        ## do not sort levels (for consistency with unitSpecific effects)
    } else { # "NegBin1"
        factor(character(nUnits))
    }
    names(index.overdisp) <- colnames(stsObj)
    dim.overdisp <- nlevels(index.overdisp)
    names.overdisp <- if (dim.overdisp == 1L) {
        "-log(overdisp)"
    } else {
        paste0("-log(", paste("overdisp", levels(index.overdisp), sep = "."), ")")
    }
  }
  environment(ddistr) <- getNamespace("stats")  # function is self-contained

  # parameter start values from fe() and ri() calls via checkFormula()
  initial <- list(
      fixed = c(unlist(all.term["initial.fe",]),
                initial.d,
                rep.int(2, dim.overdisp)),
      random = as.numeric(unlist(all.term["initial.re",])), # NULL -> numeric(0)
      sd.corr = c(unlist(all.term["initial.var",]),
                  rep.int(0, dim.corr))
  )
  # set names of parameter vectors
  names(initial$fixed) <- c(names.fe, names.d, names.overdisp)
  names(initial$random) <- names.re
  names(initial$sd.corr) <- c(names.var, head(paste("corr",1:3,sep="."), dim.corr))

  # modify initial values according to the supplied 'start' values
  initial[] <- mapply(
      FUN = function (initial, start, name) {
          if (is.null(start))
              return(initial)
          if (is.null(names(initial)) || is.null(names(start))) {
              if (length(start) == length(initial)) {
                  initial[] <- start
              } else {
                  stop("initial values in 'control$start$", name,
                       "' must be of length ", length(initial))
              }
          } else {
              ## we match by name and silently ignore additional start values
              start <- start[names(start) %in% names(initial)]
              initial[names(start)] <- start
          }
          return(initial)
      },
      initial, control$start[names(initial)], names(initial),
      SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
  
  # Done
  result <- list(response = Y,
                 terms = all.term,
                 nTime = nTime,
                 nUnits = nUnits,
                 nFE = dim.fe,
                 nd = dim.d,
                 nOverdisp = dim.overdisp,
                 nRE = dim.re,
                 rankRE = dim.re.group,
                 nVar = dim.var,
                 nCorr = dim.corr,
                 nSigma = dim.var+dim.corr,
                 nGroups = ncol(all.term),
                 namesFE = names.fe,
                 indexFE = index.fe,
                 indexRE = index.re,
                 initialTheta = c(initial$fixed, initial$random),
                 initialSigma = initial$sd.corr,
                 offset = offsets,
                 family = ddistr,
                 indexPsi = index.overdisp,
                 subset = control$subset,
                 isNA = isNA
                 )
  return(result)
}


splitParams <- function(theta, model){
  fixed <- theta[seq_len(model$nFE)]
  d <- theta[model$nFE + seq_len(model$nd)]
  overdisp <- theta[model$nFE + model$nd + seq_len(model$nOverdisp)]
  random <- theta[seq.int(to=length(theta), length.out=model$nRE)]
  list(fixed=fixed, random=random, overdisp=overdisp, d=d)
}


### compute predictor

meanHHH <- function(theta, model, subset=model$subset, total.only=FALSE)
{
  ## unpack theta
  pars <- splitParams(theta, model)
  fixed <- pars$fixed
  random <- pars$random

  ## unpack model
  term <- model$terms
  offsets <- model$offset
  offsets[[2L]] <- offsets[[2L]](pars$d) # evaluate at current parameter value
  nGroups <- model$nGroups
  comp <- unlist(term["offsetComp",])
  idxFE <- model$indexFE
  idxRE <- model$indexRE
  
  #isNA <- model$isNA

  toMatrix <- function (x, r=model$nTime, c=model$nUnits)
      matrix(x, r, c, byrow=TRUE)
  
  ## go through groups of parameters and compute predictor of each component,
  ## i.e. lambda_it, phi_it, nu_it, EXCLUDING the multiplicative offset terms,
  ## as well as the resulting component mean (=exppred * offset)
  computePartMean <- function (component)
  {
    pred <- nullMatrix <- toMatrix(0)
    #is.na(pred) <- isNA # set missing values in observed Y to NA
                         # -> FIXME: why? seems to be awkward... and
                         # incompatible with simulation for t=1 (given y.start)
    
    if(!any(comp==component)) { # component not in model -> return 0-matrix
        zeroes <- pred[subset,,drop=FALSE]
        return(list(exppred = zeroes, mean = zeroes))
    }
    
    for(i in seq_len(nGroups)[comp==component]){
      fe <- fixed[idxFE==i]
      if(term["unitSpecific",i][[1]]){
        fe <- nullMatrix
        which <- term["which",i][[1]]
        fe[,which] <- toMatrix(fixed[idxFE==i],c=sum(which))
      }
      if(term["random",i][[1]]){
        re <- random[idxRE==i]
        "%m%" <- get(term["mult",i][[1]])
        Z.re <- toMatrix(term["Z.intercept",i][[1]] %m% re)
      } else {
        Z.re <- 0
      }
      X <- term["terms",i][[1]]
      pred <- pred + X*fe + Z.re
    }

    exppred <- exp(pred[subset,,drop=FALSE])
    offset <- offsets[[component]]
    if (length(offset) > 1) offset <- offset[subset,,drop=FALSE]
    ##<- no subsetting if offset is scalar (time- and unit-independent)
    list(exppred = exppred, mean = exppred * offset)
  }
  
  ## compute component means
  ar <- computePartMean(1)
  ne <- computePartMean(2)
  end <- computePartMean(3)
  
  ## Done
  epidemic <- ar$mean + ne$mean
  endemic <- end$mean
  if (total.only) epidemic + endemic else
  list(mean=epidemic+endemic, epidemic=epidemic, endemic=endemic,
       epi.own=ar$mean, epi.neighbours=ne$mean,
       ar.exppred=ar$exppred, ne.exppred=ne$exppred, end.exppred=end$exppred)
}


### compute dispersion in dnbinom (mu, size) parametrization

sizeHHH <- function (theta, model, subset = model$subset)
{
    if (model$nOverdisp == 0L) # Poisson case
        return(NULL)

    ## extract dispersion in dnbinom() parametrization
    pars <- splitParams(theta, model)
    size <- exp(pars$overdisp)             # = 1/psi, pars$overdisp = -log(psi)

    ## return either a vector or a time x unit matrix of dispersion parameters
    if (is.null(subset)) {
        unname(size)  # no longer is "-log(overdisp)"
    } else {
        matrix(data = size[model$indexPsi],
               nrow = length(subset), ncol = model$nUnits, byrow = TRUE,
               dimnames = list(NULL, names(model$indexPsi)))
    }
}

## auxiliary function used in penScore and penFisher
## it sums colSums(x) within the groups defined by f (of length ncol(x))
.colSumsGrouped <- function (x, f, na.rm = TRUE)
{
    nlev <- nlevels(f)
    if (nlev == 1L) { # all columns belong to the same group
        sum(x, na.rm = na.rm)
    } else {
        dimx <- dim(x)
        colsums <- .colSums(x, dimx[1L], dimx[2L], na.rm = na.rm)
        if (nlev == dimx[2L]) { # each column is its own group
            colsums[order(f)]
        } else { # sum colsums within groups
            unlist(lapply(
                X = split.default(colsums, f, drop = FALSE),
                FUN = sum
                ), recursive = FALSE, use.names = FALSE)
        }
    }
}


############################################


penLogLik <- function(theta, sd.corr, model, attributes=FALSE)
{
  if(any(is.na(theta))) stop("NAs in regression parameters.", ADVICEONERROR)

  ## unpack model
  subset <- model$subset
  Y <- model$response[subset,,drop=FALSE]
  #isNA <- model$isNA[subset,,drop=FALSE]
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE
  
  ## unpack random effects
  if (dimRE > 0) {
      pars <- splitParams(theta, model)
      randomEffects <- pars$random
      sd   <- head(sd.corr, model$nVar)
      corr <- tail(sd.corr, model$nCorr)
      dimBlock <- model$rankRE[model$rankRE>0]
      Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)
  }
  
  ############################################################

  ## evaluate dispersion
  psi <- sizeHHH(theta, model,
                 subset = if (dimPsi > 1L) subset) # else scalar or NULL

  #psi might be numerically equal to 0 or Inf in which cases dnbinom (in meanHHH)
  #would return NaN (with a warning). The case size=Inf rarely happens and
  #corresponds to a Poisson distribution. Currently this case is not handled
  #in order to have the usual non-degenerate case operate faster.
  #For size=0, log(dnbinom) equals -Inf for positive x or if (x=0 and mu=0), and
  #zero if x=0 and mu>0 and mu<Inf. Thus, if there is at least one case in Y
  #(x>0, which is always true), we have that sum(ll.units) = -Inf, hence: 
  if (any(psi == 0)) return(-Inf)

  ## evaluate mean
  mu <- meanHHH(theta, model, total.only=TRUE)
  # if, numerically, mu=Inf, log(dnbinom) or log(dpois) both equal -Inf, hence:
  #if (any(is.infinite(mu))) return(-Inf)
  # however, since mu=Inf does not produce warnings below and this is a rare
  # case, it is faster to not include this conditional expression

  ## penalization term for random effects
  lpen <- if (dimRE==0) 0 else { # there are random effects
    ##-.5*(t(randomEffects)%*%Sigma.inv%*%randomEffects)
    ## the following implementation takes ~85% less computing time !
    -0.5 * crossprod(randomEffects, Sigma.inv) %*% randomEffects
  }
  lpen <- c(lpen)                       # drop 1x1 matrix dimensions
  
  ## log-likelihood
  ll.units <- .colSums(model$family(Y,mu,psi),
                       length(subset), model$nUnits, na.rm=TRUE)

  ## penalized log-likelihood
  ll <- sum(ll.units) + lpen

  ## Done
  if (attributes) {
      attr(ll, "loglik") <- ll.units
      attr(ll, "logpen") <- lpen
  }
  ll
}


penScore <- function(theta, sd.corr, model)
{
  if(any(is.na(theta))) stop("NAs in regression parameters.", ADVICEONERROR)

  ## unpack model
  subset <- model$subset
  Y <- model$response[subset,,drop=FALSE]
  isNA <- model$isNA[subset,,drop=FALSE]
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE
  term <- model$terms
  nGroups <- model$nGroups
  dimd <- model$nd

  ## unpack parameters
  pars <- splitParams(theta, model)
  if (dimRE > 0) {
      randomEffects <- pars$random
      sd   <- head(sd.corr, model$nVar)
      corr <- tail(sd.corr, model$nCorr)
      dimBlock <- model$rankRE[model$rankRE>0]
      Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)
  }

  ## evaluate dispersion
  psi <- sizeHHH(theta, model,
                 subset = if (dimPsi > 1L) subset) # else scalar or NULL

  ## evaluate mean
  mu <- meanHHH(theta, model)
  meanTotal <- mu$mean

  ############################################################
  
  ## helper function for derivatives
  derivHHH.factor <- if(dimPsi > 0L){ # NegBin
      psiPlusMu <- psi + meanTotal    # also used below for calculation of grPsi
      psiYpsiMu <- (psi+Y) / psiPlusMu
      Y/meanTotal - psiYpsiMu
  } else { # Poisson
      Y/meanTotal - 1
  }
  derivHHH <- function (dmu) derivHHH.factor * dmu

  ## go through groups of parameters and compute the gradient of each component
  computeGrad <- function(mean.comp){
  
    grad.fe <- numeric(0L)
    grad.re <- numeric(0L)
       
    for(i in seq_len(nGroups)){
      comp <- term["offsetComp",i][[1]]
      Xit<- term["terms",i][[1]] # eiter 1 or a matrix with values
      if(is.matrix(Xit)){
        Xit <- Xit[subset,,drop=FALSE]
      }
      summ <- get(term["summ",i][[1]])
      dTheta <- derivHHH(mean.comp[[comp]]*Xit)
      dTheta[isNA] <- 0   # dTheta must not contain NA's (set NA's to 0)
      
      if(term["unitSpecific",i][[1]]){
        which <- term["which",i][[1]]
        dTheta <- summ(dTheta)[ which ]
        grad.fe <- c(grad.fe,dTheta)
        
      } else if(term["random",i][[1]]){
        Z <- term["Z.intercept",i][[1]]  
        "%m%" <- get(term["mult",i][[1]])
        dRTheta <- .colSums(dTheta %m% Z, length(subset), term["dim.re",i][[1]])
        grad.re <- c(grad.re, dRTheta)
        grad.fe <- c(grad.fe, sum(dTheta))
      } else{
        grad.fe <- c(grad.fe, summ(dTheta))
      }
    }
    
    list(fe=grad.fe, re=grad.re)
  }
  
  gradients <- computeGrad(mu[c("epi.own","epi.neighbours","endemic")])

  ## gradient for parameter vector of the neighbourhood weights
  grd <- if (dimd > 0L) {
      dneOffset <- model$offset[[2L]](pars$d, type="gradient")
      ##<- this is always a list (of length dimd) of matrices
      onescore.d <- function (dneoff) {
          dmudd <- mu$ne.exppred * dneoff[subset,,drop=FALSE]
          grd.terms <- derivHHH(dmudd)
          sum(grd.terms, na.rm=TRUE)
      }
      unlist(clapply(dneOffset, onescore.d), recursive=FALSE, use.names=FALSE)
  } else numeric(0L)
  
  ## gradient for overdispersion parameter psi
  grPsi <- if(dimPsi > 0L){
      dPsiMat <- psi * (digamma(Y+psi) - digamma(psi) + log(psi) + 1
                        - log(psiPlusMu) - psiYpsiMu)
      .colSumsGrouped(dPsiMat, model$indexPsi)
  } else numeric(0L)
  
  ## add penalty to random effects gradient
  s.pen <- if(dimRE > 0) c(Sigma.inv %*% randomEffects) else numeric(0L)
  if(length(gradients$re) != length(s.pen))
      stop("oops... lengths of s(b) and Sigma.inv %*% b do not match")
  grRandom <- c(gradients$re - s.pen)

  ## Done
  res <- c(gradients$fe, grd, grPsi, grRandom)
  res
}


penFisher <- function(theta, sd.corr, model, attributes=FALSE)
{
  if(any(is.na(theta))) stop("NAs in regression parameters.", ADVICEONERROR)

  ## unpack model
  subset <- model$subset
  Y <- model$response[subset,,drop=FALSE]
  isNA <- model$isNA[subset,,drop=FALSE]
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE
  term <- model$terms
  nGroups <- model$nGroups
  dimd  <- model$nd
  dimFE <- model$nFE
  idxFE <- model$indexFE
  idxRE <- model$indexRE
  indexPsi <- model$indexPsi

  ## unpack parameters
  pars <- splitParams(theta, model)
  if (dimRE > 0) {
      randomEffects <- pars$random
      sd   <- head(sd.corr, model$nVar)
      corr <- tail(sd.corr, model$nCorr)
      dimBlock <- model$rankRE[model$rankRE>0]
      Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)
  }

  ## evaluate dispersion
  psi <- sizeHHH(theta, model,
                 subset = if (dimPsi > 1L) subset) # else scalar or NULL

  ## evaluate mean
  mu <- meanHHH(theta, model)
  meanTotal <- mu$mean
  
  ############################################################
  
  ## helper functions for derivatives:
  if (dimPsi > 0L) { # negbin
    psiPlusY <- psi + Y
    psiPlusMu <- psi + meanTotal
    psiPlusMu2 <- psiPlusMu^2
    psiYpsiMu <- psiPlusY / psiPlusMu
    psiYpsiMu2 <- psiPlusY / psiPlusMu2
    deriv2HHH.fac1 <- psiYpsiMu2 - Y / (meanTotal^2)
    deriv2HHH.fac2 <- Y / meanTotal - psiYpsiMu
    ## psi-related derivatives
    dThetadPsi.fac <- psi * (psiYpsiMu2 - 1/psiPlusMu)
    dThetadPsi <- function(dTheta){
        dThetadPsi.fac * dTheta
    }
    dPsiMat <- psi * (digamma(psiPlusY) - digamma(psi) + log(psi) + 1
                      - log(psiPlusMu) - psiYpsiMu)  # as in penScore()
    dPsidPsiMat <- psi^2 * (
        trigamma(psiPlusY) - trigamma(psi) + 1/psi - 1/psiPlusMu -
            (meanTotal-Y)/psiPlusMu2) + dPsiMat
  } else { # poisson
    deriv2HHH.fac1 <- -Y / (meanTotal^2)
    deriv2HHH.fac2 <- Y / meanTotal - 1
  }
  deriv2HHH <- function(dTheta_l, dTheta_k, dTheta_lk){
      dTheta_l * dTheta_k * deriv2HHH.fac1 + dTheta_lk * deriv2HHH.fac2
  }

  ## go through groups of parameters and compute the hessian of each component
  computeFisher <- function(mean.comp){
    # initialize hessian
    hessian.FE.FE <- matrix(0,dimFE,dimFE)
    hessian.FE.RE <- matrix(0,dimFE,dimRE)
    hessian.RE.RE <- matrix(0,dimRE,dimRE)
    
    hessian.FE.Psi <- matrix(0,dimFE,dimPsi)
    hessian.Psi.RE <- matrix(0,dimPsi,dimPsi+dimRE) # CAVE: contains PsiPsi and PsiRE

    hessian.FE.d <- matrix(0,dimFE,dimd)
    hessian.d.d <- matrix(0,dimd,dimd)
    hessian.d.Psi <- matrix(0,dimd,dimPsi)
    hessian.d.RE <- matrix(0,dimd,dimRE)
    
    ## derivatives wrt neighbourhood weight parameters d
    if (dimd > 0L) {
        phi.doff <- function (dneoff) {
            mu$ne.exppred * dneoff[subset,,drop=FALSE]
        }
        ## for type %in% c("gradient", "hessian"), model$offset[[2L]] always
        ## returns a list of matrices. It has length(pars$d) elements for the
        ## gradient and length(pars$d)*(length(pars$d)+1)/2 for the hessian.
        dneOffset <- model$offset[[2L]](pars$d, type="gradient")
        dmudd <- lapply(dneOffset, phi.doff)
        d2neOffset <- model$offset[[2L]](pars$d, type="hessian")
        d2mudddd <- lapply(d2neOffset, phi.doff)
        ## d l(theta,x) /dd dd (fill only upper triangle, BY ROW)
        ij <- 0L
        for (i in seq_len(dimd)) {
            for (j in i:dimd) {
                ij <- ij + 1L  #= dimd*(i-1) + j - (i-1)*i/2  # for j >= i
                ## d2mudddd contains upper triangle by row (=lowertri by column)
                d2ij <- deriv2HHH(dmudd[[i]], dmudd[[j]], d2mudddd[[ij]])
                hessian.d.d[i,j] <- sum(d2ij, na.rm=TRUE)
            }
        }
    }

    if (dimPsi > 0L) {
        ## d l(theta,x) /dpsi dpsi
        dPsidPsi <- .colSumsGrouped(dPsidPsiMat, indexPsi)
        hessian.Psi.RE[,seq_len(dimPsi)] <- if (dimPsi == 1L) {
            dPsidPsi
        } else {
            diag(dPsidPsi)
        }
        ## d l(theta) / dd dpsi
        for (i in seq_len(dimd)) {      # will not be run if dimd==0
            ## dPsi.i <- colSums(dThetadPsi(dmudd[[i]]),na.rm=TRUE)
            ## hessian.d.Psi[i,] <- if(dimPsi==1L) sum(dPsi.i) else dPsi.i[order(indexPsi)]
            hessian.d.Psi[i,] <- .colSumsGrouped(dThetadPsi(dmudd[[i]]), indexPsi)
        }
    }
    
    ##
    i.fixed <- function(){
      if(random.j){
        Z.j <- term["Z.intercept",j][[1]]  
        "%mj%" <- get(term["mult",j][[1]])
        hessian.FE.RE[idxFE==i,idxRE==j] <<- colSums(didj %mj% Z.j)
        ##<- didj must not contain NA's (all NA's set to 0)
        dIJ <- sum(didj,na.rm=TRUE)     # fixed on 24/09/2012
      } else if(unitSpecific.j){
        dIJ <- colSums(didj,na.rm=TRUE)[ which.j ]
      } else {
        dIJ <- sum(didj,na.rm=TRUE)
      }
      hessian.FE.FE[idxFE==i,idxFE==j] <<- dIJ  
    }
    ##
    i.unit <- function(){
      if(random.j){
        Z.j <- term["Z.intercept",j][[1]] 
        "%mj%" <- get(term["mult",j][[1]]) 
        dIJ <- colSums(didj %mj% Z.j)   # didj must not contain NA's (all NA's set to 0)
        hessian.FE.RE[idxFE==i,idxRE==j] <<- diag(dIJ)[ which.i, ] # FIXME: does not work if type="car"
        dIJ <- dIJ[ which.i ]           # added which.i subsetting in r432
      } else if(unitSpecific.j){
        dIJ <- diag(colSums(didj))[ which.i, which.j ] 
      } else {
        dIJ <- colSums(didj)[ which.i ]
      }
      hessian.FE.FE[idxFE==i,idxFE==j] <<- dIJ  
    }
    ##
    i.random <- function(){
      if(random.j){
        Z.j <- term["Z.intercept",j][[1]]  
        "%mj%" <- get(term["mult",j][[1]])
        hessian.FE.RE[idxFE==i,idxRE==j] <<- colSums(didj %mj% Z.j)
        hessian.FE.RE[idxFE==j,idxRE==i] <<- colSums(didj %m% Z.i)

        if(length(Z.j)==1 & length(Z.i)==1){
          Z <- Z.i*Z.j
          hessian.RE.RE[which(idxRE==i),idxRE==j] <<- diag(colSums( didj %m% Z))
        } else if(length(Z.j)==1 & length(Z.i)>1){         #*
          Z.j <- diag(nrow=model$nUnits)
          for(k in seq_len(ncol(Z.j))){
            Z <- Z.i*Z.j[,k]
            hessian.RE.RE[idxRE==i,which(idxRE==j)[k]] <<- colSums( didj %m% Z)
          }  
        } else if(length(Z.j)>1 & length(Z.i)==1){         #* 
          Z.i <- diag(nrow=model$nUnits)
          for(k in seq_len(ncol(Z.i))){
            Z <- Z.i[,k]*Z.j
            hessian.RE.RE[which(idxRE==i)[k],idxRE==j] <<- colSums( didj %mj% Z)
          }  
        } else {         
          for(k in seq_len(ncol(Z.j))){
            Z <- Z.i*Z.j[,k]
            hessian.RE.RE[which(idxRE==i)[k],idxRE==j] <<- colSums( didj %m% Z)
          }  
        }        
        dIJ <- sum(didj)
      } else if(unitSpecific.j){
        dIJ <- colSums(didj %m% Z.i)
        hessian.FE.RE[idxFE==j,idxRE==i] <<- diag(dIJ)[ which.j, ]
        dIJ <- dIJ[ which.j ]
      } else {
        hessian.FE.RE[idxFE==j,idxRE==i] <<- colSums(didj %m% Z.i)
        dIJ <- sum(didj)
      }
      hessian.FE.FE[idxFE==i,idxFE==j] <<- dIJ    
    }
  ##----------------------------------------------           

    for(i in seq_len(nGroups)){ #go through rows of hessian
      # parameter group belongs to which components
      comp.i <- term["offsetComp",i][[1]]      
      # get covariate value
      Xit <- term["terms",i][[1]] # eiter 1 or a matrix with values
      if(is.matrix(Xit)){
        Xit <- Xit[subset,,drop=FALSE]
      }
      m.Xit <- mean.comp[[comp.i]] * Xit
      
      random.i <- term["random",i][[1]]
      unitSpecific.i <- term["unitSpecific",i][[1]]

      ## fill psi-related entries and select fillHess function
      if (random.i) {
        Z.i <- term["Z.intercept",i][[1]]   # Z.i and %m% (of i) determined here
        "%m%" <- get(term["mult",i][[1]])   # will also be used in j's for loop
        fillHess <- i.random
        if (dimPsi > 0L) {
          dThetadPsiMat <- dThetadPsi(m.Xit)
          hessian.FE.Psi[idxFE==i,] <- .colSumsGrouped(dThetadPsiMat, indexPsi)
          dThetadPsi.i <- .colSums(dThetadPsiMat %m% Z.i, length(subset), term["dim.re",i][[1]], na.rm=TRUE)
          if (dimPsi==1L) {
              hessian.Psi.RE[,dimPsi + which(idxRE==i)] <- dThetadPsi.i
          } else {
              hessian.Psi.RE[cbind(indexPsi,dimPsi + which(idxRE==i))] <- dThetadPsi.i
              ## FIXME: does not work with type="car"
          }
        }
      } else if (unitSpecific.i) {
        which.i <- term["which",i][[1]]
        fillHess <- i.unit
        if (dimPsi > 0L) {
          dThetadPsi.i <- .colSums(dThetadPsi(m.Xit), length(subset), model$nUnits, na.rm=TRUE)
          if (dimPsi==1L) {
              hessian.FE.Psi[idxFE==i,] <- dThetadPsi.i[which.i]
          } else {
              hessian.FE.Psi[cbind(which(idxFE==i),indexPsi[which.i])] <-
                  dThetadPsi.i[which.i]
          }
        }
      } else {
        fillHess <- i.fixed
        if (dimPsi > 0L) {
          ## dPsi <- colSums(dThetadPsi(m.Xit),na.rm=TRUE)
          ## hessian.FE.Psi[idxFE==i,] <- if (dimPsi==1L) sum(dPsi) else dPsi[order(indexPsi)]
          hessian.FE.Psi[idxFE==i,] <- .colSumsGrouped(dThetadPsi(m.Xit), indexPsi)
        }
      }

      ## fill pars$d-related entries
      for (j in seq_len(dimd)) {      # will not be run if dimd==0
          didd <- deriv2HHH(dTheta_l = m.Xit, dTheta_k = dmudd[[j]],
                            dTheta_lk = if (comp.i == 2) dmudd[[j]] * Xit else 0)
          didd[isNA] <- 0
          hessian.FE.d[idxFE==i,j] <- if (unitSpecific.i) {
              colSums(didd,na.rm=TRUE)[which.i]
          } else sum(didd)
          if (random.i) hessian.d.RE[j,idxRE==i] <- colSums(didd %m% Z.i)
      }
      
      ## fill other (non-psi, non-d) entries (only upper triangle, j >= i!)
      for(j in i:nGroups){
        comp.j <- term["offsetComp",j][[1]]
        
        Xjt <- term["terms",j][[1]] # eiter 1 or a matrix with values
        if(is.matrix(Xjt)){
          Xjt <- Xjt[subset,,drop=FALSE]
        }
        # if param i and j do not belong to the same component, d(i)d(j)=0
        m.Xit.Xjt <- if (comp.i != comp.j) 0 else m.Xit * Xjt
        
        didj <- deriv2HHH(dTheta_l = m.Xit, dTheta_k = mean.comp[[comp.j]]*Xjt,
                          dTheta_lk = m.Xit.Xjt)
        didj[isNA]<-0

        random.j <- term["random",j][[1]]      
        unitSpecific.j <- term["unitSpecific",j][[1]]
        which.j <- term["which",j][[1]]
        
        fillHess()
      }
      
    }
    
    
    ######################################################### 
    ## fill lower triangle of hessians and combine them
    ########################################################
    hessian <- rbind(cbind(hessian.FE.FE,hessian.FE.d,hessian.FE.Psi,hessian.FE.RE),
                     cbind(matrix(0,dimd,dimFE),hessian.d.d,hessian.d.Psi,hessian.d.RE),
                     cbind(matrix(0,dimPsi,dimFE+dimd),hessian.Psi.RE),
                     cbind(matrix(0,dimRE,dimFE+dimd+dimPsi),hessian.RE.RE))

    hessian[lower.tri(hessian)] <- 0    # FIXME: should already be the case!
    diagHessian <- diag(hessian)
    fisher <- -(hessian + t(hessian))
    diag(fisher) <- -diagHessian   
  
    return(fisher)
  }

  fisher <- computeFisher(mu[c("epi.own","epi.neighbours","endemic")])

  ## add penalty for random effects
  pen <- matrix(0, length(theta), length(theta))
  Fpen <- if(dimRE > 0){
    thetaIdxRE <- seq.int(to=length(theta), length.out=dimRE)
    pen[thetaIdxRE,thetaIdxRE] <- Sigma.inv
    fisher + pen
  } else fisher

  ## Done
  if(attributes){
    attr(Fpen, "fisher") <- fisher
    attr(Fpen, "pen") <- pen
  }
  Fpen
}


#################################################

sqrtOf1pr2 <- function(r){ sqrt(1+r^2) }

getSigmai <- function(sd,    # vector of length dim with log-stdev's
                      correlation, # vector of length dim with correlation
                                   # parameters, 0-length if uncorrelated
                      dim
                      ){
  if(dim==0) return(NULL)
                      
  Sigma.i <- if (length(correlation) == 0L) diag(exp(2*sd), dim) else {
      D <- diag(exp(sd), dim)
      L <- diag(nrow=dim)
      L[2,1:2] <- c(correlation[1],1)/sqrtOf1pr2(correlation[1])
      if (dim==3) {
          L[3,] <- c(correlation[2:3],1)/sqrtOf1pr2(correlation[2])
          L[3,2:3] <- L[3,2:3]/sqrtOf1pr2(correlation[3])
      }
      D %*% tcrossprod(L) %*% D  # ~75% quicker than D %*% L %*% t(L) %*% D
  }
  return(Sigma.i)
}

getSigmaiInv <- function(sd,    # vector of length dim with log-stdev's
                         correlation, # vector of length dim with correlation
                                      # parameters, 0-length if uncorrelated
                         dim
                         ){
    
  if(dim==0) return(NULL) 
  
  Sigma.i.inv <- if (length(correlation) == 0L) diag(exp(-2*sd), dim) else {
      r <- correlation
      Dinv <- diag(exp(-sd), dim)
      L <- diag(nrow=dim)
      L[2,1:2] <- c(-r[1],sqrtOf1pr2(r[1]))
      if(dim==3){
          L[3,1] <- r[1]*r[3]-r[2]*sqrtOf1pr2(r[3])
          L[3,2] <- -L[2,2]*r[3]
          L[3,3] <- sqrtOf1pr2(r[2])*sqrtOf1pr2(r[3])
      }
      Dinv %*% crossprod(L) %*% Dinv  # ~75% quicker than Dinv %*% t(L) %*% L %*% Dinv
  }
  
  return(Sigma.i.inv)
}

#* allow blockdiagonal matrix blockdiag(A,B), with A=kronecker product, B=diagonal matrix?
getSigmaInv <- function(sd, correlation, dimSigma, dimBlocks, SigmaInvi=NULL){
  if(is.null(SigmaInvi)){
    SigmaInvi <- getSigmaiInv(sd,correlation,dimSigma)
  }
  if(length(unique(dimBlocks))==1){  # kronecker product formulation possible
    kronecker(SigmaInvi,diag(nrow=dimBlocks[1]))
    # the result is a symmetric matrix if SigmaInvi is symmetric
  } else { # kronecker product not possible -> correlation=0 
    diag(rep.int(diag(SigmaInvi),dimBlocks))
  }
}

getSigma <- function(sd, correlation, dimSigma, dimBlocks, Sigmai=NULL){
  if(is.null(Sigmai)){
    Sigmai <- getSigmai(sd,correlation,dimSigma)
  }
  if(length(unique(dimBlocks))==1){  # kronecker product formulation possible
    kronecker(Sigmai,diag(nrow=dimBlocks[1]))
    # the result is a symmetric matrix if Sigmai is symmetric
  } else { # kronecker product not possible -> correlation=0 
    diag(rep.int(diag(Sigmai),dimBlocks))
  }
}


## Approximate marginal likelihood for variance components
## Parameter and model unpacking at the beginning (up to the ###...-line) is
## identical in marScore() and marFisher()
marLogLik <- function(sd.corr, theta, model, fisher.unpen=NULL, verbose=FALSE){

  dimVar <- model$nVar
  dimCorr <- model$nCorr
  dimSigma <- model$nSigma
  
  if(dimSigma == 0){
    return(-Inf)
  }
  
  if(any(is.na(sd.corr))){
    # in order to avoid nlminb from running into an infinite loop (cf. bug
    # report #15052), we have to emergency stop() in this case.
    # As of R 2.15.2, nlminb() throws an error if it receives NA from
    # any of the supplied functions.
    stop("NAs in variance parameters.", ADVICEONERROR)
  }

  sd   <- head(sd.corr,dimVar)
  corr <- tail(sd.corr,dimCorr)

  pars <- splitParams(theta,model)  
  randomEffects <- pars$random
  dimRE <- model$nRE
  
  dimBlocks <- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, dimVar, dimBlocks)

  # if not given, calculate unpenalized part of fisher info 
  if(is.null(fisher.unpen)){
    fisher.unpen <- attr(penFisher(theta, sd.corr, model,attributes=TRUE), "fisher")
  }
  
  # add penalty to fisher
  fisher <- fisher.unpen
  thetaIdxRE <- seq.int(to=length(theta), length.out=dimRE)
  fisher[thetaIdxRE,thetaIdxRE] <- fisher[thetaIdxRE,thetaIdxRE] + Sigma.inv

  ############################################################
  
  # penalized part of likelihood
  # compute -0.5*log(|Sigma|) - 0.5*RE' %*% Sigma.inv %*% RE
  # where -0.5*log(|Sigma|) = -dim(RE_i)*[Sum(sd_i) -0.5*log(1+corr_i^2)]
  ##lpen <- -0.5*(t(randomEffects)%*%Sigma.inv%*%randomEffects)
  ## the following implementation takes ~85% less computing time !
  lpen <- -0.5 * crossprod(randomEffects, Sigma.inv) %*% randomEffects
  loglik.pen <- sum(-dimBlocks*sd) + c(lpen)
  if(dimCorr >0){
    loglik.pen <- loglik.pen + 0.5*dimBlocks[1]*sum(log(1+corr^2))
  }
  
  ## approximate marginal likelihood
  logdetfisher <- determinant(fisher,logarithm=TRUE)$modulus
  lmarg <- loglik.pen -0.5*c(logdetfisher)
      
  return(lmarg)
}


marScore <- function(sd.corr, theta, model, fisher.unpen=NULL, verbose=FALSE){

  dimVar <- model$nVar
  dimCorr <- model$nCorr
  dimSigma <- model$nSigma
  
  if(dimSigma == 0){
    return(numeric(0L))
  }

  if(any(is.na(sd.corr))) stop("NAs in variance parameters.", ADVICEONERROR)
  
  sd   <- head(sd.corr,dimVar)
  corr <- tail(sd.corr,dimCorr)

  pars <- splitParams(theta,model)  
  randomEffects <- pars$random
  dimRE <- model$nRE
  
  dimBlocks <- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, dimVar, dimBlocks)

  # if not given, calculate unpenalized part of fisher info 
  if(is.null(fisher.unpen)){
    fisher.unpen <- attr(penFisher(theta, sd.corr, model,attributes=TRUE), "fisher")
  }
  
  # add penalty to fisher
  fisher <- fisher.unpen
  thetaIdxRE <- seq.int(to=length(theta), length.out=dimRE)
  fisher[thetaIdxRE,thetaIdxRE] <- fisher[thetaIdxRE,thetaIdxRE] + Sigma.inv

  # inverse of penalized fisher info
  F.inv <- try(solve(fisher),silent=TRUE)
  if(inherits(F.inv,"try-error")){
    if(verbose) cat("  WARNING (in marScore): penalized Fisher is singular!\n")
    #return(rep.int(0,dimSigma))
    ## continuing with the generalized inverse often works, otherwise we would
    ## have to stop() here, because nlminb() cannot deal with NA's
    F.inv <- ginv(fisher)
  }
  F.inv.RE <- F.inv[thetaIdxRE,thetaIdxRE]

  ############################################################
  
  ## compute marginal score and fisher for each variance component
  # initialize score and fisher info
  marg.score <- rep.int(NA_real_,dimSigma)
  
  ## specify functions for derivatives
  deriv1 <- switch(dimVar, dSigma1, dSigma2, dSigma3)
  
  d1Sigma <- deriv1(sd, corr)
  Sigmai.inv <- getSigmaiInv(sd, corr, dimVar)
  
 
  # derivation of log determinant
  # -.5*tr(Sigma^-1 %*% dSigma/ds) = -R (for sd.i) 
  #                                =  R*corr.i/(corr.i^2+1) (for corr.i)
  d1logDet <- c(-dimBlocks,dimBlocks[1]*corr/(corr^2+1))
  
  # go through all variance parameters
  for(i in seq_len(dimSigma)){
    dSi <- -Sigmai.inv %*% d1Sigma[,,i] %*% Sigmai.inv # CAVE: sign
    dS.i <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=dSi)
    #dlpen.i <- -0.5* t(randomEffects) %*% dS.i %*% randomEffects
    # ~85% faster implementation using crossprod() avoiding "slow" t():
    dlpen.i <- -0.5 * crossprod(randomEffects, dS.i) %*% randomEffects
    #tr.d1logDetF <- sum(diag(F.inv.RE %*% dS.i))
    tr.d1logDetF <- sum(F.inv.RE * dS.i)   # since dS.i is symmetric
    #<- needs 1/100 (!) of the computation time of sum(diag(F.inv.RE %*% dS.i))
    marg.score[i] <- d1logDet[i] + c(dlpen.i) - 0.5 * tr.d1logDetF
  }
  
  return(marg.score)
}


marFisher <- function(sd.corr, theta, model, fisher.unpen=NULL, verbose=FALSE){
  
  dimVar <- model$nVar
  dimCorr <- model$nCorr
  dimSigma <- model$nSigma
  
  if(dimSigma == 0){
    return(matrix(numeric(0L),0L,0L))
  }

  if(any(is.na(sd.corr))) stop("NAs in variance parameters.", ADVICEONERROR)
  
  sd   <- head(sd.corr,dimVar)
  corr <- tail(sd.corr,dimCorr)

  pars <- splitParams(theta,model)  
  randomEffects <- pars$random
  dimRE <- model$nRE
  
  dimBlocks <- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, dimVar, dimBlocks)

  # if not given, calculate unpenalized part of fisher info 
  if(is.null(fisher.unpen)){
    fisher.unpen <- attr(penFisher(theta, sd.corr, model,attributes=TRUE), "fisher")
  }
  
  # add penalty to fisher
  fisher <- fisher.unpen
  thetaIdxRE <- seq.int(to=length(theta), length.out=dimRE)
  fisher[thetaIdxRE,thetaIdxRE] <- fisher[thetaIdxRE,thetaIdxRE] + Sigma.inv

  # inverse of penalized fisher info
  F.inv <- try(solve(fisher),silent=TRUE)
  if(inherits(F.inv,"try-error")){
    if(verbose) cat("  WARNING (in marFisher): penalized Fisher is singular!\n")
    #return(matrix(Inf,dimSigma,dimSigma))
    ## continuing with the generalized inverse often works, otherwise we would
    ## have to stop() here, because nlminb() cannot deal with NA's
    F.inv <- ginv(fisher)
  }
  F.inv.RE <- F.inv[thetaIdxRE,thetaIdxRE]
  ## declare F.inv.RE as a symmetric matrix?
  ##F.inv.RE <- new("dsyMatrix", Dim = dim(F.inv.RE), x = c(F.inv.RE))
  ## -> no, F.inv.RE %*% dS.i becomes actually slower (dS.i is a "sparseMatrix")
  
  ############################################################
  
  marg.hesse <- matrix(NA_real_,dimSigma,dimSigma)
  
  ## specify functions for derivatives
  deriv1 <- switch(dimVar,dSigma1, dSigma2, dSigma3)
  deriv2 <- switch(dimVar,d2Sigma1, d2Sigma2, d2Sigma3)
  
  d1Sigma <- deriv1(sd, corr)
  d2Sigma <- deriv2(sd, corr, d1Sigma)
  Sigmai.inv <- getSigmaiInv(sd, corr, dimVar)
  
 
  # 2nd derivatives of log determinant
  d2logDet <- diag(c(rep.int(0,dimVar),-dimBlocks[1]*(corr^2-1)/(corr^2+1)^2),dimSigma)

  # function to convert dS.i and dS.j matrices to sparse matrix objects
  dS2sparse <- if (dimCorr > 0) function (x) {
      forceSymmetric(as(x, "sparseMatrix")) # dS.i & dS.j are symmetric
  } else function (x) { #as(x, "diagonalMatrix")
      new("ddiMatrix", Dim = dim(x), diag = "N", x = diag(x))
  }
  
  # go through all variance parameters
  for(i in seq_len(dimSigma)){
    # compute first derivative of the penalized Fisher info (-> of Sigma^-1) 
    # with respect to the i-th element of Sigma (= kronecker prod. of Sigmai and identity matrix)
    # Harville Ch15, Eq. 8.15: (d/d i)S^-1 = - S^-1 * (d/d i) S * S^-1
    SigmaiInv.d1i <- Sigmai.inv %*% d1Sigma[,,i]
    dSi <- -SigmaiInv.d1i %*% Sigmai.inv
    dS.i <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=dSi)
    dS.i <- dS2sparse(dS.i)
    
    # compute second derivatives
    for(j in i:dimSigma){
      # compute (d/d j) S^-1
      SigmaiInv.d1j <- Sigmai.inv %*% d1Sigma[,,j]
      dSj <- -SigmaiInv.d1j %*% Sigmai.inv
      dS.j <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=dSj)
      dS.j <- dS2sparse(dS.j)
      # compute (d/di dj) S^-1
      #dS.ij <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,
      #                  Sigmai=d2Sigma[[i]][,,j])
      
      # compute second derivatives of Sigma^-1 (Harville Ch15, Eq 9.2)
      d2S <- (- Sigmai.inv %*% d2Sigma[[i]][,,j] + 
               SigmaiInv.d1i %*% SigmaiInv.d1j +
               SigmaiInv.d1j %*% SigmaiInv.d1i) %*% Sigmai.inv 
      
      dSij <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=d2S)

      #d2lpen.i <- -0.5* t(randomEffects) %*% dSij %*% randomEffects
      # ~85% faster implementation using crossprod() avoiding "slow" t():
      d2lpen.i <- -0.5 * crossprod(randomEffects, dSij) %*% randomEffects

      # compute second derivative of log-determinant of penFisher
      mpart1 <- dS.j %*% F.inv.RE  # 3 times as fast as the other way round
      mpart2 <- dS.i %*% F.inv.RE
      mpart <- mpart1 %*% mpart2
      ## speed-ups: - tr(F.inv.RE %*% dSij) simply equals sum(F.inv.RE * dSij)
      ##            - accelerate matrix product by sparse matrices dS.i and dS.j
      ##            - use cyclic permutation of trace:
      ##              tr(F.inv.RE %*% dS.j %*% F.inv.RE %*% dS.i) =
      ##              tr(dS.j %*% F.inv.RE %*% dS.i %*% F.inv.RE)
      tr.d2logDetF <- -sum(Matrix::diag(mpart)) + sum(F.inv.RE * dSij)
      
      marg.hesse[i,j] <- marg.hesse[j,i] <-
          d2logDet[i,j] + d2lpen.i - 0.5 * tr.d2logDetF
    }
    
  }  
  
  marg.Fisher <- as.matrix(-marg.hesse)
   
  return(marg.Fisher)
}


## first and second derivatives of the covariance matrix
dSigma1 <- function(sd,corr){
  derivs <- array(2*exp(2*sd), c(1,1,1)) 
  return(derivs)
}

#d1: result of dSigma1
d2Sigma1 <- function(sd,corr,d1){
  return(list(dsd1=2*d1))
}

dSigma2 <- function(sd,corr){
  derivs <- array(0,c(2,2,3))
  
  dSigma <- diag(2*exp(2*sd))
  if(length(corr)>0){
    dSigma[1,2] <- dSigma[2,1] <- exp(sum(sd[1:2]))*corr[1]/sqrtOf1pr2(corr[1])
    
    # derivative of corr_1
    derivs[2,1,3] <- derivs[1,2,3] <- exp(sum(sd[1:2]))/(sqrtOf1pr2(corr[1])^3)
  }
  
  derivs[,,1:2] <- dSigma
  # derivative of sd_1
  derivs[2,2,1] <- 0
  # derivative of sd_2
  derivs[1,1,2] <- 0

  return(derivs)
}

d2Sigma2 <- function(sd,corr, d1){
  
  derivs <- array(0,c(2,2,3))
  result <- list(dsd1=d1, dsd2=derivs, dcorr1=derivs)
  result$dsd1[1,1,1] <- 2*d1[1,1,1]
  result$dsd1[2,2,2] <- 0
  
  result$dsd2[,,2:3]<- d1[,,2:3]
  result$dsd2[2,2,2] <- 2*d1[2,2,2]
    
    
  if(length(corr)>0){
    result$dcorr1[2,1,3] <- result$dcorr1[1,2,3] <-
        -(3*corr[1]*exp(sum(sd[1:2])))/(sqrtOf1pr2(corr[1])^5)
  }
  
  return(result)
}


dSigma3 <- function(sd,corr){
  derivs <- array(0,c(3,3,6))
  
  dSigma <- diag(2*exp(2*sd)) #
  if(length(corr)>0){
    dSigma[1,2] <- dSigma[2,1] <- exp(sum(sd[1:2]))*corr[1]/sqrtOf1pr2(corr[1]) #
    dSigma[1,3] <- dSigma[3,1] <- exp(sum(sd[c(1,3)]))*corr[2]/sqrtOf1pr2(corr[2]) #
    dSigma[2,3] <- dSigma[3,2] <- exp(sum(sd[c(2,3)]))*(corr[1]*corr[2]*sqrtOf1pr2(corr[3])+corr[3])/prod(sqrtOf1pr2(corr[1:3]))#
    
    # derivative of corr_1
    derivs[2,1,4] <- derivs[1,2,4] <- exp(sum(sd[1:2]))/(sqrtOf1pr2(corr[1])^3)
    derivs[3,2,4] <- derivs[2,3,4] <-(exp(sum(sd[2:3]))*(corr[2]*sqrtOf1pr2(corr[3])-prod(corr[c(1,3)])))/ (prod(sqrtOf1pr2(corr[2:3]))*(sqrtOf1pr2(corr[1])^3))#
    
    # derivative of corr_2
    derivs[3,1,5] <- derivs[1,3,5] <- exp(sum(sd[c(3,1)]))/(sqrtOf1pr2(corr[2])^3)#
    derivs[3,2,5] <- derivs[2,3,5] <- (exp(sum(sd[2:3]))*(corr[1]*sqrtOf1pr2(corr[3])-prod(corr[c(2,3)])))/ (prod(sqrtOf1pr2(corr[c(1,3)]))*(sqrtOf1pr2(corr[2])^3)) #
    
    # derivative of corr_3
    derivs[3,2,6] <- derivs[2,3,6] <- exp(sum(sd[2:3]))/ (prod(sqrtOf1pr2(corr[c(1,2)]))*(sqrtOf1pr2(corr[3])^3))
  }
  
  derivs[,,1:3] <- dSigma
  # derivative of sd_1
  derivs[2:3,2:3,1] <- 0
  # derivative of sd_2
  derivs[1,c(1,3),2] <- derivs[3,c(1,3),2] <- 0
  # derivative of sd_3
  derivs[1:2,1:2,3] <- 0
  
  return(derivs)

}

d2Sigma3 <- function(sd,corr, d1)
{  
  derivs <- array(0,c(3,3,6))
  result <- list(dsd1=d1, dsd2=derivs, dsd3=derivs, dcorr1=derivs, dcorr2=derivs, dcorr3=derivs)
  
  result$dsd1[1,1,1] <- 2*d1[1,1,1]
  result$dsd1[2,2:3,2] <- result$dsd1[3,2,2] <- 0
  result$dsd1[2:3,2:3,3] <-  0 #
  
  result$dsd2[,,2]<- d1[,,2]
  result$dsd2[2,2,2] <- 2*d1[2,2,2]
  result$dsd2[3,2,3] <- result$dsd2[2,3,3] <- d1[3,2,3]#
    
  result$dsd3[,,3]<- d1[,,3]
  result$dsd3[3,3,3] <- 2*d1[3,3,3]#
  
  if (length(corr)>0)
  {
    result$dsd1[2:3,2:3,4] <-  0
    result$dsd1[2:3,2:3,5] <-  0
    result$dsd1[,,6] <-  0
    
    result$dsd2[,,c(4,6)] <- d1[,,c(4,6)]
    result$dsd2[3,2,5] <- result$dsd2[2,3,5] <- d1[3,2,5]
    
    result$dsd3[3,2,4] <- result$dsd3[2,3,4] <- d1[3,2,4]
    result$dsd3[,,c(5,6)] <- d1[,,c(5,6)]

    # derivative of corr_1
    result$dcorr1[2,1,4] <- result$dcorr1[1,2,4] <- -(exp(sum(sd[1:2]))*3*corr[1])/(sqrtOf1pr2(corr[1])^5) #
    result$dcorr1[3,2,4] <- result$dcorr1[2,3,4] <- -(exp(sum(sd[2:3]))*(corr[1]*(3*corr[2]*sqrtOf1pr2(corr[3])-2*prod(corr[c(1,3)])) + corr[3]) )/ (prod(sqrtOf1pr2(corr[2:3]))*(sqrtOf1pr2(corr[1])^5)) #
    
    result$dcorr1[3,2,5] <- result$dcorr1[2,3,5] <- (exp(sum(sd[2:3]))*(sqrtOf1pr2(corr[3])+prod(corr[1:3])))/ (prod(sqrtOf1pr2(corr[c(1,2)])^3)*sqrtOf1pr2(corr[3]))

    result$dcorr1[3,2,6] <- result$dcorr1[2,3,6] <- -(exp(sum(sd[2:3]))*corr[1])/ (prod(sqrtOf1pr2(corr[c(1,3)])^3)*sqrtOf1pr2(corr[2]))
      
    # derivative of corr_2
    result$dcorr2[3,1,5] <- result$dcorr2[1,3,5] <- -(exp(sum(sd[c(3,1)]))*3*corr[2])/(sqrtOf1pr2(corr[2])^5)
    result$dcorr2[3,2,5] <- result$dcorr2[2,3,5] <- -(exp(sum(sd[2:3]))*(corr[2]*(3*corr[1]*sqrtOf1pr2(corr[3])-2*prod(corr[c(2,3)])) + corr[3]) )/ (prod(sqrtOf1pr2(corr[c(1,3)]))*(sqrtOf1pr2(corr[2])^5))
    
    result$dcorr2[3,2,6] <- result$dcorr2[2,3,6] <-
        -exp(sum(sd[2:3]))*corr[2] / # SM @ 14/05/13: formula fixed, marFisher()
                                     # and hhh4()$Sigma.cov[5,6] are now correct
            (prod(sqrtOf1pr2(corr[c(2,3)])^3)*sqrtOf1pr2(corr[1]))
    
    # derivative of corr_3
    result$dcorr3[3,2,6] <- result$dcorr3[2,3,6] <- -(exp(sum(sd[2:3]))*3*corr[3])/ (prod(sqrtOf1pr2(corr[c(1,2)]))*sqrtOf1pr2(corr[3])^5)
  }

  return(result)
}


### Various optimizers

updateParams_nlminb <- function (start, ll, sc, fi, ..., control)
{
    lower <- control[["lower"]]; control$lower <- NULL
    upper <- control[["upper"]]; control$upper <- NULL
    scale <- control[["scale"]]; control$scale <- NULL
    negll <- function (x, ...) -ll(x, ...)
    negsc <- function (x, ...) -sc(x, ...)
    ## run the optimization
    res <- nlminb(start, negll, gradient=negsc, hessian=fi, ...,
                  scale=scale, control=control, lower=lower, upper=upper)
    if (any(is.finite(c(lower, upper)))) checkParBounds(res$par, lower, upper)
    ## Done
    list(par=res$par, ll=-res$objective,
         rel.tol=getRelDiff(res$par, start),
         convergence=res$convergence, message=res$message)
}

updateParams_nr <- function (start, ll, sc, fi, ..., control)
{
    ## objective function
    llscfi <- function (x, ...) {
        loglik <- ll(x, ...)
        attr(loglik, "score") <- sc(x, ...)
        attr(loglik, "fisher") <- fi(x, ...)
        loglik
    }
    ## run the optimization
    res <- newtonRaphson(start, llscfi, ...,
                         control=control, verbose=control$verbose)
    ## Done
    list(par=res$coefficients, ll=res$loglikelihood,
         rel.tol=getRelDiff(res$coefficients, start),
         convergence=res$convergence, message=res$message)
}

updateParams_nlm <- function (start, ll, sc, fi, ..., control)
{
    ## objective function
    negllscfi <- function (x, ...) {
        negloglik <- -ll(x, ...)
        attr(negloglik, "gradient") <- -sc(x, ...)
        attr(negloglik, "hessian") <- fi(x, ...)
        negloglik
    }
    ## run the optimization
    res <- do.call("nlm", args=c(alist(p=start, f=negllscfi, ...), control))
    ## Done
    list(par=res$estimate, ll=-res$minimum,
         rel.tol=getRelDiff(res$estimate, start),
         convergence=as.numeric(res$code>2), message=res$message)
    ## nlm returns convergence status in $code, 1-2 indicate convergence,
    ## 3-5 indicate non-convergence
}

updateParams_optim <- function (start, ll, sc, fi, ..., control)
{
    ## Note: "fi" is not used in optim
    method <- control[["method"]]; control$method <- NULL
    lower <- control[["lower"]]; control$lower <- NULL
    upper <- control[["upper"]]; control$upper <- NULL
    res <- optim(start, ll, sc, ...,    # Note: control$fnscale is negative
                 method=method, lower=lower, upper=upper, control=control)
    if (any(is.finite(c(lower, upper)))) checkParBounds(res$par, lower, upper)
    ## Done
    list(par=res$par, ll=res$value,
         rel.tol=getRelDiff(res$par, start),
         convergence=res$convergence, message=res$message)
}


## Calculate relative parameter change criterion.
## We use a weaker criterion than the maximum relative parameter change
## max(abs(sd.corr.new/sd.corr - 1))
getRelDiff <- function (final, start)
    max(abs(final - start)) / max(abs(start))

checkParBounds <- function (par, lower, upper)
{
    if (is.null(names(par))) names(par) <- seq_along(par)
    if (any(atl <- par <= lower))
        cat("  WARNING: parameters reached lower bounds:",
            paste(names(par)[atl], par[atl], sep="=", collapse=", "), "\n")
    if (any(atu <- par >= upper))
        cat("  WARNING: parameters reached upper bounds:",
            paste(names(par)[atu], par[atu], sep="=", collapse=", "), "\n")
}


## default control arguments for updates
defaultOptimControl <- function (method = "nlminb", lower = -Inf, upper = Inf,
                                 iter.max = NULL, verbose = 0)
{
    if (is.null(iter.max)) iter.max <- 20 + 280*(method=="Nelder-Mead")
    lowVerbose <- verbose %in% 0:2
    luOptimMethod <- method %in% c("Brent", "L-BFGS-B")
    defaults.nr <- list(scoreTol=1e-5, paramTol=1e-7, F.inc=0.01, stepFrac=0.5,
                        niter=iter.max, verbose=verbose)
    defaults.nlminb <- list(iter.max=iter.max, scale=1, lower=lower, upper=upper,
                            trace=if(lowVerbose) c(0,0,5)[verbose+1] else 1)
    defaults.nlm <- list(iterlim=iter.max, check.analyticals=FALSE,
                         print.level=if(lowVerbose) c(0,0,1)[verbose+1] else 2)
    defaults.optim <- list(maxit=iter.max, fnscale=-1, trace=max(0,verbose-1),
                           lower=if (luOptimMethod) lower else -Inf,
                           upper=if (luOptimMethod) upper else Inf)
    switch(method, "nr" = defaults.nr, "nlm" = defaults.nlm,
           "nlminb" = defaults.nlminb, defaults.optim)
}

setOptimControl <- function (method, control, ...)
{
    defaults <- defaultOptimControl(method, ...)
    cntrl <- modifyList(defaults, control)
    ## ensure fnscale < 0 (optim performs minimization)
    if (!is.null(cntrl$fnscale)) {      # i.e., using optim()
        cntrl$method <- method          # append method to control list
        if (cntrl$fnscale > 0) cntrl$fnscale <- -cntrl$fnscale
    }
    cntrl
}


## fitHHH is the main workhorse where the iterative optimization is performed
fitHHH <- function(theta, sd.corr, model,
                   cntrl.stop=list(tol=1e-5, niter=100),
                   cntrl.regression=list(method="nlminb"),
                   cntrl.variance=list(method="nlminb"),
                   verbose=0, shrinkage=FALSE)
{  
  dimFE.d.O <- model$nFE + model$nd + model$nOverdisp
  dimRE <- model$nRE

  getUpdater <- function (cntrl, start, ...) {
      method <- cntrl$method; cntrl$method <- NULL
      if (length(start) == 1 && method == "Nelder-Mead") {
          method <- "Brent"
          message("Switched optimizer from \"Nelder-Mead\" to \"Brent\"",
                  " (dim(", deparse(substitute(start)), ")=1)")
      }
      list(paste("updateParams",
                 if (method %in% c("nlminb", "nlm", "nr"))
                 method else "optim", sep="_"),
           control = setOptimControl(method, cntrl, ...))
  }

  ## ## artificial lower bound on intercepts of epidemic components
  ## reg.lower <- rep.int(-Inf, length(theta))
  ## reg.lower[grep("^(ar|ne)\\.(1|ri)", model$namesFE)] <- -20
  
  ## set optimizer for regression parameters
  updateRegressionControl <- getUpdater(cntrl.regression, theta,
                                        ## lower=reg.lower,
                                        iter.max=if(dimRE==0) 100,
                                        verbose=verbose+(dimRE==0))
  updateRegression <- function (theta, sd.corr)
      do.call(updateRegressionControl[[1]],
              alist(theta, penLogLik, penScore, penFisher,
                    sd.corr=sd.corr, model=model,
                    control=updateRegressionControl[[2]]))

  ## set optimizer for variance parameters
  updateVarianceControl <- getUpdater(cntrl.variance, sd.corr,
                                      lower=-5, upper=5, verbose=verbose)
  updateVariance <- function (sd.corr, theta, fisher.unpen)
      do.call(updateVarianceControl[[1]],
              alist(sd.corr, marLogLik, marScore, marFisher,
                    theta=theta, model=model,
                    fisher.unpen=fisher.unpen, verbose=verbose>1,
                    control=updateVarianceControl[[2]]))

  ## Let's go
  if (verbose>0) {
      cat(as.character(Sys.time()), ":",
          if (dimRE == 0) "Optimization of regression parameters" else
          "Iterative optimization of regression & variance parameters", "\n")
  }
  
  if (dimRE == 0) { # optimization of regression coefficients only
      parReg <- updateRegression(theta, sd.corr)
      theta <- parReg$par
      if ((convergence <- parReg$convergence) != 0 && !is.null(parReg$message))
          cat("! Non-convergence message from optimizer:", parReg$message, "\n")
  } else { # swing between updateRegression & updateVariance
      convergence <- 99
      i <- 0
      while(convergence != 0 && (i < cntrl.stop$niter)){
          i <- i+1
          if (verbose>0) cat("\n")
          
          ## update regression coefficients
          parReg <- updateRegression(theta, sd.corr)
          theta <- parReg$par
          fisher.unpen <- attr(penFisher(theta, sd.corr, model, attributes=TRUE),
                               "fisher")

          if(verbose>0)
              cat("Update of regression parameters: ",
                  "max|x_0 - x_1| / max|x_0| =", parReg$rel.tol, "\n")
          
          if(parReg$convergence != 0) {
              if (!is.null(parReg$message))
                  cat("! Non-convergence message from optimizer:",
                      parReg$message, "\n")
              cat("Update of regression coefficients in iteration ",
                  i, " unreliable\n")
          }

          if(parReg$convergence > 20 && shrinkage){
              cat("\n\n***************************************\nshrinkage",
                  0.1*theta[abs(theta)>10],"\n")
              theta[abs(theta)>10] <- 0.1*theta[abs(theta)>10]
              diag(fisher.unpen) <- diag(fisher.unpen)+1e-2
          }
          
          ## update variance parameters
          parVar <- updateVariance(sd.corr, theta, fisher.unpen)

          if(verbose>0)
              cat("Update of variance parameters:  max|x_0 - x_1| / max|x_0| =",
                  parVar$rel.tol, "\n")
          
          if(parVar$convergence!=0) {
              if (!is.null(parVar$message)) print(parVar$message)
              cat("Update of variance parameters in iteration ", i, " unreliable\n")
          }

          ## NA values in sd.corr cause a stop() already in marLogLik()
          ## if(any(is.na(parVar$par))){
          ##     updateVarianceControl[[1]] <- "updateParams_optim"
          ##     updateVarianceControl[[2]]$method <-
          ##         if (length(sd.corr) == 1L) "Brent" else "Nelder-Mead"
          ##     cat("  WARNING: at least one updated variance parameter is not a number\n",
          ##         "\t-> NO UPDATE of variance\n",
          ##         "\t-> SWITCHING to robust", dQuote(updateVarianceControl[[2]]$method),
          ##         "for variance updates\n")
          ## } else
          sd.corr <- parVar$par

          ## overall convergence ?
          if( (parReg$rel.tol < cntrl.stop$tol) && (parVar$rel.tol < cntrl.stop$tol)
             && (parReg$convergence==0) && (parVar$convergence==0) ) 
              convergence <- 0

          ## exit loop if no more change in parameters (maybe false convergence)
          if (parReg$rel.tol == 0 && parVar$rel.tol == 0)
              break
      }
  }
  
  if(verbose > 0) {
    cat("\n")
    cat(as.character(Sys.time()), ":", if (convergence==0)
        "Optimization converged" else "Optimization DID NOT CONVERGE", "\n\n")
  }
  
  ll <- penLogLik(theta=theta,sd.corr=sd.corr,model=model)
  fisher <- penFisher(theta=theta,sd.corr=sd.corr,model=model)
  fisher.var <- marFisher(sd.corr=sd.corr, theta=theta, model=model,
                          fisher.unpen=fisher.unpen)
  
  list(fixef=head(theta,dimFE.d.O), ranef=tail(theta,dimRE), sd.corr=sd.corr,
       loglik=ll, fisher=fisher, fisherVar=fisher.var,
       convergence=convergence, dim=c(fixed=dimFE.d.O,random=dimRE))
}



## check analytical score functions and Fisher informations for
## a given model (the result of interpretControl(control, stsObj))
## and given parameters theta (regression par.) and sd.corr (variance par.).
## This is a wrapper around functionality of the numDeriv and maxLik packages.
checkAnalyticals <- function (model,
                              theta = model$initialTheta,
                              sd.corr = model$initialSigma,
                              methods = c("numDeriv","maxLik"))
{
    cat("\nPenalized log-likelihood:\n")
    resCheckPen <- sapply(methods, function(derivMethod) {
        if (requireNamespace(derivMethod)) {
            do.call(paste("checkDerivatives", derivMethod, sep="."),
                    args=alist(penLogLik, penScore, penFisher, theta,
                    sd.corr=sd.corr, model=model))
        }
    }, simplify=FALSE, USE.NAMES=TRUE)
    if (length(resCheckPen) == 1L) resCheckPen <- resCheckPen[[1L]]
    
    resCheckMar <- if (length(sd.corr) == 0L) list() else {
        cat("\nMarginal log-likelihood:\n")
        fisher.unpen <- attr(penFisher(theta, sd.corr, model, attributes=TRUE),
                             "fisher")
        resCheckMar <- sapply(methods, function(derivMethod) {
            if (requireNamespace(derivMethod)) {
                do.call(paste("checkDerivatives", derivMethod, sep="."),
                        args=alist(marLogLik, marScore, marFisher, sd.corr,
                        theta=theta, model=model, fisher.unpen=fisher.unpen))
            }
        }, simplify=FALSE, USE.NAMES=TRUE)
        if (length(resCheckMar) == 1L) resCheckMar[[1L]] else resCheckMar
    }
    
    list(pen = resCheckPen, mar = resCheckMar)
}
