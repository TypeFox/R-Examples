`checkDist` <- function(distname.y, scale.OT = 1) {

  ##===========================================================================
  ## Check the provided distribution name, and find the suitable function
  ## (root) name as well as parameter names.
  ##===========================================================================

  special.y <- TRUE
  
  if (distname.y == "exponential") {
    funname.y <- "exp"
    parnames.y <- "rate"
    parLower.y <- 0
    parUpper.y <- Inf
    names(parUpper.y) <- names(parLower.y) <- parnames.y
    scale.y <- 1/scale.OT
  } else if (distname.y == "weibull") {
    funname.y <- "weibull"
    parnames.y <- c("shape", "scale")
    parLower.y <- c(0, 0)
    parUpper.y <- c(Inf, Inf)
    names(parUpper.y) <- names(parLower.y) <- parnames.y
    scale.y <- c(1, scale.OT)
  } else if (distname.y == "gpd") {
    ## At the time, evd is necessary
    funname.y <- "gpd"
    parnames.y <- c("scale", "shape")
    parLower.y <- c(0, -2.0)
    parUpper.y <- c(Inf, 2.0)
    names(parUpper.y) <- names(parLower.y) <- parnames.y
    scale.y <- c(scale.OT, 1)
  } else if (distname.y == "GPD") {
    ## At the time, evd is necessary
    funname.y <- "GPD"
    parnames.y <- c("scale", "shape")
    parLower.y <- c(0, -2.0)
    parUpper.y <- c(Inf, 2.0)
    names(parUpper.y) <- names(parLower.y) <- parnames.y
    scale.y <- c(scale.OT, 1)
  } else if (distname.y %in% c("log-normal", "lognormal")) {
    distname.y <- "log-normal"
    funname.y <- "lnorm"
    parnames.y <- c("meanlog", "sdlog")
    parLower.y <- c(-Inf, 0)
    parUpper.y <- c(Inf, Inf)
    names(parUpper.y) <- names(parLower.y) <- parnames.y
    scale.y <- c(1, 1)
  } else if (distname.y == "gamma"){
    funname.y <- "gamma"
    parnames.y <- c("shape", "scale")
    parLower.y <- c(0, 0)
    parUpper.y <- c(Inf, Inf)
    names(parUpper.y) <- names(parLower.y) <- parnames.y
    scale.y <- c(1, scale.OT)
  } else if (distname.y == "lomax"){
    funname.y <- "lomax"
    parnames.y <- c("shape", "scale")
    parLower.y <- c(0, 0)
    parUpper.y <- c(Inf, Inf)
    names(parUpper.y) <- names(parLower.y) <- parnames.y
    scale.y <- c(1, scale.OT)
  }  else if (distname.y == "maxlo"){
    funname.y <- "maxlo"
    parnames.y <- c("shape", "scale")
    parLower.y <- c(0, 0)
    parUpper.y <- c(Inf, Inf)
    names(parUpper.y) <- names(parLower.y) <- parnames.y
    scale.y <- c(1, scale.OT)
  } else if (distname.y %in% c("MixExp2", "mixexp2")) {
    distname.y <- "mixexp2"
    funname.y <- "mixexp2"
    parnames.y <- c("prob1", "rate1", "delta")
    parLower.y <- c(0, 0, 0)
    parUpper.y <- c(1, Inf, Inf)
    names(parUpper.y) <- names(parLower.y) <- parnames.y
    scale.y <- c(1, 1/scale.OT, 1/scale.OT)
  } else{
    warning("warning: distribution not in target list. Still EXPERIMENTAL")
    special.y <- FALSE
    funname.y <- distname.y
    parnames.y <- NULL
    parLower.y <- NULL
    parUpper.y <- NULL
    scale.y <- NULL
  }

  list(special.y = special.y,
       distname.y = distname.y,
       funname.y = funname.y,
       parnames.y = parnames.y,
       parLower.y = parLower.y,
       parUpper.y = parUpper.y,
       scale.y = scale.y)

}

##*****************************************************************************

`transFuns` <-
  function(trans.y,
           distname.y) {
    
    ##=========================================================================
    ## prepare transforms if wanted/possible
    ##=========================================================================
    
    if (!is.null(trans.y)) {
      if( !is.character(trans.y) || !(trans.y %in% c("square", "log")) ) 
        stop("trans.y must be NULL or be character in c( \"square\', \"log\")")
      else if (distname.y != "exponential") {
        stop("non-null value for 'trans.y' is only allowed when distname.y == \"exponential\"") 
      } else {
        transFlag <- TRUE
        if (trans.y == "square") {
          transfun <- function(x) x^2 
          invtransfun <- get("sqrt", mode = "function")   
        }
        if (trans.y == "log") {
          transfun <- get("log", mode = "function")   
          invtransfun <- get("exp", mode = "function")   
        }
      }
    } else {
      transFlag <- FALSE
      transfun <- NULL
      invtransfun <- NULL
    }

    list(transFlag = transFlag,
         transfun = transfun,
         invtransfun = invtransfun)
  
}

##*****************************************************************************

`makeFuns` <-
  function(funname.y,
           parnames.y,
           fixed.par.y = list(),
           trace = 1) {

    ##=========================================================================
    ## Analysis of provided values for parameters (FIXED parameters)
    ##
    ## 'fixed.y'    is a flag indicating which parameter is fixed.
    ##              It must have the same length and the same order as
    ##              'parnames.y'
    ## 'p.y'        is the number of parameters for y in estimation
    ##              1 <= p.y <= parnb.y.
    ## 'pf.y'       is the number of fixed parameters.
    ##
    ##=========================================================================
    
    fixed.par.y <- unlist(fixed.par.y)
    
    m <- match(names(fixed.par.y), table = parnames.y)
    
    if ( any(is.na(m)) ) stop("some names not understood in 'fixed.par.y'")
    
    if ( length(fixed.par.y) && any(is.na(fixed.par.y)) )
      stop("fixed.par.y must contain only non-missing values")
    
    fixed.y <- rep(FALSE, length(parnames.y))
    names(fixed.y) <- parnames.y
    fixed.y[m] <- TRUE
    
    p.y <- length(parnames.y[!fixed.y])
    if(p.y == 0) warning("No parameter to estimate")
    
    if ( any(fixed.y)) {
      pf.y <- sum(fixed.y)
      if (trace) {
        cat("o Fixed parameters for y\n")
        print(fixed.y)
      }
    } else pf.y <- 0
    
    if (trace) {
      cat("o Number of param for exceedances to estimate from data\n")
      cat("  p.y = ", p.y, "\n")
    }
    
    ##=========================================================================
    ## Make a list with the distributions functions
    ## 
    ## In this part, formals of functions are changed following the
    ## ideas in "fitdistr" of the MASS package. See therein.
    ##
    ## Note: the code could be rearranged using a loop on function names
    ## 
    ##=========================================================================
    
    ## reorder arguments to densfun and co
    dfun.y <- get(paste("d", funname.y, sep = ""), mode = "function")
    pfun.y <- get(paste("p", funname.y, sep = ""), mode = "function")
    qfun.y <- get(paste("q", funname.y, sep = ""), mode = "function")
    
    fms <- formals(dfun.y)
    args <- names(fms)
    m <- match(parnames.y, args)
    
    if(any(is.na(m)))
      stop("'parnames.y' specifies names which are not arguments to 'densfun'")
    
    ## 'x' is maintened in pole position '1', then come the params in parnames,
    ## then the other if any (e.g. surrogate parameters)
    formals(dfun.y) <- c(fms[c(1L, m)], fms[-c(1L, m)])
    
    ## Caution: in function call, remember to use  log = TRUE
    
    logf.y <- function(parm, x) dfun.y(x, parm, log = TRUE)
    
    ## Same thing for 'pfun'. Note that although the main arg of
    ## distributions functions is usually names "q", we force it 
    ## to be "x" here because f.y and F.y are usaed in the same
    ## manner in the log-likelihood!
    
    fms <- formals(pfun.y)
    args <- names(fms)
    m <- match(parnames.y, args)
    if(any(is.na(m)))
      stop("parnames.y specifies names which are not arguments to 'pfun'")
    
    formals(pfun.y) <- c(fms[c(1L, m)], fms[-c(1L, m)])
    
    F.y <- function(parm, x) pfun.y(x, parm)
    
    ## reorder arguments to densfun
    fms <- formals(qfun.y)
    args <- names(fms)
    m <- match(parnames.y, args)
    if(any(is.na(m)))
      stop("'parnames.y' specifies names which are not arguments to 'qfun'")
    
    formals(qfun.y) <- c(fms[c(1L, m)], fms[-c(1L, m)])
    
    q.y <- function(parm, p) qfun.y(p, parm)
    
    ##=========================================================================
    ## Hack formals and body for the wrapper functions
    ## The case p.y == 0 (no estimation) is not possible at the time
    ##=========================================================================
    
    if (p.y >= 1L) {
      
      ## to remove later!!!
      str <- paste(paste("parm[", 1:p.y, "]", collapse = ", ", sep = ""), ")")
      
      nfn <- parnames.y[!fixed.y]
      
      str <- paste(paste(paste(nfn, paste("parm[\"", nfn, "\"]", sep = ""), sep = " = "),
                         collapse = ", "), sep = "")
      
    } else {
      if (pf.y == 0L)
        stop("no parameter to estimate and no parameter fixed. Check distribution")
    }
    
    ## Add fixed values if necessary
    ## This is done by modifying the body of the functions by calling
    ## the relevant function in it with suitable NAMED args. 
    
    if (pf.y >= 1L) {
      ## modif du 2010-01-25 for fixed par
      if (p.y >= 1L) {
        strf <- paste(paste(paste(names(fixed.par.y), fixed.par.y, sep = " = "),
                            collapse = ", "), sep = "")
        str <- paste(str, strf, sep = ", ")
      } else {
        str <- paste(paste(paste(names(fixed.par.y), fixed.par.y, sep = " = "),
                           collapse = ", "), sep = "")
      }
      
    } 
    
    body(logf.y) <- parse(text = paste("dfun.y(x,", str,", log = TRUE)") )     
    body(F.y) <- parse(text = paste("pfun.y(x,", str, ")"))
    body(q.y) <- parse(text = paste("qfun.y(p,", str, ")"))
    
    ##=========================================================================
    ## a list of functions  to be exported.
    ## Note that the definition of 'dfun.y', 'qfun.y' and 'pfun.y'
    ## is taken from the environment where the funs are defined, i.e.
    ## here (lexical scoping).
    ##
    ## CAUTION
    ##
    ## Do not re-use the functions 'logf.y', 'q.y' ot 'F.y' in any
    ## environment where the definition of 'dfun.y', 'pfun.y' or 'qfun.y'
    ## could be different!
    ##
    ##=========================================================================
    
    funs <- list(fixed.y = fixed.y,
                 pf.y = pf.y,
                 p.y = p.y,
                 dfun.y = dfun.y,
                 pfun.y = pfun.y,
                 qfun.y = qfun.y,
                 logf.y = logf.y,
                 q.y = q.y, 
                 F.y = F.y)
  }
