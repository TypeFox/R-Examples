## =============================================================================
## MCMC auxilliary functions: logP_typeX gets as input f, parameters,
## the parameter prior, and 1/sigma and returns -log likelihood
## Depending on f(), and on settings for model variance,
## There are 6 different types: 1, 2, 2b, 3, 4, 5
## =============================================================================

# type 1. f() returns -2 * log likelhood, there is no model variance
logP_type1 <- function(f, p, PPnew, divsigma, ...) {
  SSnew  <- f(p, ...)
  list(pnew = 0.5*(SSnew + PPnew), SSnew = SSnew)     #-log(p_model*p_params)
}

# type 2. f() returns residuals, model variance is specified, only one value
logP_type2 <- function(f, p, PPnew, divsigma, ...) {
  SSnew  <- f(p, ...)
  SSnew  <- sum(SSnew^2)
  list(pnew = 0.5*(sum(SSnew*divsigma) + PPnew), SSnew = SSnew)
}

# type 2b. f() returns residuals, model variance is specified, many values
logP_type2b <- function(f, p, PPnew, divsigma, ...) {
  SSnew  <- f(p, ...)
  SSnew  <- SSnew^2
  list(pnew = 0.5*(sum(SSnew*divsigma) + PPnew), SSnew = SSnew)
}

# type 3. f() returns instance of class modcost, model variance not specified
logP_type3 <- function(f, p, PPnew, divsigma, ...) {
  SSnew  <- f(p,...)$minlogp    # select -2*log(probability)
  list(pnew=0.5*(SSnew + PPnew), SSnew = SSnew)     #-log(p_model*p_params)
}

# type 4. f() returns instance of class modcost, model variance specified,
# one value per data point
logP_type4 <- function(f, p, PPnew, divsigma, ...) {
  SSnew  <- f(p, ...)
  SSnew  <- (SSnew$residuals$res)^2
  list(pnew=0.5*(sum(SSnew*divsigma) + PPnew), SSnew = SSnew)
}

# type 5. f() returns instance of class modcost, model variance specified,
# one value per observed variable
logP_type5 <- function(f, p, PPnew, divsigma, ...) {
   SSnew  <- f(p, ...)
   SSnew  <- SSnew$var$SSR.unweighted
   list(pnew=0.5*(sum(SSnew*divsigma) + PPnew), SSnew = SSnew)
}

## =============================================================================
## The acceptance functions ...
## =============================================================================

Test <- function(fnew, fold, SSnew, SSold, useSigma, divsigma, PPnew, PPold) {
    if (is.na(fnew) | is.infinite(fnew))
       return(0)
    if (!useSigma)
       test <- exp(fold-fnew)
    else                  # uses SSR, scaled by sigma
       test <- exp(-0.5*sum((SSnew - SSold)*divsigma) - 0.5*(PPnew - PPold))
    return(test)
}

Accept <- function(tst) {
    if (is.nan(tst)|is.na(tst))
      accept <- FALSE
    else if (tst<0)
      accept<-FALSE
    else if (tst>1)
      accept <- TRUE
    else accept <- (runif(1) < tst)

    return(accept)
}

## =============================================================================
## Markov Chain Monte Carlo - main function
## Part of this R-code (the DR-part) is based on the matlab code (Marco Laine)
## as available from http://www.helsinki.fi/~mjlaine/mcmc/
## =============================================================================

modMCMC <- function (f, p, ..., jump = NULL, lower = -Inf, upper = +Inf,
                     prior = NULL, var0 = NULL, wvar0 = NULL, n0 = NULL,
                     niter = 1000, outputlength = niter, burninlength = 0,
                     updatecov = niter, covscale = 2.4^2/length(p),
                     ntrydr = 1, drscale  = NULL, verbose = TRUE) {

  ##----------------------------------------------------------------------------
  ## 1. check input
  ##----------------------------------------------------------------------------
  npar   <- length(p)
  pnames <- names(p)

  if (length(lower)!= npar & length(lower)!= 1)
    stop("length of 'lower' should be either 1 or = number of parameters")
  if (length(upper)!= npar & length(upper)!=1)
    stop("length of 'upper' should be either 1 or = number of parameters")

  lower[is.na(lower)] <- -Inf
  upper[is.na(upper)] <- Inf
  limits   <- (any(lower != -Inf) || any(upper != Inf))
  if (any(p < lower) | any(p > upper))
    stop("values of 'p' should be inbetween 'lower' and 'upper'")

  ## iterations, burnin, output
  if (outputlength > (niter-burninlength))
    outputlength <- niter-burninlength

  ## recalculate the number of iterations ~ burnin and output...
  if (burninlength >= niter) stop("burninlenght is larger than niter")

  ou <- ceiling((niter - burninlength)/outputlength)  # output interval
  niter <- niter - (niter - burninlength) %% ou
  outputlength <- (niter - burninlength) %/% ou
  ou1 <- burninlength + ou                          # first output iteration nr

  ## delayed rejection input
  ntrydr <- max(1, ntrydr)
  if(is.null (drscale)) drscale=c(1/5, 1/4, 1/3)
  ldr <- length(drscale)
  if (ldr < ntrydr) drscale <- c(drscale, rep(drscale[ldr], ntrydr - ldr))

  ## adaptive metropolis input
  if (updatecov == 0) updatecov <- niter
  isR       <- updatecov < niter  # if TRUE: will update the jump covariances
  nupdate   <- updatecov          # first step to update covariance matrix

  ##----------------------------------------------------------------------------
  ## 2. Function to generate new parameter values
  ##----------------------------------------------------------------------------

  ## NewParsMn used if parameters are generated with multiD normal distribution
  ## (this may not be the case at first, but will be if updatecov < niter
  NewParsMN <- function(p, R) {
    pp <- as.vector(p + rnorm(npar)%*%R) # R is cholesky decomposition
    names(pp) <- pnames
    return(pp)
  }

  ## Check jump input (to generate new parameters)
  if (is.null(jump)) {                     # default is a sd of 0.1
    jump  <- 0.1*abs(p)
    NewPars <- function(p, ...) p + rnorm(npar, sd = jump)
    R <- diag(jump)
  } else if (is.matrix(jump)) {            # covariance matrix
    R       <- chol(jump)
    NewPars <- NewParsMN
  } else if (is.numeric(jump)) {           # a vector or value
    NewPars <- function(p, ...) p + rnorm(npar, sd=jump)
    ii <- which (jump<0)
    jump[ii] <- 0.1*abs(p[ii])
    R <- diag(nrow = npar, jump)
  } else {
    NewPars <- jump                        # jump is a function
    if (ntrydr > 1) stop ("cannot combine jump function with delayed rejection")
  }

  ##----------------------------------------------------------------------------
  ## 3. parameter prior, absent or function that returns -2 * log(prior par prob)
  ##----------------------------------------------------------------------------

  if (is.null(prior))
    Prior <- function(p) return(0)
  else   # uninformative
    Prior <- function(p) return(prior(p))

  PPnew  <- Prior(p)

  ##----------------------------------------------------------------------------
  ## 4. initialise model variance
  ##----------------------------------------------------------------------------

  ## if 'var0' has a value, it is used to scale the SSR; divsigma = 1/var
  useSigma <- !is.null(var0)

  if (useSigma)
    divsigma <- 1/var0
  else
    divsigma <- 1

  ## if 'wvar0' or n0 has a value, then sigma will be updated
  updateSigma <- (useSigma & !is.null(wvar0) | (useSigma & !is.null(n0)))

  ## Number of elements in var0
  lenvar0 <- length(var0)

  ##----------------------------------------------------------------------------
  ## 5. initialise the log likelihood function (logP) and the initial SSR (SSold)
  ## There are 6 types, depending on model variance and function 'f' returning
  ## either -2*log Probability, the model residuals, or a list of class modFit.
  ##----------------------------------------------------------------------------

  ## Call function for first time ...
  SSnew <- f(p, ...)

  ## type I: numeric value, no model variance
  if (is.numeric(SSnew) & !useSigma) {
    N <- length(SSnew)
    if (N>1)
      stop("if 'var0' is NULL, 'f' should either return -2*log of the model probability, OR an instance of class 'modCost'")
    logP  <- function (...) logP_type1(...)
    SSold <- SSnew
  }
  ## type II: numeric value, + model variance -> f returns model residuals
  else if (is.numeric(SSnew)& useSigma) {
    N <- length(SSnew)
    if (N == 1)
      stop("if 'var0' has a value, 'f' should return the model residuals OR an instance of class 'modCost'")
    if (lenvar0 != 1 && lenvar0 != N)
      stop("length of 'var0' should either 1 or = length of model residuals as returned by 'f'")

    if (lenvar0 == 1) {
      logP <- function(...)  logP_type2(...)
      SSold <- sum(SSnew^2)
    } else  { # = N
      logP <- function(...)  logP_type2b(...)
      SSold <- SSnew^2
    }
  } else if (class(SSnew) != "modCost")
    stop("'f' should either return the -2*log of the model probability OR an instance of class 'modCost'")
  ## type III: class modCost & ! usesigma
  else if (!useSigma)  {
    N <- nrow(SSnew$residuals)     # total number of data points
    logP <- function(...)  logP_type3(...)
    SSold <- SSnew$minlogp
  } else {
  ## type IV and V: class modCost & sigma. Either uses residuals or variable SSR.
    if (lenvar0 == nrow(SSnew$residuals)) { # use data residuals
      N <- nrow(SSnew$residuals)   # total number of data points
      logP <- function(...)  logP_type4(...)
      SSold <- (SSnew$residuals$res)^2
    } else
    if (lenvar0 != 1 && lenvar0 != nrow(SSnew$var))
      stop("function 'f' is not compatible with length of var0")
    else {
      N <- SSnew$var$N             # number of data points per variable
      logP <- function(...)  logP_type5(...)
      SSold <- SSnew$var$SSR.unweighted
    }
  }

  ## parameter used in gamma draw.
  if (is.null(n0)) n0 <- wvar0*N

  if (updateSigma)
      divsigma <- rgamma(lenvar0, shape = 0.5*(n0 + N),
         rate = 0.5*(n0*var0 + SSold))

  ##----------------------------------------------------------------------------
  ## initialise some quantities
  ##----------------------------------------------------------------------------

  parold  <- p         # parameter values
  PPold   <- PPnew     # prior parameter prob
  sigold  <- divsigma  # 1/model variance
  funpold <- 0.5 * (sum(SSold*divsigma) + PPold)
  SSnew   <- Inf       # sum of squared residuals

  ##----------------------------------------------------------------------------
  ## The delayed rejection procedure ...  (based on matlab code)
  ##----------------------------------------------------------------------------

  A_count  <- 0       # counter to number of alpha entrances
  dr_steps <- 0       # counter to number of dr steps...


  AlphaFun <- function(arg) {
  #----------------------------------------------------------------------------
  # recursive acceptance function for delayed rejection
  # arg$pri: the prior for the parameters
  # arg$ss : the sum of squares
  #----------------------------------------------------------------------------

    stage <- length(arg)-1      # the stage we are in
    A_count <<- A_count+1

    ## recursively compute past alphas
    a1 <- 1
    a2 <- 1

    if (stage > 1) {
      for (k in 1:(stage-1)) {
        a1 <- a1*(1-AlphaFun(arg[1:(k+1)]))

        a2 <- a2*(1-AlphaFun(arg[seq(from=stage+1, to=stage+1-k, by=-1)]))

        if  (is.null(a2) | is.na(a2))
          return (0)
        if  (a2==0)
          return (0)
      }
    }
    ss2  <- arg[[stage+1]]$ss
    ss1  <- arg[[1]]$ss
    pri2 <- arg[[stage+1]]$pri
    pri1 <- arg[[1]]$pri

    y  <- -0.5*(sum((ss2 - ss1)*divsigma) + (pri2 - pri1))

    for (k in 1:stage)
      y <- y + qfun(k, arg)

    y <- min(1, exp(y)*a2/a1)

    return(y)
  }

  ## gaussian n-th stage log proposal ratio
  ## log of q_i(yn, ... yn-j)/q_i(x, y1, ...yj)

  qfun <- function(iq, arg) {
    LL <- function(x) sum(x*x)

    stage <- length(arg)-1
    if (stage == iq)              # we are symmetric...
      return(0)
    else {
      iR <- invR[[iq]]
      y1 <- arg[[1]]$p            # y1
      y2 <- arg[[iq+1]]$p         # yi
      y3 <- arg[[stage+1]]$p      # yn
      y4 <- arg[[stage-iq+1]]$p   # yn - i
      z <- -0.5*(LL((y4-y3) %*% iR) - LL((y2-y1) %*% iR))
    }
    return(z)
  }

  ##----------------------------------------------------------------------------
  ## the MCMC jumps...
  ##----------------------------------------------------------------------------

  ## matrices/vectors with results
  pars      <- matrix(nrow = npar, ncol = outputlength)    # parameter values
  SSpars    <- vector(length = outputlength)               # SS value
  priorpars <- vector(length = outputlength)               # parameter priors

  if (useSigma)
    sig     <- matrix(nrow = outputlength, ncol = lenvar0) # model variances
  else sig<-NULL

  ## best function and parameter values will be kept
  bestPar   <- parold
  bestfunp  <- funpold

  ## keep track of ...
  naccepted <- 0               # number of saved parameters
  ii        <- 0               # counter to saved parameters
  naccsave  <- 0               # number of saved accepted parameters
  ipos      <- 0
  icov      <- 0
  Rold      <- 0               # these two to know when to update the covariances...
  Rnew      <- 1

  # keep initial runs during burnin, if adaptive metropolis algorithm (isR)
  if (isR & burninlength>0) {
    SaveIni      <- TRUE
    naccsave2    <- 0
    burnupdate   <- updatecov
  } else SaveIni <- FALSE

  # if delayed rejection
  if (ntrydr > 1) {
    Rw   <- list()      # The cholesky decompositions used
    invR <- list()      # The inverses of the choleski decompositions
    trypath <- list()   # The parameter tries
    x <- list()         # previously accepted
    y <- list()         # rejected
    z <- list()         # newly tried
  }

  # loop over iterations ....
  for ( i in 1:niter) {

    accept <- FALSE

    parnew <- NewPars(parold, R)       # new parameter set
    if (any(parnew < lower) | any(parnew > upper)) {
      accept <- FALSE
      funpnew <- Inf

    } else {
      PPnew   <- Prior(parnew)
      logfun  <- logP(f, parnew, PPnew, divsigma, ...)
      funpnew <- logfun$pnew           # probability of new parameter set
      SSnew   <- logfun$SSnew
        if (is.infinite(funpnew)){
          alpha  <- 0
          accept <- FALSE
        } else  {
          if (! is.na(funpnew))
            if (funpnew < bestfunp)  { # new best value - save it!
               bestfunp <- funpnew
               bestPar  <- parnew
            }
          alpha  <- Test(funpnew, funpold, SSnew, SSold, useSigma, divsigma, PPnew, PPold)
          accept <- Accept(alpha)
          if (accept) {
            SSold  <- SSnew
            PPold  <- PPnew
            if (useSigma) sigold <- divsigma
          }
        }
      }

    # start delayed rejection
    itry <- 1

    if (ntrydr > 1 & !accept) {
      x$p   <- parold
      x$ss  <- SSold
      x$pri <- PPold

      y$p   <- parnew
      y$ss  <- SSnew
      y$pri <- PPnew

      trypath[[1]] <- x
      trypath[[2]] <- y

      ## Create new R's and the inverse of R, scaled - ONLY IF R HAS CHANGED
      if (Rnew != Rold) {
        invR[[1]] <- ginv(R)
        Rw[[1]]   <- R
        for (j in 2:ntrydr) {
          Rw[[j]]   <- Rw[[j-1]] * drscale[j-1]
          invR[[j]] <- ginv(Rw[[j]])
        }
        Rold <- Rnew
      }

      while (! accept & itry < ntrydr) {
        itry <- itry + 1

        parnew  <- NewPars(parold, Rw[[itry]])
        z$p <- parnew

        if (any(parnew < lower) | any(parnew > upper)) {
          funpnew <- Inf
          z$pri   <- 0
          z$ss    <- Inf
        } else {
          PPnew   <- Prior(parnew)
          logfun  <- logP(f, parnew, PPnew, divsigma, ...)
          funpnew <- logfun$pnew        # probability of new parameter set
          SSnew   <- logfun$SSnew
          z$pri   <- PPnew
          z$ss    <- SSnew
        }
        trypath[[itry+1]] <- z

        if (is.infinite(funpnew)) {
          accept <- FALSE
          break
        }
        else {
          alpha <- AlphaFun(trypath[1:(itry + 1)])
          dr_steps <- dr_steps +1
          accept   <- Accept(alpha)
          if (accept) {
            SSold <- SSnew
            PPold <- PPnew
            if (useSigma) sigold <- divsigma
          }
        }
      }
    } # end delayed rejection

    ## saving if accepted or during burnin in case of adaptive metropolis...
    if (accept) {
      naccepted <- naccepted +1
      parold    <- parnew           # new parameter set replaces the old one...
      funpold   <- funpnew
    }

    if(updateSigma)
      divsigma <- rgamma(lenvar0, shape = 0.5*(n0 + N), rate = 0.5*(n0*var0 + SSold))

    ## during burnin: SaveIni will be true
    if (SaveIni) {  # use burnin to update covariance...
      ipos <- ipos + 1
      if (ipos > outputlength) {
        ipos <- 1
        icov <- outputlength
      }
      pars[,ipos] <-parold

      if (accept) naccsave2 <- naccsave2 + 1

      if (i > burnupdate & naccsave2 > 5 ) {  #  5 is arbitrary...
        jj <- max(ipos, icov)
        if (npar > 1)
          Covar <- (cov(t(pars[,1:jj])) + diag(1e-16,npar)) * covscale
        else Covar <- var(pars[,1:jj]) + 1e-16 * covscale
        RR    <- try(chol(Covar))
        if (is.numeric(RR) ) {
          R <- RR
          NewPars    <- NewParsMN   # this function becomes the default update function
          burnupdate <- i + updatecov
          Rnew  <- Rnew + 1
        }
      }
    }

    ## past burnin: current parameter set saved every nupdate steps ...
    if (i == ou1) {
      if (SaveIni)
        isR <- FALSE
      SaveIni       <- FALSE
      ii            <- ii + 1
      pars[,ii]     <- parold
      SSpars[ii]    <- sum(SSold)
      priorpars[ii] <- PPold
      if (accept)
        naccsave <- naccsave + 1
      ou1 <- i + ou

      if (useSigma)
        sig[ii,] <- 1/sigold

      ## update jump parameters by estimating param covariance.. only if > 5 accepted ones!
      if (isR & ii > nupdate & naccsave > 5) {
        ie <- max(ipos, icov, ii)
        if (npar > 1)
          Covar <- (cov(t(pars[, 1:ie])) + diag(1e-16, npar)) * covscale
        else Covar <- var(pars[, 1:ie]) + 1e-16 * covscale
        RR     <- try(chol(Covar))

        if (is.numeric(RR)) {
          R <- RR
          NewPars <- NewParsMN   # this function becomes the default update function
          nupdate <- ii + updatecov
          Rnew <- Rnew + 1
        }
      }
    }
  } # end of loop

  ##----------------------------------------------------------------------------
  ## 3. finalisation
  ##----------------------------------------------------------------------------

  if (is.null(pnames))
    pnames <- paste("p", 1:npar, sep="")
  if (useSigma) {
    if (length(var0) == ncol(sig))
      colnames(sig) <- names(var0)
    else
      colnames(sig) <- paste("var", 1:ncol(sig), del = "")
  }

  if (updateSigma)
     settings <- data.frame(var0 = var0, n0 = n0, N = N)
  else if (useSigma)
     settings <- paste("constant variance =", var0)
  else
     settings <- "Constant variance = 1"

  rownames(pars) <- pnames

  if (verbose)
     cat("number of accepted runs: ", naccepted,
            " out of ", niter, " (", 100 * naccepted/niter,
              "%) \n", sep = "")

  count <- c(dr_steps = dr_steps, Alfasteps = A_count, num_accepted = naccepted,
             num_covupdate = Rnew - 1)

  if (! is.null(sig)) {
    if (is.null(colnames(sig))) colnames(sig) <- rep("model", ncol(sig))
    colnames(sig) <- paste("var_", colnames(sig), sep = "")
  }
  res <- list(pars = t(pars), SS = SSpars, naccepted = naccepted, sig = sig,
             bestpar = bestPar, bestfunp = bestfunp, prior = priorpars,
             count = count, settings = settings)

  class(res) <- "modMCMC"
  return(res)
}

## =============================================================================
## S3 methods of modMCMC
## =============================================================================

pairs.modMCMC <- function (x, Full = FALSE, which = 1:ncol(x$pars),
                           remove = NULL, nsample = NULL, ...) {

  panel.main <- function(x, y, ...)
    points(x[ii], y[ii], ...)

  var <- colnames(x$pars)
  which <- selectvar(which, var, "x$pars", Nall = FALSE)

  X <- x$pars[, which]


  if (is.vector(X)) X<- as.matrix(X)
  if (is.null(nsample))
    ii <- 1:nrow(X) else
    ii <- sample((1:nrow(X)), nsample)
  if (Full)
    X <- cbind(X, SSR = x$SS)
  if (Full & !is.null(x$sig))
    X <- cbind(X, x$sig)

  if (! is.null (remove)) {
    if (max(remove) > nrow(X))
      stop("too many runs should be removed from modMCMC object")
    if (min(remove) < 1)
      stop("cannot remove negative runs from modMCMC object")
    X <- X[-remove,]
  }


  labels <- colnames(X)
  dots <- list(...)
  dotnames <- names(dots)
  
  if(! "diag.panel" %in% dotnames)
    dots$diag.panel  <- panel.hist
  if(! "lower.panel" %in% dotnames)
    dots$lower.panel <- panel.cor
  if (! "upper.panel" %in% dotnames)
    dots$upper.panel <- panel.main
    
  dots$gap <- if(is.null(dots$gap)) 0 else dots$gap
  dots$labels <- if(is.null(dots$labels)) labels else dots$labels

  do.call("pairs", c(alist(X), dots))
}

## -----------------------------------------------------------------------------

cumuplot.modMCMC <- function (x, Full = FALSE, which = 1:ncol(x$pars),
                              remove = NULL, ...) {

  var <- colnames(x$pars)
  which <- selectvar(which, var, "x$pars", Nall = FALSE)

  mcmc <- x$pars[, which]
  if (! is.null (remove)) {
    if (max(remove) > nrow(mcmc))
      stop("too many runs should be removed from modMCMC object")
    if (min(remove) < 1)
      stop("cannot remove negative runs from modMCMC object")
    mcmc <- mcmc[-remove,]
  }

  if (Full) mcmc <- cbind(mcmc, x$SS)
  if (Full & !is.null(x$sig))
      mcmc <- cbind(mcmc, x$sig)
  cumuplot(as.mcmc(mcmc), ...)
}

## -----------------------------------------------------------------------------

plot.modMCMC <- function (x, Full = FALSE, which = 1:ncol(x$pars), trace = TRUE,
                          remove = NULL, ask = NULL, ...) {

  var <- colnames(x$pars)
  which <- selectvar(which, var, "x$pars", Nall = FALSE)

  np <- NP <- length(which)
  if (Full)
    np <- np +1
  if (Full & !is.null(x$sig))
    np <- np +ncol(x$sig)

  dots <- list(...)
  nmdots <- names(dots)

  ## Set par(mfrow) and ask.
  ask <- setplotpar(nmdots, dots, np, ask)

  ## interactively wait if there are remaining figures
  if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
  }

  mcmc <- x$pars
  if (! is.null(remove)) {
    if (max(remove) > nrow(mcmc))
      stop("too many runs should be removed from modMCMC object")
    if (min(remove) < 1)
      stop("cannot remove negative runs from modMCMC object")
    mcmc <- mcmc[-remove,]
  }
  Main <- is.null(dots$main)

  dots$xlab <- if(is.null(dots$xlab)) "iter" else dots$xlab
  dots$ylab <- if(is.null(dots$ylab)) "" else dots$ylab
  dots$type <- if(is.null(dots$type)) "l" else dots$type

  for(i in which) {
    if (Main) dots$main <- colnames(mcmc)[i]
    do.call("plot", c(alist(mcmc[,i]), dots))
    if (trace) lines(lowess(mcmc[,i]), col = "darkgrey", lwd = 2)
  }

  if (Full) {
    plot(x$SS, type="l", main = "SSR", xlab = "iter", ylab = "")
    if (trace) lines(lowess(x$SS), col = "darkgrey", lwd = 2)
  }

  if (Full & !is.null(x$sig))
    for ( i in 1:ncol(x$sig)) {
       Name <- if(!is.null(colnames(x$sig)))colnames(x$sig)[i]
         else
           { if (ncol(x$sig) == 1) "Model" else paste("variable", i)}
       plot(x$sig[, i], type = "l", main = Name,
                    xlab = "iter", ylab = "variance", log = "y")
       if (trace)
         lines(lowess(x$sig[, i]), col = "darkgrey", lwd = 2)
    }
}

## -----------------------------------------------------------------------------

hist.modMCMC <- function (x, Full=FALSE, which=1:ncol(x$pars),
                          remove=NULL, ask = NULL, ...) {

  iswhat <- TRUE
  if (length(which) > 0)
    if (any(is.na(which))) iswhat <- FALSE
  if (iswhat & is.null(which)) iswhat <- FALSE
  np <- NP <- 0
  var <- colnames(x$pars)

  if (iswhat)  {
    which <- selectvar(which, var, " x$pars", Nall = FALSE)
    np <- NP <- length(which)
  }
  if (Full)
    np <- np +1
  if (Full & !is.null(x$sig))
    np <- np +ncol(x$sig)

  dots <- list(...)
  nmdots <- names(dots)

  ## Set par mfrow and ask.
  ask <- setplotpar(nmdots, dots, np, ask)

  ## interactively wait if there are remaining figures
  if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
  }

  Main <- is.null(dots$main)
  dots$xlab    <- if(is.null(dots$xlab))    ""    else dots$xlab
  dots$freq    <- if(is.null(dots$freq))    FALSE else dots$freq
  dots$ylab    <- if(is.null(dots$ylab))    "-"   else dots$ylab
  dots$breaks  <- if (is.null(dots$breaks)) 100   else dots$breaks

  mcmc <- x$pars
  if (! is.null (remove)) {
    if (max(remove) > nrow(mcmc))
      stop("too many runs should be removed from modMCMC object")
    if (min(remove) < 1)
      stop("cannot remove negative runs from modMCMC object")
    mcmc <- mcmc[-remove,]
  }

  if (iswhat)
    for(i in which) {
      if (Main) dots$main <- colnames(mcmc)[i]
      do.call("hist", c(alist(mcmc[, i]), dots))
  }
  if (Full) {
    dots$main <- "SSR"
    do.call("hist", c(alist(x$SS), dots))
  }

  if (Full & !is.null(x$sig)) {
    for ( i in 1:ncol(x$sig))  {
      if (Main) dots$main <- paste("error", colnames(x$sig)[i])
      dots$ylab <- colnames(x$sig)[i]
      do.call("hist", c(alist(x$sig[, i]), dots))
    }
  }
}

## -----------------------------------------------------------------------------

summary.modMCMC <- function (object, remove = NULL, ...) {

  mcmc <- object$pars
  if (! is.null (remove)) {
    if (max(remove) > nrow(mcmc))
      stop("too many runs should be removed from modMCMC object")
    if (min(remove) < 1)
      stop("cannot remove negative runs from modMCMC object")
    mcmc <- mcmc[-remove,]
  }
  Res <- data.frame(rbind(
    mean = apply(mcmc, 2, mean),
    sd   = apply(mcmc, 2, sd),
    min  = apply(mcmc, 2, min),
    max  = apply(mcmc, 2, max),
    q025 = apply(mcmc, 2, quantile, probs = 0.25),
    q050 = apply(mcmc, 2, quantile, probs = 0.5),
    q075 = apply(mcmc, 2, quantile, probs = 0.75)
  ))

  if (!is.null(object$sig))
    Res <- data.frame(Res,
      sig = rbind(
        mean = apply(object$sig, 2, mean),
        sd   = apply(object$sig, 2, sd),
        min  = apply(object$sig, 2, min),
        max  = apply(object$sig, 2, max),
        q025 = apply(object$sig, 2, quantile, probs = 0.25),
        q050 = apply(object$sig, 2, quantile, probs = 0.5),
        q075 = apply(object$sig, 2, quantile, probs = 0.75)
    ))
  return(Res)
}
