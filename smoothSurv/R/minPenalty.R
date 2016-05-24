#############################################
#### AUTHOR:    Arnost Komarek           ####
####            (23/07/2004)             ####
####                                     ####
#### FILE:      minPenalty.R             ####
####                                     ####
#### FUNCTIONS: minPenalty.R             ####
#############################################

### ====================================================================================
### minPenalty: minimize the penalty term under the constraints 
### ====================================================================================
minPenalty <- function(knots = NULL,
                       dist.range = c(-6, 6),
                       by.knots = 0.3,
                       sdspline = NULL,
                       difforder = 3,
                       init.c,
                       maxiter = 200,                       
                       rel.tolerance = 1e-10,
                       toler.chol = 1e-15,
                       toler.eigen = 1e-3,
                       maxhalf = 10,
                       debug = 0,
                       info = TRUE)
{
  est.c <- TRUE
    
  ### Create knots and other parameters that control the fit
  ### --------------------------------------------------------
  pars <- smoothSurvReg.control(est.c = est.c, est.scale = FALSE, maxiter = maxiter, firstiter = 0,
                                rel.tolerance = rel.tolerance, toler.chol = toler.chol, toler.eigen = toler.eigen, maxhalf = maxhalf,
                                debug = debug, info = info, lambda.use = 1.0, sdspline = sdspline, difforder = difforder,
                                dist.range = dist.range, by.knots = by.knots, knots = knots, nsplines = NULL, last.three = NULL)

  
  ### Initial values for c coefficients
  ### ----------------------------------
  if (missing(init.c) || is.null(init.c)){
     ## try to approximate normal distribution
     best.dens <- "dnorm"
     init.c <- find.c(pars$knots, pars$sdspline, best.dens)
     if (!is.null(init.c)){
        i1 <- which.max(init.c)
        i2 <- which.max(init.c[-i1]); i2 <- ifelse(i2 < i1, i2, i2 + 1)
        i3 <- which.max(init.c[-c(i1, i2)]); i3 <- ifelse(i3 < min(i1, i2), i3, ifelse(i3 < max(i1, i2) - 1, i3 + 1, i3 + 2))
        last.three.temp <- c(i1, i2, i3)
        init.c <- give.c(pars$knots, pars$sdspline, last.three.temp, init.c[-last.three.temp])
        init.c[init.c < 1e-5] <- 1e-5
     }
     else{
        ## USE ANOTHER METHOD ==> LATER ON (MAYBE)
        stop("Singularity when computing initial c's, try to give your own initial c's or a's  ")
     }
  }
  else{
     if(length(init.c) != pars$nsplines) stop("Incorrect length of the vector of initial c coefficients. ")
     if((sum(init.c) < 0.99) || (sum(init.c) > 1.01)) stop("Sum of initial c coefficients is not 1. ")
     if((sum(init.c < 0) > 0) || (sum(init.c > 1) > 0)) stop("All c coefficients must be between 0 and 1. ")
     init.c <- give.c(pars$knots, pars$sdspline, pars$last.three, init.c[-pars$last.three])
     init.c[init.c < 1e-5] <- 1e-5
  }
  acoef <- c2a(init.c, pars$last.three[1])  

  ### Optimize the penalty
  ### ---------------------
  nknots <- pars$nsplines
  
    ## Dimension of d's (g - 3 or 0) to be estimated
  nUa <- ifelse(est.c, nknots - 3, 0)

    ## Number of parameters to be estimated
  nparam <- nUa

    ## Size of dCdD
  ndcdd <- ifelse(est.c, nknots * (nknots - 3), 1)

    ## Size of matrices used to compute df
  ndfm <- ifelse(est.c, nknots - 1, 1)

  fit <- .C("smoothSurvReg84",
                  as.integer(0),                 # n
                  as.integer(0),                 # nY
                  as.integer(0),                 # nX
                  as.integer(0),                 # nZ
                  as.integer(pars$nsplines),
                  as.double(0),                  # matrix X
                  as.double(0),                  # matrix Y 
                  as.double(0),                  # offset
                  as.double(0),                  # matrix Z
                  as.double(pars$knots),         # original sequence of knots (also on output)
                  as.double(pars$sdspline),
   lastThree =    as.integer(pars$last.three - 1),  # C++ indeces of a coefficients which are expressed as the function of the remaining ones
                  as.integer(pars$est.scale),
                  as.integer(pars$est.c),
   beta =         as.double(0),                  # initial beta
   logscale =     as.double(0),                  # initial gamma
   acoef =        as.double(acoef),              # on OUTPUT: all a coefficients (ZERO's included)
   ccoef =        double(pars$nsplines),         # on OUTPUT: c's corresponding to a's
   penalloglik =  double(1),
   loglik =       double(1),
                  as.double(0),                  # correction to likelihood
   penalty =      double(1),
   H =            double(nparam*nparam),         # minus Hessian of the penalized log-likelihood
   I =            double(nparam*nparam),         # minus Hessian of the un-penalized log-likelihood
   G =            double(nparam*nparam),         # minus Hessian of the penalty term => H = I + G
   U =            double(nparam),                # score vector at the convergence
   dCdD =         double(ndcdd),                 # s of c's w.r.t. d's
   Ha =           double(ndfm*ndfm),
   Ia =           double(ndfm*ndfm),
   Ga =           double(ndfm*ndfm),
   dCon =         double(2 * ndfm),
                  as.double(pars$lambda.use),    # lambda
                  as.integer(pars$difforder),
   iter =         as.integer(pars$maxiter),
                  as.integer(pars$firstiter),
                  as.double(pars$rel.tolerance),
                  as.double(pars$toler.chol),
                  as.double(pars$toler.eigen),
                  as.integer(pars$maxhalf),
                  as.integer(pars$info),
                  as.integer(pars$debug),
   fail =         integer(1),
   nonPosDefH =   integer(1),
  PACKAGE = "smoothSurv"
  )
  
  warn <- ""
  if (fit$fail >= 99){
        warn <- "No fit is produced "
        temp <- list(fail = fit$fail)
        return(temp)
  }


## Warnings concerning the convergence
  warn.num <- fit$fail %% 10
  warn <- switch(warn.num + 1,
             "OK",
             "OK",
             "OK",
             "H not positive definite and eigen value decomposition failed",
             "Not converging, not able to increase the objective function",
             "Not possible to find the reference knots",
             "Ran out of iterations and did not converge"
          )

## Print warnings
  if (warn != "OK"){
      warning(paste(warn, " ", sep = ""))
      warn <- paste(warn, ".", sep = "")
  }
  fail.num <- warn.num

  warn.all <- data.frame(c(warn))
  rownames(warn.all) <- c("Convergence")
  colnames(warn.all) <- "Warnings"

  
## Labels for c coefficients
  ind.d <- (1:nknots)[-fit$lastThree]  
  cname <- paste("c(",pars$knots,")", sep="")                         ## all c's
  aname <- paste("a(",pars$knots,")", sep="")                         ## all a's
  if (est.c) dname <- aname[ind.d]
  else       dname <- NULL
  names(fit$ccoef) <- cname                        # these are possibly fixed c's
  names(fit$acoef) <- aname                        # these are possibly fixed c's

  knotname <- paste("knot[",1:nknots,"]", sep = "")

## Basis spline SD (normal density)
  sd.spline <- rep(sdspline, nknots)


## Put all spline information into a dataframe
  ccoef <- fit$ccoef
  acoef <- fit$acoef
  names(ccoef) <- cname
  spline <- data.frame(Knot = pars$knots, SD.spline = sd.spline,
                       c.coef = ccoef,
                       a.coef = acoef)
  colnames(spline) <- c("Knot", "SD basis", "c coef.", "a coef.")
  rownames(spline) <- knotname

## Resulting object
  temp <- list(spline = spline,
               penalty = fit$penalty,
               iter = fit$iter,
               warning = warn.all,
               fail = fail.num
               )
  
  return(temp)      
}  
