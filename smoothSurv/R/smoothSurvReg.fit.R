###########################################
#### AUTHOR:    Arnost Komarek         ####
####            29/04/2004             ####
####                                   ####
#### FILE:      smoothSurvReg.fit.R    ####
####                                   ####
#### FUNCTIONS: smoothSurvReg.fit      ####
####            MP.pseudoinv           ####
###########################################

### =========================================================================
### smoothSurvReg.fit: Survival regression with smoothed error distribution,
###                      fitter used inside smoothSurvReg
### =========================================================================
## on OUTPUT: fail = 0,   everything OK
##                 = 3,   not converging because I am not able to make H positive definite
##                 = 4,   not converging because of too many half-steps
##                 = 5,   not possible to find reference knots
##                 = 6,   not converging because of too many iterations
##                 +10     final H is not positive definite
##                 +20     df <= 0 or not defined because I do not have H^{-1}
smoothSurvReg.fit <- function(x, z, y, offset = NULL, correctlik,
                              init, controlvals, common.logscale)
{
    ## Get a list of control values for iteration process.
    controlvals <- do.call("smoothSurvReg.control", controlvals)
    est.c <- controlvals$est.c
    est.scale <- controlvals$est.scale
    maxiter <- controlvals$maxiter
    firstiter <- controlvals$firstiter
    eps <- controlvals$rel.tolerance
    tolChol <- controlvals$toler.chol
    tolEigen <- controlvals$toler.eigen
    maxhalf <- controlvals$maxhalf
    debug <- controlvals$debug
    info <- controlvals$info
    lambda.use <- controlvals$lambda.use
    sdspline <- controlvals$sdspline
    difforder <- controlvals$difforder
    knots <- controlvals$knots
    nknots <- controlvals$nsplines
    last.three <- controlvals$last.three

    ## Initial values
    ## !!! I do not check whether correctly given (check is performed in smoothSurvReg())
    ## !!! I do not assume this function will be used directly by the user
    ## !!! If the user wishes to use it, it is his/her responsibility to give a proper list of initials
    beta <- init$beta
    gama <- init$logscale
    ccoef <- init$ccoef
    acoef <- c2a(init$ccoef, last.three[1])

    ## Design matrix (usually created in smoothSurvReg())
    if (!is.matrix(x)) stop("Invalid x matrix ")
    n <- nrow(x)
    namesx <- dimnames(x)[[2]]
    nvarx <- ncol(x)

    ## Design matrix for log(scale) (usually created in smoothSurvReg())
    if (!is.matrix(z)) stop("Invalid z matrix ")
    nz <- nrow(x)
    if (nz != n) stop("x and z matrices have different number of rows ")
    namesz <- dimnames(z)[[2]]
    nvarz <- ncol(z)
    
    ## This will determine a type of response (later).
    ## 3 columns for interval censored data, 2 columns for only right/left censored data.
    if (!is.matrix(y)) stop("Invalid y matrix ")
    ny <- ncol(y)
    if (dim(y)[1] != n) stop("Invalid y matrix ")

    ## Offset term.
    if (is.null(offset)) offset <- rep(0,n)

    ## Boolean => integer
    eest.scale <- 1*est.scale
    eest.c <- 1*est.c
    iinfo <- 1*info

    ## Dimensions of regression parameters (beta + scale) to be estimated
    nUregres <- ifelse(est.scale, nvarx + nvarz, nvarx)
    nestScale <- ifelse(est.scale, nvarz, 0)

    ## Dimension of d's (g - 3 or 0) to be estimated
    nUa <- ifelse(est.c, nknots - 3, 0)

    ## Number of parameters to be estimated
    nparam <- nUregres + nUa

    ## Size of dCdD
    ndcdd <- ifelse(est.c, nknots * (nknots - 3), 1)

    ## Size of matrices used to compute df
    ndfm <- ifelse(est.c, nknots - 1, 1)

    if (nparam == 0) stop("Nothing to be estimated... ")   # this should never occure but one never knows...

    ## Fit the model
    fit <- .C("smoothSurvReg84",
                     as.integer(n),
                     as.integer(ny),
                     as.integer(nvarx),
                     as.integer(nvarz),
                     as.integer(nknots),
                     as.double(x),
                     as.double(y),
                     as.double(offset),
                     as.double(z),
                     as.double(knots),              # original sequence of knots (also on output)
                     as.double(sdspline),
      lastThree =    as.integer(last.three - 1),    # C++ indeces of a coefficients which are expressed as the function of the remaining ones
                     as.integer(eest.scale),
                     as.integer(eest.c),
      beta =         as.double(beta),
      logscale =     as.double(gama),
      acoef =        as.double(acoef),              # on OUTPUT: all a coefficients (ZERO's included)
      ccoef =        double(nknots),                # on OUTPUT: c's corresponding to a's
      penalloglik =  double(1),
      loglik =       double(1),
                     as.double(correctlik),
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
                     as.double(lambda.use),
                     as.integer(difforder),
      iter =         as.integer(maxiter),
                     as.integer(firstiter),
                     as.double(eps),
                     as.double(tolChol),
                     as.double(tolEigen),
                     as.integer(maxhalf),
                     as.integer(iinfo),
                     as.integer(debug),
      fail =         integer(1),
      nonPosDefH =   integer(1),
    PACKAGE = "smoothSurv"
    )

    warn <- ""
    warn2 <- ""
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

## Warnings concerning positive definitness of H
    fit$nonPosDefH <- fit$nonPosDefH > 0
    if (fit$nonPosDefH)  warn2 <- "Final H is not positive definite"
    else                 warn2 <- "OK"


    noPosH <- fit$nonPosDefH
    nodf <- FALSE


## Print warnings
    if (warn != "OK"){
        warning(paste(warn, " ", sep = ""))
        warn <- paste(warn, ".", sep = "")
    }
    if (warn2 != "OK"){
        warning(paste(warn2, " ", sep = ""))
        warn2 <- paste(warn2, ".", sep = "")
    }

    fail.num <- warn.num + 10*noPosH


## Minus Hessian matrix and its components, score vector
    if (fail.num != 5){
       H <- matrix(fit$H[1:(nparam^2)], nrow = nparam)
       I <- matrix(fit$I[1:(nparam^2)], nrow = nparam)
       G <- matrix(fit$G[1:(nparam^2)], nrow = nparam)
       U <- fit$U[1:nparam]
       singH <- FALSE
    }
    else{
       H <- matrix(rep(NA, nparam^2), nrow = nparam)
       I <- matrix(rep(NA, nparam^2), nrow = nparam)
       G <- matrix(rep(NA, nparam^2), nrow = nparam)
       U <- rep(NA, nparam)
       singH <- TRUE
    }


## Compute inversion of H matrix (if it's possible)
    if (!singH){
       Hinv <- qr(H, tol = 1e-07)
       if (Hinv$rank == ncol(Hinv$qr))  Hinv <- solve(Hinv)
       else                             singH <- TRUE
    }

    if (singH)  Hinv <- matrix(rep(NA, nparam^2), nrow = nparam)


## Compute variance matrices
    if (est.c){
       var <- Hinv
       var2 <- Hinv %*% I %*% Hinv
    }
    else{
       var <- Hinv
       var2 <- Hinv
    }


## Compute degrees of freedom
    if (!est.c){
       df <- nparam
    }
    else{
       Ha <- matrix(fit$Ha[1:(ndfm^2)], nrow = ndfm)
       Ia <- matrix(fit$Ia[1:(ndfm^2)], nrow = ndfm)
       Ga <- matrix(fit$Ga[1:(ndfm^2)], nrow = ndfm)
       dCon <- matrix(fit$dCon[1:(2*ndfm)], nrow = ndfm)

    ## Basis for projection and projected Hessians
       qr.K <- qr(dCon)
       q.K <- qr.Q(qr.K, complete = TRUE)
       y.K <- q.K[, 1:ncol(dCon)]
       z.K <- q.K[, (ncol(dCon)+1):ncol(q.K)]
       r.K <- qr.R(qr.K, complete = FALSE)
       Ha.proj <- t(z.K) %*% Ha %*% z.K
       Ia.proj <- t(z.K) %*% Ia %*% z.K

    ## Inversion of the penalized projected Hessian
       singHa <- FALSE
       Hainv <- qr(Ha.proj, tol = 1e-07)
       if (Hainv$rank == ncol(Hainv$qr))  Hainv <- solve(Hainv)
       else                               singHa <- TRUE

    ## Finally, degrees of freedom
       if (singHa){
          nodf <- TRUE
          df <- NA
       }
       else{
          dfMat <- Ia.proj %*% Hainv
          dfspline <- sum(diag(dfMat))
          df <- nUregres + dfspline
          if (df <= 0) nodf <- TRUE
       }
    }

    if (debug > 0 && !nodf){
       cat("Diagonal for spline df: \n")
       cat(diag(dfMat))
       cat("\n")
    }


## Compute degrees of freedom (method 2)
    nodf2 <- FALSE
    if (FALSE){
        if (!est.c){
           df2 <- nparam
        }
        else{
           Hspline <- H[(nUregres+1):nparam, (nUregres+1):nparam]
           Ispline <- I[(nUregres+1):nparam, (nUregres+1):nparam]

    ## Inversion of the spline part of the Hessian
           singHspline <- FALSE
           Hsplineinv <- qr(Hspline, tol = 1e-07)
           if (Hsplineinv$rank == ncol(Hsplineinv$qr))  Hsplineinv <- solve(Hsplineinv)
           else                                         singHspline <- TRUE

    ## Finally, degrees of freedom
           if (singHspline){
              nodf2 <- TRUE
              df2 <- NA
           }
           else{
              dfspline2 <- sum(diag(Ispline %*% Hsplineinv))
              df2 <- nUregres + dfspline2
              if (df2 <= 0) nodf2 <- TRUE
           }
       }
    }

## Warnings about degrees of freedom
    if (nodf){
       fail.num <- fail.num + 20
       warn.df <- "Non positive degrees of freedom"
    }
    else
       warn.df <- "OK"


## Put all warnings together
    warn.all <- data.frame(c(warn, warn2, warn.df))
    rownames(warn.all) <- c("Convergence     ", "Final minus Hessian     ", "df     ")
    colnames(warn.all) <- "Warnings"


## Compute variances of c's (Delta method, one by one), do not compute their cross covariance
## and variance of d's
    var.a <- rep(NA, nknots)
    var2.a <- rep(NA, nknots)
    ind.d <- (1:nknots)[-fit$lastThree]

    if (est.c){
      Hinv.d <- var[(nUregres+1):(nparam), (nUregres+1):(nparam)]
      Hinv2.d <- var2[(nUregres+1):(nparam), (nUregres+1):(nparam)]
      var.d <- diag(Hinv.d)
      var2.d <- diag(Hinv2.d)
      var.d[abs(var.d) < 1e-5] <- 0
      var2.d[abs(var2.d) < 1e-5] <- 0
      var.d[var.d < 0] <- NaN
      var2.d[var2.d < 0] <- NaN

      var.a[ind.d] <- var.d
      var2.a[ind.d] <- var2.d

      dCdD <- matrix(fit$dCdD, nrow = nknots - 3, ncol = nknots)
      Hinv.c <- t(dCdD) %*% Hinv.d %*% dCdD
      Hinv2.c <- t(dCdD) %*% Hinv2.d %*% dCdD

      var.c <- diag(Hinv.c)
      var2.c <- diag(Hinv2.c)
      var.c[var.c < 0] <- NaN
      var2.c[var2.c < 0] <- NaN
    }
    else{
      var.c <- rep(NA, nknots)
      var2.c <- rep(NA, nknots)
      dCdD <- NA
    }

    sd.a <- sqrt(var.a)
    sd2.a <- sqrt(var2.a)

    sd.c <- sqrt(var.c)
    sd2.c <- sqrt(var2.c)


## Variance of the regression part
    var.regres <- diag(var)[1:nUregres]
    var2.regres <- diag(var2)[1:nUregres]
    var.regres[var.regres < 0] <- NaN
    var2.regres[var2.regres < 0] <- NaN

    sd.regres <- sqrt(var.regres)
    sd2.regres <- sqrt(var2.regres)


## Labels for beta's and log(scale) and their sd
    if (!is.null(x)){
         regresname <- namesx
         if (is.null(regresname)) regresname <- paste("x", 1:ncol(x), sep="")
         regresname2 <- regresname
         regres.est <- fit$beta
         names(regres.est) <- regresname2

         if (est.scale){
            if (is.null(namesz)) namesz <- paste("z", 1:nvarz, sep = "")
            if (common.logscale){
              regresname <- c(regresname, "Log(scale)")
              regresname2 <- c(regresname, "Scale")
              regres.est <- c(regres.est, fit$logscale, exp(fit$logscale))
              sd.regres <- c(sd.regres, NA)
              sd2.regres <- c(sd2.regres, NA)
              names(regres.est) <- regresname2
              names(sd.regres) <- regresname2
              names(sd2.regres) <- regresname2              
            }
            else{
              regresname <- c(regresname, paste("LScale.", namesz, sep = ""))
              regres.est <- c(regres.est, fit$logscale)
              names(regres.est) <- regresname
              names(sd.regres) <- regresname
              names(sd2.regres) <- regresname                            
            }                          
         }
         else{
            names(regres.est) <- regresname
            names(sd.regres) <- regresname
            names(sd2.regres) <- regresname
         }
    }


## Beta and log(scale) pars. estimates with sd
    if(nUregres > 0){
      regres <- data.frame(Value = regres.est,
                          Std.Error = sd.regres,
                          Std.Error2 = sd2.regres)
      colnames(regres) <- c("Value", "Std.Error", "Std.Error2")
    }
    else regres <- NULL     ## never possible with this version of the program


## Labels for c coefficients
    cname <- paste("c(",knots,")", sep="")                         ## all c's
    aname <- paste("a(",knots,")", sep="")                         ## all a's
    if (est.c) dname <- aname[ind.d]
    else       dname <- NULL
    names(fit$ccoef) <- cname                        # these are possibly fixed c's
    names(fit$acoef) <- aname                        # these are possibly fixed c's

    knotname <- paste("knot[",1:nknots,"]", sep = "")

    names(sd.c) <- cname
    names(sd2.c) <- cname
    names(sd.a) <- aname
    names(sd2.a) <- aname


## Basis spline SD (normal density)
    sd.spline <- rep(sdspline, nknots)


## Put all spline information into a dataframe
    ccoef <- fit$ccoef
    acoef <- fit$acoef
    names(ccoef) <- cname
    spline <- data.frame(Knot = knots, SD.spline = sd.spline,
                       c.coef = ccoef,
                       Std.Error.c = sd.c,
                       Std.Error2.c = sd2.c,
                       a.coef = acoef,
                       Std.Error.a = sd.a,
                       Std.Error2.a = sd2.a
                       )
    colnames(spline) <- c("Knot", "SD basis",
                  "c coef.", "Std.Error.c", "Std.Error2.c",
                  "a coef.", "Std.Error.a", "Std.Error2.a")
    rownames(spline) <- knotname


## Names for Hessian matrix etc.
    allname <- c(regresname, dname)
    dimnames(H) <- list(allname, allname)
    dimnames(I) <- list(allname, allname)
    dimnames(G) <- list(allname, allname)
    dimnames(var) <- list(allname, allname)
    dimnames(var2) <- list(allname, allname)
    names(U) <- allname
    if (est.c) dimnames(dCdD) <- list(dname, cname)


## Loglikelihood, penalty etc.
    loglik <- data.frame(Log.Likelihood = fit$loglik,
                         Penalty = fit$penalty,
                         Penalized.Log.Likelihood = fit$penalloglik)
    colnames(loglik) <- c("Log Likelihood", "Penalty",
                          "Penalized Log Likelihood")
    rownames(loglik) <- "     "


## Degree of smooth
    degree.smooth <- data.frame(lambda.use, log(lambda.use), df, nparam, 
                                nUregres-nestScale, nestScale, nUa)
    colnames(degree.smooth) <- c("Lambda", "Log(Lambda)", "df",
                  "Number of parameters", "Mean param.", "Scale param.", "Spline param.")
    rownames(degree.smooth) <- "  "

#    degree.smooth <- data.frame(lambda.use, df, df2, nparam, 
#                                nUregres-nestScale, nestScale, nUa)
#    colnames(degree.smooth) <- c("lambda", "df", "df2",
#                  "Number of parameters", "Mean param.", "Scale param.", "Spline param.")
#    rownames(degree.smooth) <- "  "


## AIC
    aic <- fit$loglik - df


## indicators of estimated components
    estimated <- c(est.scale, est.c, common.logscale)
    names(estimated) <- c("Scale", "ccoef", "common.logscale")

    temp <- list(regres = regres,
                 spline = spline,
                 loglik = loglik,
                 aic = aic,
                 degree.smooth = degree.smooth,
                 var = var,
                 var2 = var2,
                 dCdD = dCdD,
                 iter = fit$iter,
                 estimated = estimated,
                 warning = warn.all,
                 fail = fail.num
                 )

#    if (debug > 0){
       temp$H <- H
       temp$I <- I
       temp$G <- G
       temp$U <- U
#    }

    return(temp)

}


### ==========================================================
### MP.pseudoinv: Moore-Penrose pseudoinverse of the matrix x
### ==========================================================
##  * using the eigen-values decomposition
##  * eigen values lower than or equal to toler are considered to be 0
##  * !!! x is assumed to be symmetric
MP.pseudoinv <- function(x, toler = 1e-7){
   #eigen.x <- La.eigen(x, symmetric = TRUE)  ## From some point, La.eigen is no longer part of R as it is no longer needed.
                                              ## It became a default for eigen.
   eigen.x <- eigen(x, symmetric = TRUE)  
   nonZeroEV <- abs(eigen.x$values) > toler

   llambda <- sum(nonZeroEV)
   if (!llambda) return(NA)
   lambda <- eigen.x$values[nonZeroEV]
   
   ev.x <- diag(lambda, nrow=llambda, ncol=llambda)
   evec.x <- eigen.x$vectors[, nonZeroEV]
   x.inv <- evec.x %*% diag(1/lambda, nrow=llambda, ncol=llambda) %*% t(evec.x)    ## Moore-Penrose pseudoinverse of x
   
   return(x.inv)
}

