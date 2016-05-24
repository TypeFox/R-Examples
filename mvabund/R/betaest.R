# check if in line 178 an adjustment similar to the one at the beginning, taking into
# account if adjustment is necessay, needs to be done

################################################################################
## BETAEST Estimate the generalised linear model parameters, given nuisance   ##
## for a multivariate response                                                ##
################################################################################
### OLD: CHOICE: 1=Poisson, 2=negative binomial
##
## INPUT parameters:
## ABUNDANCES   = the response matrix
## X            = the matrix of the independent variables
## WEIGHTS      = vector of weights, weighted least squares is used  
##                (that is, minimizing sum(w*e^2))
## THETA        = the parameter of the neg bin model,
##                the model is Poisson if phi==0 
## TOLLEVEL     = the sensitivity in calculations near 0
## ITER         = maximum number of iterations
## MUSTART      = the initial value of mu
## LINKFUN, LINKINV, VARIANCE, MU.ETA, FAMILY.CHAR:
##              = functions depending on mu and model.parameters (phi)  
##                through the choice of the family and the link function
##                family(link)$linkfun, ... should be passed here
##                the DEFAULT values for these parameters are those  
##                appropriate for the poisson(link="log") family
##                For the VARIANCE of the neg. binomial family, 
##                use negative.binomial(phi=phi,link=...)$variance 
##                note that either phi or theta must be passed ! ? prob in multidims!
##
##                (using the log-link, it is better for computational efficiency to use
##                mu.eta <- function(eta){ pm <- exp(eta)
##                pm[pm<.Machine$double.eps] <- .Machine$double.eps
##                pm }
##                instead of the family$mu.eta function
##
## OUTPUT parameters: 
## beta, mu     = the model parameter estimates
## vrB          = a list of iteratively reweighted least squares
################################################################################

betaest <- function( abundances, x, weights=rep(1, times=nRows), 
    offset=matrix(0, nRows, nVars), phi, tollevel= 1.e-4, 
    iter=100, mustart=abundances+tollevel, family=poisson(), linkfun=family$linkfun,
    linkinv=family$linkinv, variance=family$variance, mu.eta=family$mu.eta,
    family.char=family$family[1], trace=FALSE, etastart=NULL, start=NULL ){

    abundances <- muOld <- mu <- tol <- as.matrix(abundances)  # eta  <-

    nRows         <- nrow(abundances)
    nVars         <- ncol(abundances)
    x             <- as.matrix(x)
    xnames        <- dimnames(x)[[2]]
    nParams       <- ncol(x)

    beta          <- matrix(0, nrow=nParams , ncol=nVars) 
    offset        <- as.matrix(offset)
    betaStart   <- NULL

   if (!is.null(etastart)){
       eta <- as.matrix(etastart)
       if(nrow(eta)!=nRows | ncol(eta)!=nParams)
          stop("'eta' should have the dimensions ", nRows, "x",nVars)
       mu <- linkinv(eta)
       if (!is.null(start)){
       if(is.vector(start)){
          betaStart <- matrix(rep(start, times=nVars),nParams, nVars)
       } else if(ncol(start)!=nVars){
          stop("columns of 'start' should equal", nVars,
            "and correspond to the response variables")
       } else betaStart <- start
       if (nrow(betaStart) != nParams){
            stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
              nParams, paste(deparse(xnames), collapse = ", ")),
              domain = NA)
        }
       }
    } else if (!is.null(start)){
       if(is.vector(start)){
          betaStart <- matrix(rep(start, times=nVars),nParams, nVars)
       } else if(ncol(start)!=nVars){
          stop("columns of 'start' should equal", nVars,
            "and correspond to the response variables")
       } else betaStart <- start

        if (nrow(betaStart) != nParams){
            stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
              nParams, paste(deparse(xnames), collapse = ", ")),
              domain = NA)
        } else {
            eta <-  offset + x %*% betaStart
        }
        mu <- linkinv(eta)
   } else {
     eta <- linkfun(mustart)
     mu  <- mustart
   }
   etaOld <- eta
   
    if(! family$validmu(mu)) {
        validmu.char <- unlist(strsplit(deparse(body(family$validmu)), " && "))
        if( any( "all(mu > 0)" == validmu.char  ))
            mu[mu<=0]     <- tollevel # adjust 0 counts
        if(any( "all(mu < 1)" == validmu.char ))    
            mu[mu >=1]    <- 1 - tollevel
    }

    if(regexpr("quasipoisson", family.char[1])!=-1){
        isPoiss <- family$phi == 1
    } else if(family.char[1]=="poisson") {
        isPoiss <- rep(TRUE, times = nVars) 
    } else if(regexpr("Negative Binomial", family.char[1])!=-1){
        isPoiss <- family$phi == 0
    } else isPoiss <- rep(FALSE, times = nVars) 

    tol[]           <- 1
    
    if(tollevel>=1) {
        stop("choose 'tollevel'<<1")
    } else converged      <- TRUE #  initialised converged

    step          <- 0
    warningParam   <- warningBeta <- c()
    Wi    <- vrB   <- list()
    # the Wi may have different lengths due to different number of 'good' data
    # in the first iteration of the loop, inProgress = 1:nVars and every instance of Wi
    # and vrB is initialised
    good <- weights > 0
    if (all(!good)) {
        stop("no observations informative")
    }
    
#####################
 
    if(family.char[1]=="gaussian" & family$link == "identity"){

        qrx <- qr(x[good,, drop = FALSE]*sqrt(weights[good]))
        vrB <- chol2inv( qrx$qr[1:qrx$rank,1:qrx$rank, drop=FALSE] )
        z   <- abundances + offset
        if( nParams > 0) {
        beta   <- vrB %*% t((x)[good,, drop = FALSE]) %*%
                                     z[good,, drop = FALSE]
        } else beta <- matrix(nrow=0, ncol=nVars)
        mu    <- (x %*% beta + offset)
        vrBl   <- list()
        for(i in 1:nVars)  vrBl[[i]]  <- vrB
        vrB   <- vrBl
        
        if (trace) {
          dev  <-  c(matrix(1, nrow=1,ncol=nRows) %*%
              family$dev.resids(abundances, mu, weights))
          cat("Deviance ="); print( dev); cat("\n")
        }
        step <- 1
#####################
  
    } else if (all(isPoiss) & family$link=="log") {

        # eta[]           <- linkfun(mu)
        muOld[]         <- 0
        etaOld[]        <- 0
        offsetplus1     <- offset + 1
        dev.resids      <- family$dev.resids
        dev             <- numeric(nVars)

        while (  max(c(tollevel-1.e-4, tol), na.rm=TRUE) > tollevel & step < iter )  {

            inProgress <- which(sapply(1:nVars, function(x) {
            (max(tol[,x])   > tollevel)  } ) )
             # find the variables ( max(tol,[],1) > tollevel )

            step      <- step + 1
            z         <- eta - offsetplus1 +  abundances/ mu
            wvrb      <-  mu*weights
            wvrb[wvrb > 1e+50] <- 1e+50

            if( nParams > 0) {
            for (i in inProgress) {

                txw <- t((x)[good,, drop = FALSE]* wvrb[good,i] )
                vrB[[i]]   <- try (chol2inv(chol( txw %*% x[good,, drop = FALSE])),
                  silent = TRUE)
                if (inherits(vrB[[i]] , "try-error") | any(is.na(vrB[[i]])) ) {
                    inProgress    <- inProgress[inProgress!=i]
                    warningParam  <- c(warningParam, i)
                    warningBeta   <- c(warningBeta, i)
                    converged     <- FALSE

                    vrB[[i]] <- matrix(NA, nrow=nParams, ncol=nParams)
                    beta[,i] <- mu[,i] <- NA
                } else beta[,i]  <- vrB[[i]] %*% txw %*% z[good,i, drop = FALSE]
            }
            } else {
                beta <- matrix(nrow=0, ncol=nVars)
            }
            etaOld[,inProgress] <- eta[,inProgress]
            eta                 <- x %*% beta + offset

            muOld[,inProgress] <- mu[,inProgress]
            mu[,inProgress]    <- (linkinv(eta[,inProgress]) )
            dev[inProgress]    <- c(matrix(1, nrow=1,ncol=nRows) %*%
                      dev.resids(abundances, mu, weights))[inProgress]
                      
            if (trace) {cat("Deviance ="); print( dev); cat("Iterations -", iter, "\n")}

            if(any(!is.finite(dev[inProgress]))) {
                if (is.null(betaStart))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                warning("step size truncated due to divergence", call. = FALSE)
                ii <- 1
                while(any(!is.finite(dev[inProgress]))) {
                  which.na <- (inProgress)[which(!is.finite(dev[inProgress]))]
                  ii <- ii + 1
                  beta[,which.na] <- (beta[,which.na] + betaStart[,which.na])/2
                  eta[,which.na] <- x %*% beta[,which.na] + offset[,which.na]
                  # (x %*% start)
                  mu[,which.na] <- linkinv(eta[,which.na])
                    # eval(linkinv.eta)[which(!is.finite(dev[inProgress]))]
                    # linkinv(eta <- eta + offset)
                  dev[which.na] <- c(matrix(1, nrow=1,ncol=nRows) %*%
                      dev.resids(abundances, mu, weights))[which.na]
                  if(ii > iter & any(!is.finite(dev[inProgress]))) {
                    which.na <- (inProgress)[which(!is.finite(dev[inProgress]))]
                    warning("inner loop 1; cannot correct step size for variable ",
                      which.na)
                    beta[,which.na] <- eta[,which.na] <- mu[,which.na] <- NA
                    dev[which.na] <- tol[,which.na] <- NA
                    inProgress <- inProgress[! inProgress %in% which.na]
                  }
                }
             if (trace) {cat("Step halved: new deviance ="); print(dev); cat( "\n")}
            }
            
            tol       <- abs ( mu - muOld )

            if(any(is.na(tol)) | any((abundances-mu)> (1e+5*max(abundances))) |
                any(!is.finite(mu)) ) {
                    which.na <- which(sapply(1:nVars, function(x) {
                    any(is.na(tol[,x])) | any((abundances-mu)[,x] > (1e+5*max(abundances[,x]))) |
                        any(!is.finite(mu)  ) } ) )
                    if(step>1){
                      mu[,which.na]  <- muOld[,which.na]
                      eta[,which.na] <- etaOld[,which.na]
                      beta[,which.na] <- betaStart[,which.na]
                    }
                    # exclude which.na values from further inProgress vectors
                    tol[,which.na]   <- NA
                    warningParam  <- unique(c(warningParam, which.na))
                    warningBeta   <- unique(c(warningBeta, which.na) )
                    inProgress <- inProgress[! inProgress %in% which.na]
            }
            betaStart   <- beta

        }
        if (step==iter) {
              warningParam <- c(warningParam, inProgress)
              converged <- FALSE
        }

 #####################  all other families and links
 
    } else {

         if(family.char== "Negative Binomial" & family$link == "varstab"){
              linkinv.eta <- expression(linkinv(eta)[,inProgress] )
         } else linkinv.eta <- expression(linkinv(eta[,inProgress]) )

         # eta[]           <- linkfun(mu)
         muOld[]         <- 0
         etaOld[]        <- 0
         mat0            <- muOld
         eye.nRows       <- matrix(0, nrow=nRows, ncol=nRows)
         eye.nRows[1+0:(nRows-1)*(nRows+1)]    <- 1
         dev.resids      <- family$dev.resids
         dev             <- numeric(nVars)

        varmu <- mu.eta.val <- mat0

        while (  max(c(tollevel-1.e-4, tol), na.rm=TRUE) > tollevel & step < iter )  {
         #  to avoid warning if all tol are NA

            inProgress <- which( sapply(1:nVars, function(x) {
            max(tol[,x]) > tollevel } ) )
            varmu[]       <-  variance(mu)
            mu.eta.val[]  <-  mu.eta(eta)
            good2 <- (mu.eta.val != 0)

            # if( family.char== "Negative Binomial"){
            #      dev.resids <- negative.binomial(phi=phi[inProgress],
            #         link=family$link)$dev.resids
            # } else if(family.char== "quasipoisson" )
            #      dev.resids <- quasipoisson(phi=phi[inProgress],
            #         link=family$link)$dev.resids
                  
            step      <- step + 1

            z         <- (eta - offset) + ( abundances - mu ) / mu.eta.val
            
            # var     <- variance(mu) 
            # for the neg bin family that cannot be calculated for one variable only 
              # in the multidim case because theta and phi are vectors

            wvrb <- weights * (mu.eta.val)^2/varmu
            wvrb[wvrb > 1e+50] <- 1e+50
            # is generally necessary
            if(any(is.na(varmu[good,inProgress, drop=FALSE]))){
                tmp <- is.na(colSums(varmu[good,inProgress, drop=FALSE]))
                warning("NAs in V(mu) in variable ", inProgress[tmp], " in step ", step)
                warningBeta   <- c(warningBeta, inProgress[tmp])
                inProgress    <- inProgress[- which(tmp)]
                converged     <- FALSE
            }
            if(any(varmu[good,inProgress, drop=FALSE] == 0)) {
                tmp <- (colSums(varmu[good,inProgress, drop=FALSE] == 0)>0)
                warning("0s in V(mu) in variable ", inProgress[tmp], " in step ", step )
                warningBeta  <- c(warningBeta, inProgress[tmp])
                inProgress   <- inProgress[- which(tmp)]
                converged    <- FALSE
            }
            if(any(is.na(mu.eta.val[good,inProgress, drop=FALSE]))){
                tmp <- is.na(colSums(mu.eta.val[good,inProgress, drop=FALSE]))
                warning("NAs in d(mu)/d(eta) in variable ", inProgress[tmp], " in step ", step)
                warningBeta  <- c(warningBeta, inProgress[tmp])
                inProgress   <- inProgress[- which(tmp)]
                converged    <- FALSE
            }

            # betaOld <- beta
            if( nParams > 0) {
            for (i in inProgress){
                    
                goodi <- good & good2[,i]
                if (all(! goodi)) {  
                    converged <- FALSE
                    warning(paste("no observations informative at iteration", 
                      step, "for y variable", i) )
                    tol[,i]  <- NA  # so that i will be excluded from further inProgress vectors
                    # for i the values of beta and subsequently mu, eta 
                      # remain as in the last iteration
                    next    # go to the next i in inProgress
                }  

                txw <- t(x[goodi,, drop = FALSE] * wvrb[goodi,i])
                vrB[[i]]   <- try (chol2inv(chol(txw  %*%
                               x[goodi,, drop = FALSE] )), silent = TRUE)
                muOld[,inProgress] <- mu[,inProgress]

                if (inherits(vrB[[i]] , "try-error") | any(is.na(vrB[[i]])) )  {
                    inProgress    <- inProgress[inProgress!=i]
                    warningParam  <- c(warningParam, i)
                    converged     <- FALSE
                    vrB[[i]] <- matrix(NA, nrow=nParams, ncol=nParams)
                    beta[,i] <- mu[,i] <- NA
                } else {
                    beta[,i]   <- vrB[[i]] %*% txw %*% z[goodi,i, drop = FALSE]
                }
            }
            } else beta <- matrix(nrow=0, ncol=nVars)
            etaOld[,inProgress] <- eta[,inProgress]
            eta         <- x %*% beta + offset
            # only change for inProgress so that the other values stay the same
                # --> important for stopping calc when iterative estimation goes wrong  
            mu[,inProgress]    <- eval(linkinv.eta)

            if ( inherits(family, "family.mvabund") |  regexpr("poisson",family.char)==1){
              dev[inProgress] <- c(matrix(1, nrow=1,ncol=nRows) %*%
                  dev.resids(abundances, mu, weights))[inProgress]
            } else {
            for (i in inProgress)
                dev[i] <- c(sum(dev.resids(abundances[,i], mu[,i], weights)))
            }
            if (trace) {cat("Deviance ="); print( dev); cat("Iterations -", iter, "\n")}

            if(any(!is.finite(dev[inProgress]))) {
                if (is.null(betaStart))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                warning("step size truncated due to divergence", call. = FALSE)
                ii <- 1
                while(any(!is.finite(dev[inProgress]))) {
                  which.na <- (inProgress)[which(!is.finite(dev[inProgress]))]
                  ii <- ii + 1
                  beta[,which.na] <- (beta[,which.na] + betaStart[,which.na])/2
                  eta[,which.na] <- x %*% beta[,which.na] + offset[,which.na]
                  # (x %*% start)
                  mu[,which.na] <- eval(linkinv.eta)[which(!is.finite(dev[inProgress]))]

                  if ( inherits(family, "family.mvabund") |
                    regexpr("poisson",family.char)==1){
                    dev[which.na] <- c(matrix(1, nrow=1,ncol=nRows) %*%
                      dev.resids(abundances, mu, weights))[which.na]
                  } else {
                    for (i in which.na)
                      dev[i] <- c(sum(dev.resids(abundances[,i], mu[,i], weights)))
                  }
                  if(ii > iter & any(!is.finite(dev[inProgress]))) {
                    which.na <- (inProgress)[which(!is.finite(dev[inProgress]))]
                    warning("inner loop 1; cannot correct step size for variable ",
                      which.na)
                    beta[,which.na] <- eta[,which.na] <- mu[,which.na] <- NA
                    dev[which.na] <- tol[,which.na] <- NA
                    inProgress <- inProgress[! inProgress %in% which.na]
                  }
                }
             if (trace) {cat("Step halved: new deviance ="); print(dev); cat( "\n")}
           }
           tol  <- abs ( mu - muOld )

            if(any(is.na(tol)) |
                any(!is.finite(mu)) ) {
                    which.na <- which(sapply(1:nVars, function(x) {
                    any(is.na(tol[,x])) | any((abundances-mu)[,x] > (1e+5*max(abundances[,x]))) |
                        any(!is.finite(mu[,x])) } ) )
                    if(step>1){
                      mu[,which.na]  <- muOld[,which.na]
                      eta[,which.na] <- etaOld[,which.na]
                      beta[,which.na] <- betaStart[,which.na]
                    }
                    # exclude which.na values from further inProgress vectors
                    tol[,which.na]   <- NA
                    warningParam  <- unique(c(warningParam, which.na))
                    warningBeta   <- unique(c(warningBeta, which.na) )
                    inProgress <- inProgress[! inProgress %in% which.na]
            }
            betaStart <- beta

            }
            if (step==iter) {
                warningParam <- c(warningParam, inProgress)
                converged <- FALSE
            }
     }
        
        if (!converged) {
            if(family.char[1]=="quasipoisson" |
                substring(family.char[1], first=1, last =17)=="Negative Binomial") {
                # warning("not convergence in beta estimation reached for variable", paste(inProgress, collapse =", "))
               } else warning(paste("algorithm did not converge for variables",
                paste(warningParam, collapse =", ") ))
               # warning(paste("estimates for variables", paste(inProgress, collapse =", "), "did not converge in IRLS"))
        }   
        # if (boundary) warning("algorithm stopped at boundary value")
        eps <- 10 * .Machine$double.eps
        if (family.char[1] == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("fitted probabilities numerically 0 or 1 occurred")
        }
        if (family.char[1] == "poisson") {
            if (any(mu < eps)) 
                warning("fitted rates numerically 0 occurred")
        }   

    return(list(beta=beta, mu=mu, vrb=vrB, converged=converged, step = step,
        warningBeta=warningBeta))

}

