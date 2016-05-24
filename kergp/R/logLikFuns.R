##============================================================================
## This is the log-likelihood CONCENTRATED with respect to beta for the
## model
##
##   y = F beta + eta
##
## where eta is Norm(0, C) with variance C := [ k(x_i, x_j) ]_{i,j}
##
## This function IS NOT maximised as such but copied during the call to
## the estimation function.
##
## By using a name begining by ".", 'R CMD check' will not complain about
## missing doc for this function.
##
##============================================================================

.logLikFun0 <- function(par,
                        object,
                        y, X, F, # = NULL,
                        gradEnv, # = NULL,
                        compGrad, # = TRUE,
                        noise,
                        trace # = FALSE) {
                        ) {
  
    n <- length(y)
    lpar <- length(par)
    
    if (trace) print(par)
    
    if (noise){
        lparNN <- lpar - 1L
        ## 'if' added by Yves on 2013-08-26
        if ( any( par[-lpar] < coefLower(object) ) ||
            any( par[-lpar] > coefUpper(object) )  ||
            par[lpar] < 0 ) {
            return(NA) 
        }
        coef(object) <- par[-lpar] ## MMM METHOD HERE
        C <- covMat(object, X, compGrad = compGrad, checkNames = FALSE, index = 1L) 
        ## MMM METHOD HERE
        diag(C) <- diag(C) + par[lpar]
    } else {
        ## 'if' added by Yves on 2013-08-26
        if ( any( par < coefLower(object) ) ||
            any( par > coefUpper(object) ) ) {
            return(NA) 
        }
        lparNN <- lpar
        coef(object) <- par ## MMM METHOD HERE
        C <- covMat(object, X, compGrad = compGrad, checkNames = FALSE, index = 1L) 
        ## MMM METHOD HERE
    }
    ## cat("object changed\n")
    
    ## ## XXX augment diagonal???
    ##   dC <- mean(diag(C))
    ##   md <- mean(dC)
    ##   diag(C) <- md * (1 + sqrt(sqrt(.Machine$double.eps)))
    
    L <- try(chol(C))
    if (inherits(L, "try-error")) {
        if (trace) cat ("chol error\n")
        return(NA)
    }
    L <- t(L)
    
    ## yStar <- backsolve(L, y, upper.tri = FALSE)
    
    if (!is.null(F)) {
        pF <- NCOL(F)
        p1 <- pF + 1L
        FStarPlus <- forwardsolve(L, cbind(F, y))
        qrFStarPlus <- qr(FStarPlus)
        dStar <- qr.R(qrFStarPlus)[p1, p1]
        eStar  <- qr.Q(qrFStarPlus)[ , p1] * dStar
        sseStar <- dStar^2
    } else {
        eStar <- forwardsolve(L, y)
        sseStar <- crossprod(eStar)
    }
    
    logLik <-  -0.5 * ( n * log(2*pi) + 2 * sum(log(diag(L))) + sseStar )
    
    if (compGrad) {
        
        ##=====================================================================
        ## partial derivative with respect to parameters.
        ## We make use of the 'scores' method with a precomputed W matrix.
        ##=====================================================================
        logLik.derivative <- array(0, dim = c(lpar, 1L),
                                   dimnames = list(names(par), NULL))
        u <- backsolve(t(L), eStar)   ## u := L^(-T) %*% eStar = C^(-1) * error
        ## Compute the 'W' vector of weights here
        W <- u %*% t(u) - chol2inv(t(L))
        dW <- diag(W)
        diag(W) <- dW / 2
        W <- as.numeric(W[lower.tri(W, diag = TRUE)])
        
        logLik.derivative[1L:lparNN] <- scores(object, X, weights = W)
        
        ## noise derivation: the derivative is wrt sigma^2 := varNoise  
        if (noise){
            logLik.derivative[lpar] <- sum(dW) / 2
        }
        
        if (!is.null(gradEnv)) {
            assign("par", par, envir = gradEnv)
            assign("LLgrad", logLik.derivative, envir = gradEnv)
        }
        
    }
    
    logLik <- as.numeric(logLik)
    
    if (compGrad) {
        attr(logLik, "gradient") <- logLik.derivative
    }
    
    if (trace) {
        cat("logLik = ", logLik[1], "\n")
        if (compGrad) cat("*** grad ***", attr(logLik, "gradient"), "***\n")
    }
    
    return(logLik)
    
} ## end logLikFun

## ## ---------------------------------------------------------------------
## ##  OLD VERSION (WITHOUT SCORES), corresponding to: if(compGrad){...} ##
##    
##       logLik.derivative <- array(0, dim = c(lpar, 1L), dimnames = c(names(par), NULL))
##     
##       ## Auxiliary computation for the derivative
##       ## Will be changed in the future: vector of weights passed
##       ## to C code
##       
##       u <- backsolve(t(L), eStar)   ## u := L^(-T) %*% eStar = C^(-1) * error
##       Cinv <- chol2inv(t(L))        ## Invert C from given L
##       
##       Cinv.upper <- Cinv[upper.tri(Cinv)]
##       Cinv.diag <- diag(Cinv)
##       uu <- u %*% t(u)
##       uu.upper <- uu[upper.tri(uu)]
##       uu.diag <- diag(uu)
##       
##       ## TO DO: Compute the 'W' vector of weights here
##       
##       ##========================================================================
##       ## partial derivative with respect to parameters
##       ## The loop could be done at the C level for efficiency
##       ##========================================================================
##       
##       for (k in 1L:lparNN) {
##         
##         if (k == 1L) {
##           gradC.k <- attr(C, "gradient")
##         } else {
##           ## XXX to improve : not useful to recompute 'C' here
##           aux2 <- covMat(object, X, compGrad = compGrad, checkNames = FALSE, index = k)    
##           ## MMM METHOD HERE
##           gradC.k <- attr(aux2, "gradient")
##         }
##         
##         gradC.k.upper <- gradC.k[upper.tri(gradC.k)]
##         gradC.k.diag <- diag(gradC.k)
##         term1 <- - 2 * sum(uu.upper * gradC.k.upper) 
##         term1 <- term1 - sum(uu.diag * gradC.k.diag)
##         #### economic computation of - t(x)%*%gradC.k%*%x
##         term2 <- 2*sum(Cinv.upper * gradC.k.upper) 
##         term2 <- term2 + sum(Cinv.diag * gradC.k.diag)
##         ## economic computation of trace(Cinv%*%gradC.k)           
##         logLik.derivative[k] <- -0.5 * (term1 + term2)
##         
##       }
##   
##       #### noise derivation: the derivative is wrt sigma^2 := varNoise 
##       if (noise){
##         term1 <- -sum(uu.diag)
##         term2 <- sum(Cinv.diag)
##         logLik.derivative[lpar] <- -0.5 * (term1 + term2)
##       }

