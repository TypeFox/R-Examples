tobit2fit <- function(YS, XS, YO, XO, start,
                      weights = NULL, print.level=0,
                      maxMethod="Newton-Raphson",
                      ...) {
### The model is as follows (Amemiya 1985):
### The latent variables are:
### YS* = XS'g + u1
### y2* = X'b2 + u2
### y3* = X'b3 + u3
### The observables are:
###      / 1  if  YS* > 0
### YS = \ 0  if  YS* <= 0
###      / y2*  if  YS = 0
### y2 = \ 0    if  YS = 1
###      / 0    if  YS = 0
### y3 = \ y3*  if  YS = 1
### 
###  YS      binary or logical vector, 0 (FALSE) corresponds to observable y2*,
###          1 (TRUE) to observable y3*
###  y2, y3  numeric vector, outcomes
###  X       explanatory variables for outcomes
###  XS              -"-                selection, should include exclusion restriction
###  ...     additional parameters for maxLik
### Result:
### Object of class 'tobit2', derived from 'maxLik'.
### Includes all the components of maxLik and additionally
### twoStep   Results for Heckman two-step estimator, used for initial values
### 
   loglik <- function( beta) {
      g <- beta[ibetaS]
      b <- beta[ibetaO]
      sigma <- beta[isigma]
      if(sigma < 0) return(NA)
      rho <- beta[irho]
      if( ( rho < -1) || ( rho > 1)) return(NA)
      XS.g <- XS %*% g
      XO.b <- XO %*% b
      u2 <- YO - XO.b
      r <- sqrt( 1 - rho^2)
      B <- (XS.g + rho/sigma*u2)/r
      ll <- w * ifelse(YS == 0,
                   (pnorm(-XS.g, log.p=TRUE)),
                   -1/2*log(2*pi) - log(sigma) +
                   (pnorm(B, log.p=TRUE) - 0.5*(u2/sigma)^2)
                   )
      ll
   }
    gradlik <- function(beta) {
       g <- beta[ibetaS]
       b <- beta[ibetaO]
       sigma <- beta[isigma]
       if(sigma < 0) return(matrix(NA, nObs, nParam))
       rho <- beta[irho]
       if( ( rho < -1) || ( rho > 1)) return(matrix(NA, nObs, nParam))
       XS0.g <- as.numeric(XS0 %*% g)
       XS1.g <- as.numeric(XS1 %*% g)
       XO1.b <- as.numeric(XO1 %*% b)
                                        #      u2 <- YO1 - XO1.b
       u2 <- YO1 - XO1.b
       r <- sqrt( 1 - rho^2)
                                        #      B <- (XS1.g + rho/sigma*u2)/r
       B <- (XS1.g + rho/sigma*u2)/r
       lambdaB <- exp( dnorm( B, log = TRUE ) - pnorm( B, log.p = TRUE ) )
       gradient <- matrix(0, nObs, nParam)
       gradient[YS == 0, ibetaS] <- - w0 * XS0 *
          exp( dnorm( -XS0.g, log = TRUE ) - pnorm( -XS0.g, log.p = TRUE ) )
       gradient[YS == 1, ibetaS] <- w1 * XS1 * lambdaB/r
       gradient[YS == 1, ibetaO] <- w1 * XO1 * (u2/sigma^2 - lambdaB*rho/sigma/r)
       gradient[YS == 1, isigma] <- w1 * ( (u2^2/sigma^3 - lambdaB*rho*u2/sigma^2/r) - 1/sigma )
       gradient[YS == 1, irho] <- w1 * (lambdaB*(u2/sigma + rho*XS1.g))/r^3
       gradient
    }
    hesslik <- function(beta) {
       g <- beta[ibetaS]
       b <- beta[ibetaO]
       sigma <- beta[isigma]
       if(sigma < 0) {
         return( matrix( NA, nrow = nParam, ncol = nParam ) )
       }
       rho <- beta[irho]
       if( ( rho < -1) || ( rho > 1)) {
         return( matrix( NA, nrow = nParam, ncol = nParam ) )
       }
       XS0.g <- as.vector(XS0 %*% g)
       XS1.g <- as.vector(XS1 %*% g)
       XO1.b <- as.vector(XO1 %*% b)
       u2 <- YO1 - XO1.b
       r <- sqrt( 1 - rho^2)
       B <- (XS1.g + rho/sigma*u2)/r
       lambdaB <- exp( dnorm( B, log = TRUE ) - pnorm( B, log.p = TRUE ) )
       C <- ifelse(B > -500,
                   -exp(dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))*B -
                   exp(2 * (dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))),
                   -1)
                                        # recommended by Dimitrios Rizopoulos, KULeuven
                                        # This is a hack in order to avoid numerical problems.  How to do
                                        # it better?  How to prove the limit value?
       hess <- matrix(0, nParam, nParam)
       a <- ifelse( XS0.g < 500,
          -exp(dnorm(-XS0.g, log=TRUE) - pnorm(-XS0.g, log.p=TRUE))*XS0.g +
          ( exp( dnorm(-XS0.g, log=TRUE) - pnorm(-XS0.g, log.p=TRUE)))^2, 1 )
       hess[ibetaS,ibetaS] <- -t(XS0) %*% (w0*XS0*a) + t(XS1) %*% (w1*XS1*C)/r^2
       hess[ibetaS,ibetaO] <- -t(XS1) %*% (w1*XO1*C)*rho/r^2/sigma
       hess[ibetaO,ibetaS] <- t(hess[ibetaS,ibetaO])
       hess[ibetaS,isigma] <- -rho/sigma^2/r^2*t(XS1) %*% (w1*C*u2)
       hess[isigma,ibetaS] <- t(hess[ibetaS,isigma])
       hess[ibetaS,irho] <- t(XS1) %*%
           (w1*(C*(u2/sigma + rho*XS1.g)/r^4 + lambdaB*rho/r^3))
       hess[irho,ibetaS] <- t(hess[ibetaS,irho])
       hess[ibetaO,ibetaO] <- t(XO1) %*%
           (w1*(XO1 * ((rho/r)^2*C - 1)))/sigma^2
       hess[ibetaO,isigma] <- t(XO1) %*%
           (w1*(C*rho^2/sigma^3*u2/r^2 +
            rho/sigma^2*lambdaB/r - 2*u2/sigma^3))
       hess[isigma,ibetaO] <- t(hess[ibetaO,isigma])
       hess[ibetaO,irho] <- t(XO1) %*%
           (w1*(-C*(u2/sigma + rho*XS1.g)/r^4*rho -
            lambdaB/r^3))/sigma
       hess[irho,ibetaO] <- t(hess[ibetaO,irho])
       hess[isigma,isigma] <- sum( w1 *
                                  ( -3*u2*u2/sigma^4
                                  +2*lambdaB* u2/r *rho/sigma^3
                                  +rho^2/sigma^4 *u2*u2/r^2 *C) ) +
                                      sum(w1) / sigma^2
       hess[isigma,irho] <- hess[irho,isigma] <-
           -sum(w1*(C*rho*(u2/sigma + rho*XS1.g)/r + lambdaB)*
                u2/sigma^2)/r^3
       hess[irho,irho] <-
           sum(w1*(C*((u2/sigma + rho*XS1.g)/r^3)^2 +
               lambdaB*(XS1.g*(1 + 2*rho^2) + 3*rho*u2/sigma) / r^5 ))
       return( hess )
    }
    ## ---------------
    NXS <- ncol( XS)
    if(is.null(colnames(XS)))
        colnames(XS) <- rep("XS", NXS)
    NXO <- ncol( XO)
    if(is.null(colnames(XO)))
        colnames(XO) <- rep("XO", NXO)
    Nparam <- NXS + NXO + 2
                                        # Total # of parameters
   nObs <- length( YS)
    NO <- length( YS[YS > 0])
   nParam <- NXS + NXO + 2
   ## parameter indices
   ibetaS <- 1:NXS
   ibetaO <- seq(tail(ibetaS, 1)+1, length=NXO)
   isigma <- tail(ibetaO, 1) + 1
   irho <- tail(isigma, 1) + 1
   ## output, if asked for it
   if( print.level > 0) {
      cat( "Yo observed:", NO, "times; not observed:", nObs - NO, "times\n")
      cat( "Initial values:\n")
      cat("Selection\n")
      print(start[ibetaS])
      cat("Outcome\n")
      print(start[ibetaO])
      cat("sigma1\n")
      print(start[isigma])
      cat("rho1\n")
      print(start[irho])
   }
   ## pre-calculate a few values:
   XS0 <- XS[YS==0,,drop=FALSE]
   XS1 <- XS[YS==1,,drop=FALSE]
   YO[is.na(YO)] <- 0
   YO1 <- YO[YS==1]
   XO1 <- XO[YS==1,,drop=FALSE]
   N0 <- sum(YS==0)
   N1 <- sum(YS==1)
   
   if( !is.null( weights ) ) {
      w  <- weights
      w0 <- weights[ YS == 0 ]
      w1 <- weights[ YS == 1 ]
   } else {
      w  <- rep( 1, N0 + N1 )
      w0 <- rep( 1, N0 )
      w1 <- rep( 1, N1 )
   }
   
   # browser()
   # compareDerivatives(loglik, gradlik, t0=start )
   # range(numericGradient(loglik, t0=start)-gradlik(start))
   # compareDerivatives( function(k)sum(loglik(k)), function(k)colSums(gradlik(k)), hesslik, t0=start)
   # hesslik( start ) - numericHessian( function(k)sum(loglik(k)), function(k)colSums(gradlik(k)), t0=start )
 
   ## estimate
   result <- maxLik(loglik, grad=gradlik, hess=hesslik,
                    start=start,
                    method=maxMethod,
                    print.level=print.level, ...)
   result$tobitType <- 2
   result$method <- "ml"
   class( result ) <- c( "selection", class( result ) )
   return( result )
}
