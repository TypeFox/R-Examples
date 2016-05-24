tobitTfit <- function(YS, XS, YO, XO, start,
                      print.level=0,
                      maxMethod="Newton-Raphson",
                      index=NULL,
                      binaryOutcome=FALSE,
                      ...) {
### Tobit treatment models:
### The latent variable is:
### YS* = XS'g + u
### The observables are:
###      / 1  if  YS* > 0
### YS = \ 0  if  YS* <= 0
### YO = X'b + YS bT + v
### u, v are correlated
### 
### Arguments:
### 
###  YS        binary or logical vector, 0 (FALSE) and 1 (TRUE)
###  XS              -"-                selection, should include
###              exclusion restriction
###  YO        numeric vector, outcomes
###  XO        explanatory variables for outcomes, should include YS
###  index     individual parameter indices in the parameter vector.
###            Should always be supplied but can generate here for
###            testing purposes
###  ...       additional parameters for maxLik
###
   loglik <- function( beta) {
      betaS <- beta[iBetaS]
      betaO <- beta[iBetaO]
      sigma <- beta[iSigma]
      if(sigma <= 0) return(NA)
      rho <- beta[iRho]
      if( ( rho < -1) || ( rho > 1)) return(NA)
                           # check the range
      XS0.betaS <- XS0%*%betaS
                           # denoted by 'z' in the vignette
      XS1.betaS <- XS1%*%betaS
      v0 <- YO0 - XO0%*%betaO
      v1 <- YO1 - XO1%*%betaO
      sqrt1r2 <- sqrt( 1 - rho^2)
      B0 <- (-XS0.betaS - rho/sigma*v0)/sqrt1r2
      B1 <- (XS1.betaS + rho/sigma*v1)/sqrt1r2
      loglik <- numeric(nObs)
      loglik[i0] <- -1/2*log( 2*pi) - log( sigma) -
          0.5*( v0/sigma)^2 + pnorm( B0, log.p=TRUE) 
      loglik[i1] <- -1/2*log( 2*pi) -log( sigma) -
          0.5*( v1/sigma)^2 + pnorm( B1, log.p=TRUE) 
      #sum(loglik)
      loglik
   }
   gradlik <- function(beta) {
      ## gradient is nObs x nParam matrix
      betaS <- beta[iBetaS]
      betaO <- beta[iBetaO]
      sigma <- beta[iSigma]
      if(sigma <= 0) return(NA)
      rho <- beta[iRho]
      if( ( rho < -1) || ( rho > 1)) return(NA)
                           # check the range
      XS0.betaS <- XS0%*%betaS
                           # denoted by 'z' in the vignette
      XS1.betaS <- XS1%*%betaS
      v0 <- drop(YO0 - XO0%*%betaO)
      v1 <- drop(YO1 - XO1%*%betaO)
      sqrt1r2 <- sqrt( 1 - rho^2)
      B0 <- (-XS0.betaS - rho/sigma*v0)/sqrt1r2
      B1 <- (XS1.betaS + rho/sigma*v1)/sqrt1r2
      lambda0 <- drop(lambda(B0))
      lambda1 <- drop(lambda(B1))
      ## now the gradient itself
      gradient <- matrix(0, nObs, nParam)
      gradient[i0, iBetaS] <- -lambda0*XS0/sqrt1r2
      gradient[i1, iBetaS] <- lambda1*XS1/sqrt1r2
      gradient[i0,iBetaO] <- (lambda0*rho/sigma/sqrt1r2
                              + v0/sigma^2)*XO0
      gradient[i1,iBetaO] <- (-lambda1*rho/sigma/sqrt1r2
                              + v1/sigma^2)*XO1
      gradient[i0,iSigma] <- (-1/sigma + v0^2/sigma^3
                              + lambda0*rho/sigma^2*v0/sqrt1r2)
      gradient[i1,iSigma] <- (-1/sigma + v1^2/sigma^3
                              - lambda1*rho/sigma^2*v1/sqrt1r2)
      gradient[i0,iRho] <- -lambda0*(v0/sigma + rho*XS0.betaS)/
          sqrt1r2^3
      gradient[i1,iRho] <- lambda1*(v1/sigma + rho*XS1.betaS)/
          sqrt1r2^3
#      colSums(gradient)
      gradient
   }
   hesslik <- function(beta) {
                           # This is a hack in order to avoid numeric problems
      ## gradient is nObs x nParam matrix
      betaS <- beta[iBetaS]
      betaO <- beta[iBetaO]
      sigma <- beta[iSigma]
      if(sigma <= 0) return(NA)
      rho <- beta[iRho]
      if( ( rho < -1) || ( rho > 1)) return(NA)
                           # check the range
      XS0.betaS <- XS0%*%betaS
                           # denoted by 'z' in the vignette
      XS1.betaS <- XS1%*%betaS
      v0 <- drop(YO0 - XO0%*%betaO)
      v1 <- drop(YO1 - XO1%*%betaO)
      sqrt1r2 <- sqrt( 1 - rho^2)
      B0 <- (-XS0.betaS - rho/sigma*v0)/sqrt1r2
      B1 <- (XS1.betaS + rho/sigma*v1)/sqrt1r2
      lambda0 <- drop(lambda(B0))
      lambda1 <- drop(lambda(B1))
      CB0 <- drop(CB(B0))
      CB1 <- drop(CB(B1))
      hess <- array(0, c( nParam, nParam))
      hess[,] <- NA
      hess[iBetaS,iBetaS] <-
         t( XS0) %*% ( XS0 * CB0)/sqrt1r2^2 +
             t( XS1) %*% ( XS1 * CB1)/sqrt1r2^2
      hess[iBetaS,iBetaO]  <-
         - t( XS0) %*% ( XO0 * CB0)*rho/sqrt1r2^2/sigma -
             t( XS1) %*% ( XO1 * CB1)*rho/sqrt1r2^2/sigma
      hess[iBetaO,iBetaS] <- t(hess[iBetaS,iBetaO])
      hess[iBetaS,iSigma] <-
         -rho/sigma^2/sqrt1r2^2*t( XS0) %*% ( CB0*v0) -
             rho/sigma^2/sqrt1r2^2*t( XS1) %*% ( CB1*v1)
      hess[iSigma,iBetaS] <- t(hess[iBetaS,iSigma])
      hess[iBetaS,iRho] <- 
         (t(XS0) %*% (CB0*(v0/sigma + rho*XS0.betaS)/sqrt1r2^4
                      - lambda0*rho/sqrt1r2^3) 
          +t(XS1) %*% (CB1*(v1/sigma + rho*XS1.betaS)/sqrt1r2^4
                       + lambda1*rho/sqrt1r2^3)
          )
      hess[iRho,iBetaS] <- t(hess[iBetaS,iRho])
      ##
      hess[iBetaO,iBetaO] <- 
         t( XO0) %*% (XO0*((rho/sqrt1r2)^2*CB0 - 1))/sigma^2 +
             t( XO1) %*% (XO1*( (rho/sqrt1r2)^2 * CB1 - 1))/sigma^2
      hess[iBetaO,iSigma] <-
         (t( XO0) %*% (CB0*rho^2/sigma^3*v0/sqrt1r2^2
                       - rho/sigma^2*lambda0/sqrt1r2 
                       - 2*v0/sigma^3) 
          + t( XO1) %*% (CB1*rho^2/sigma^3*v1/sqrt1r2^2 
                         + rho/sigma^2*lambda1/sqrt1r2
                         - 2*v1/sigma^3)
          )
      hess[iSigma,iBetaO] <- t(hess[iBetaO,iSigma])
      hess[iBetaO,iRho] <-
         (t(XO0) %*% (-CB0*(v0/sigma + rho*XS0.betaS)/sqrt1r2^4*rho
                      + lambda0/sqrt1r2^3)/sigma
          + t(XO1) %*% (-CB1*(v1/sigma + rho*XS1.betaS)/sqrt1r2^4*rho
                        - lambda1/sqrt1r2^3)/sigma
          )
      hess[iRho,iBetaO] <- t(hess[iBetaO,iRho])
      ##
      hess[iSigma,iSigma] <-
         (sum(1/sigma^2
             -3*v0*v0/sigma^4
             + v0*v0/sigma^4*rho^2/sqrt1r2^2*CB0
             -2*lambda0*v0/sqrt1r2*rho/sigma^3)
          + sum(1/sigma^2
                -3*v1*v1/sigma^4
                +rho^2/sigma^4*v1*v1/sqrt1r2^2*CB1
                +2*lambda1*v1/sqrt1r2*rho/sigma^3)
          )
      hess[iSigma,iRho] <- 
         (sum((-CB0*rho*(v0/sigma + rho*XS0.betaS)/sqrt1r2 + lambda0)
              *v0/sigma^2)/sqrt1r2^3
          - sum(
              (CB1*rho*(v1/sigma + rho*XS1.betaS)/sqrt1r2 + lambda1)
              *v1/sigma^2)/sqrt1r2^3
          )
      hess[iRho,iSigma] <- t(hess[iSigma,iRho])
      hess[iRho,iRho] <-
         (sum(CB0*( (v0/sigma + rho*XS0.betaS)/sqrt1r2^3)^2
              -lambda0*(XS0.betaS*(1 + 2*rho^2) + 3*rho*v0/sigma)/
                  sqrt1r2^5
              )
          + sum(CB1*( (v1/sigma + rho*XS1.betaS)/sqrt1r2^3)^2
                +lambda1*( XS1.betaS*( 1 + 2*rho^2) + 3*rho*v1/sigma) /
              sqrt1r2^5
                )
          )
      ## l.s2x3 is zero
      hess
   }
   ## ---------------
   NXS <- ncol( XS)
   if(is.null(colnames(XS)))
      colnames(XS) <- rep("XS", NXS)
   NXO <- ncol( XO)
   if(is.null(colnames(XO)))
      colnames(XO) <- rep("XO", NXO)
   nObs <- length( YS)
   i0 <- YS==0
   i1 <- YS==1
   NO1 <- length( YS[i0])
   NO2 <- length( YS[i1])
   ## indices in for the parameter vector
   if(is.null(index)) {
      iBetaS <- 1:NXS
      iBetaO <- max(iBetaS) + seq(length=NXO)
      if(!binaryOutcome) {
         iSigma <- max(iBetaO) + 1
         iRho <- max(iSigma) + 1
      }
      else
         iRho <- max(iBetaO) + 1
      nParam <- iRho
   }
   else {
      iBetaS <- index$betaS
      iBetaO <- index$betaO
      iSigma <- index$errTerms["sigma"]
      iRho <- index$errTerms["rho"]
      nParam <- index$nParam
   }
   ## split the data by selection
   XS0 <- XS[i0,,drop=FALSE]
   XS1 <- XS[i1,,drop=FALSE]
   YO0 <- YO[i0]
   YO1 <- YO[i1]
   XO0 <- XO[i0,,drop=FALSE]
   XO1 <- XO[i1,,drop=FALSE]
   ##
   if(print.level > 0) {
      cat( "Non-participants: ", NO1,
          "; participants: ", NO2, "\n", sep="")
      cat( "Initial values:\n")
      cat("selection equation betaS:\n")
      print(start[iBetaS])
      cat("Outcome equation betaO\n")
      print(start[iBetaO])
      cat("Variance sigma\n")
      print(start[iSigma])
      cat("Correlation rho\n")
      print(start[iRho])
   }
   result <- maxLik(loglik,
                    grad=gradlik,
                    hess=hesslik,
                    start=start,
                    print.level=print.level,
                    method=maxMethod,
                    ...)
   ## compareDerivatives(#loglik,
   ##     gradlik,
   ##     hesslik,
   ##                    t0=start)
   result$tobitType <- "treatment"
   result$method <- "ml"
   class( result ) <- c( "selection", class( result ) )
   return( result )
}
