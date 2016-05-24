tobit5fit <- function(YS, XS, YO1, XO1, YO2, XO2, start,
                      print.level=0,
                      maxMethod="Newton-Raphson",
                      ...) {
### Tobit-5 models are defined as follows (Amemiya 1985):
### The latent variables are:
### YS* = XS'g + u1
### YO1* = X'b1 + u2
### YO2* = X'b2 + u3
### The observables are:
###      / 1  if  YS* > 0
### YS = \ 0  if  YS* <= 0
###       / y2*  if  YS = 0
### YO1 = \ 0    if  YS = 1
###       / 0    if  YS = 0
### YO2 = \ y3*  if  YS = 1
###
### Arguments:
### 
###  YS        binary or logical vector, 0 (FALSE) corresponds to observable y2*,
###            1 (TRUE) to observable y3*
###  YO1, YO2  numeric vector, outcomes
###  X         explanatory variables for outcomes
###  XS              -"-                selection, should include exclusion restriction
###  ...        additional parameters for maxLik
### 
    loglik <- function( beta) {
       betaS <- beta[iBetaS]
        b1 <- beta[iBetaO1]
        sigma1 <- beta[iSigma1]
       if(sigma1 <= 0) return(NA)
        rho1 <- beta[iRho1]
        if( ( rho1 < -1) || ( rho1 > 1)) return(NA)
        b2 <- beta[iBetaO2]
        sigma2 <- beta[iSigma2]
       if(sigma2 <= 0) return(NA)
        rho2 <- beta[iRho2]
        if((rho2 < -1) || (rho2 > 1)) return(NA)
                                        # check the range
        XS0.betaS <- XS0%*%betaS
        XS1.betaS <- XS1%*%betaS
        u1 <- YO1 - XO1%*%b1
        u2 <- YO2 - XO2%*%b2 
        sqrt1r22 <- sqrt( 1 - rho1^2)
        sqrt1r32 <- sqrt( 1 - rho2^2)
        B1 <- -(XS0.betaS + rho1/sigma1*u1)/sqrt1r22
        B2 <- (XS1.betaS + rho2/sigma2*u2)/sqrt1r32
        loglik <- numeric(nObs)
        loglik[YS == 0] <- -log( sigma1) - 0.5*( u1/sigma1)^2 + pnorm( B1, log.p=TRUE) - 1/2*log( 2*pi)
       loglik[YS == 1] <- -log( sigma2) - 0.5*( u2/sigma2)^2 + pnorm( B2, log.p=TRUE) - 1/2*log( 2*pi)
       loglik
}
    gradlik <- function(beta) {
       ## gradient is nObs x nParam matrix
        ## components of the gradient are ordered as: g b_2, s_2, r_2, b_3,
        ## s_3, r_3
        betaS <- beta[iBetaS]
        b1 <- beta[iBetaO1]
        sigma1 <- beta[iSigma1]
       if(sigma1 <= 0) return(NA)
        rho1 <- beta[iRho1]
        if( ( rho1 < -1) || ( rho1 > 1)) return(NA)
        b2 <- beta[iBetaO2]
        sigma2 <- beta[iSigma2]
       if(sigma2 <= 0) return(NA)
        rho2 <- beta[iRho2]
        if((rho2 < -1) || (rho2 > 1)) return(NA)
                                        # check the range
        XS0.betaS <- XS0%*%betaS
        XS1.betaS <- XS1%*%betaS
        XO1.b <- XO1%*%b1
        XO2.b <- XO2%*%b2
        u1 <- as.vector(YO1 - XO1.b)
        u2 <- as.vector(YO2 - XO2.b)
        sqrt1r22 <- sqrt( 1 - rho1^2)
        sqrt1r32 <- sqrt( 1 - rho2^2)
        B1 <- -(XS0.betaS + rho1/sigma1*u1)/sqrt1r22
        B2 <- (XS1.betaS + rho2/sigma2*u2)/sqrt1r32
        lambda1 <- as.vector(ifelse(B1 > -30, dnorm(B1)/pnorm(B1), -B1))
        lambda2 <- as.vector(ifelse(B2 > -30, dnorm(B2)/pnorm(B2), -B2))
                                        # This is a hack in order to avoid numeric problems
        ## now the gradient itself
        gradient <- matrix(0, nObs, nParam)
        gradient[YS == 0, iBetaS] <- -lambda1*XS0/sqrt1r22
        gradient[YS == 1, iBetaS] <- lambda2*XS1/sqrt1r32
        gradient[YS == 0,iBetaO1] <- (lambda1*rho1/sigma1/sqrt1r22 + u1/sigma1^2)*XO1
        gradient[YS == 1,iBetaO2] <- (-lambda2*rho2/sigma2/sqrt1r32 + u2/sigma2^2)*XO2
        gradient[YS == 0,iSigma1] <- (-1/sigma1 + u1^2/sigma1^3
                 +lambda1*rho1/sigma1^2*u1/sqrt1r22)
        gradient[YS == 1,iSigma2] <- (-1/sigma2 + u2^2/sigma2^3
                 -lambda2*rho2/sigma2^2*u2/sqrt1r32)
        gradient[YS == 0,iRho1] <- -lambda1*( u1/sigma1 + rho1*XS0.betaS)/sqrt1r22^3
        gradient[YS == 1,iRho2] <- lambda2*( u2/sigma2 + rho2*XS1.betaS)/sqrt1r32^3
        gradient
    }
    hesslik <- function(beta) {
        betaS <- beta[iBetaS]
        b1 <- beta[iBetaO1]
        sigma1 <- beta[iSigma1]
        if(sigma1 <= 0) {
            return( matrix( NA, nrow = nParam, ncol = nParam ) )
        }
        rho1 <- beta[iRho1]
        if( abs( rho1 ) > 1 ) {
            return( matrix( NA, nrow = nParam, ncol = nParam ) )
        }
        b2 <- beta[iBetaO2]
        sigma2 <- beta[iSigma2]
        if(sigma2 <= 0) {
            return( matrix( NA, nrow = nParam, ncol = nParam ) )
        }
        rho2 <- beta[iRho2]
        if( abs( rho2 ) > 1 ) {
            return( matrix( NA, nrow = nParam, ncol = nParam ) )
        }
        XS0.bS <- XS0%*%betaS
        XS1.bS <- XS1%*%betaS
        u1 <- YO1 - XO1%*%b1
        u2 <- YO2 - XO2%*%b2
        sqrt1r22 <- sqrt( 1 - rho1^2)
        sqrt1r32 <- sqrt( 1 - rho2^2)
        B1 <- -(XS0.bS + rho1/sigma1*u1)/sqrt1r22
        B2 <- (XS1.bS + rho2/sigma2*u2)/sqrt1r32
        lambda1 <- ifelse(B1 > -30, dnorm(B1)/pnorm(B1), -B1)
        lambda2 <- ifelse(B2 > -30, dnorm(B2)/pnorm(B2), -B2)
                                        # This is a hack in order to avoid numeric problems
        CB1 <- as.vector(ifelse(B1 > -500,
                                -exp(dnorm(B1, log = TRUE) - pnorm(B1, log.p = TRUE))*B1 -
                                exp(2 * (dnorm(B1, log = TRUE) - pnorm(B1, log.p = TRUE))),
                                -1))
        CB2 <- as.vector(ifelse(B2 > -500,
                                -exp(dnorm(B2, log = TRUE) - pnorm(B2, log.p = TRUE))*B2 -
                                exp(2 * (dnorm(B2, log = TRUE) - pnorm(B2, log.p = TRUE))),
                                -1))
                                        # recommended by Dimitrios Rizopoulos, KULeuven
                                        # This is a hack in order to avoid numerical problems.  How to do
                                        # it better?  How to prove the limit value?
        l.gg <- t( XS0) %*% ( XS0 * CB1)/sqrt1r22^2 +
            t( XS1) %*% ( XS1 * CB2)/sqrt1r32^2
        l.gb1 <- -t( XS0) %*%
            ( XO1 * CB1)*rho1/sqrt1r22^2/sigma1
        l.gs2 <- -rho1/sigma1^2/sqrt1r22^2*
            t( XS0) %*% ( CB1*u1)
        l.gr2 <- t( XS0) %*%
            ( CB1*( u1/sigma1 + rho1*XS0.bS)/sqrt1r22^4 -
             lambda1*rho1/sqrt1r22^3)
        l.gb2 <- -t( XS1) %*%
            ( XO2 * CB2)*rho2/sqrt1r32^2/sigma2
        l.gs3 <- -rho2/sigma2^2/sqrt1r32^2*
            t( XS1) %*% ( CB2*u2)
        l.gr3 <- t( XS1) %*%
            ( CB2*( u2/sigma2 + rho2*XS1.bS)/sqrt1r32^4 +
             lambda2*rho2/sqrt1r32^3)
        l.b1b1 <- t( XO1) %*%
            (XO1 * ( (rho1/sqrt1r22)^2 * CB1 - 1))/sigma1^2
        l.b1s2 <- t( XO1) %*%
            ( CB1*rho1^2/sigma1^3*u1/sqrt1r22^2 -
             rho1/sigma1^2*lambda1/sqrt1r22 -
             2*u1/sigma1^3)
        l.b1r2 <- t( XO1) %*%
            ( -CB1*( u1/sigma1 + rho1*XS0.bS)/sqrt1r22^4*rho1 +
             lambda1/sqrt1r22^3)/sigma1
        ## l.b1x3 is zero
        l.s2s2 <- sum(
                      1/sigma1^2
                      -3*u1*u1/sigma1^4
                      + u1*u1/sigma1^4 *rho1^2/sqrt1r22^2 *CB1
                      -2*lambda1* u1/sqrt1r22 *rho1/sigma1^3)
        l.s2r2 <- sum(
                      ( -CB1*rho1*(u1/sigma1 + rho1*XS0.bS)/sqrt1r22 +
                       lambda1)
                      *u1/sigma1^2)/sqrt1r22^3
        ## l.s2x3 is zero
        l.r2r2 <- sum(
                      CB1*( ( u1/sigma1 + rho1*XS0.bS)/sqrt1r22^3)^2
                      -lambda1*( XS0.bS*( 1 + 2*rho1^2) + 3*rho1*u1/sigma1) /
                      sqrt1r22^5
                      )
        ## l.r2x3 is zero
        l.b2b2 <- t( XO2) %*%
            (XO2 * ( (rho2/sqrt1r32)^2 * CB2 - 1))/sigma2^2
        l.b2s3 <- t( XO2) %*%
            ( CB2*rho2^2/sigma2^3*u2/sqrt1r32^2 +
             rho2/sigma2^2*lambda2/sqrt1r32 - 2*u2/sigma2^3)
        l.b2r3 <- t( XO2) %*%
            ( -CB2*( u2/sigma2 + rho2*XS1.bS)/sqrt1r32^4*rho2 -
             lambda2/sqrt1r32^3)/sigma2
        l.s3s3 <- sum(
                      1/sigma2^2
                      -3*u2*u2/sigma2^4
                      +2*lambda2* u2/sqrt1r32 *rho2/sigma2^3
                      +rho2^2/sigma2^4 *u2*u2/sqrt1r32^2 *CB2)
        l.s3r3 <- -sum(
                       ( CB2*rho2*(u2/sigma2 + rho2*XS1.bS)/sqrt1r32 +
                        lambda2)
                       *u2/sigma2^2)/sqrt1r32^3
        l.r3r3 <- sum(
                      CB2*( ( u2/sigma2 + rho2*XS1.bS)/sqrt1r32^3)^2
                      + lambda2*( XS1.bS*( 1 + 2*rho2^2) + 3*rho2*u2/sigma2) /
                      sqrt1r32^5
                      )
        hess <- array(NA, c( nParam, nParam))
        hess[iBetaS,iBetaS] <- l.gg
        hess[iBetaS,iBetaO1] <- l.gb1; hess[iBetaO1,iBetaS] <- t( l.gb1)
        hess[iBetaS,iSigma1] <- l.gs2; hess[iSigma1,iBetaS] <- t( l.gs2)
        hess[iBetaS,iRho1] <- l.gr2; hess[iRho1,iBetaS] <- t( l.gr2)
        hess[iBetaS,iBetaO2] <- l.gb2; hess[iBetaO2,iBetaS] <- t( l.gb2)
        hess[iBetaS,iSigma2] <- l.gs3; hess[iSigma2,iBetaS] <- t( l.gs3)
        hess[iBetaS,iRho2] <- l.gr3; hess[iRho2,iBetaS] <- t( l.gr3)
        hess[iBetaO1,iBetaO1] <- l.b1b1
        hess[iBetaO1,iSigma1] <- l.b1s2; hess[iSigma1,iBetaO1] <- t( l.b1s2)
        hess[iBetaO1,iRho1] <- l.b1r2; hess[iRho1,iBetaO1] <- t( l.b1r2)
        hess[iBetaO1,iBetaO2] <- 0; hess[iBetaO2,iBetaO1] <- 0
        hess[iBetaO1,iSigma2] <- 0; hess[iSigma2,iBetaO1] <- 0
        hess[iBetaO1,iRho2] <- 0; hess[iRho2,iBetaO1] <- 0
        hess[iSigma1,iSigma1] <- l.s2s2
        hess[iSigma1,iRho1] <- l.s2r2; hess[iRho1,iSigma1] <- l.s2r2
        hess[iSigma1,iBetaO2] <- 0; hess[iBetaO2,iSigma1] <- 0
        hess[iSigma1,iSigma2] <- 0; hess[iSigma2,iSigma1] <- 0
        hess[iSigma1,iRho2] <- 0; hess[iRho2,iSigma1] <- 0
        hess[iRho1,iRho1] <- l.r2r2
        hess[iRho1,iBetaO2] <- 0; hess[iBetaO2,iRho1] <- 0
        hess[iRho1,iSigma2] <- 0; hess[iSigma2,iRho1] <- 0
        hess[iRho1,iRho2] <- 0; hess[iRho2,iRho1] <- 0
        hess[iBetaO2,iBetaO2] <- l.b2b2
        hess[iBetaO2,iSigma2] <- l.b2s3; hess[iSigma2,iBetaO2] <- t( l.b2s3)
        hess[iBetaO2,iRho2] <- l.b2r3; hess[iRho2,iBetaO2] <- t( l.b2r3)
        hess[iSigma2,iSigma2] <- l.s3s3
        hess[iSigma2,iRho2] <- l.s3r3; hess[iRho2,iSigma2] <- t( l.s3r3)
        hess[iRho2,iRho2] <- l.r3r3
        hess
    }
    ## ---------------
    NXS <- ncol( XS)
    if(is.null(colnames(XS)))
        colnames(XS) <- rep("XS", NXS)
    NXO1 <- ncol( XO1)
    if(is.null(colnames(XO1)))
        colnames(XO1) <- rep("XO1", NXO1)
    NXO2 <- ncol( XO2)
    if(is.null(colnames(XO2)))
        colnames(XO2) <- rep("XO2", NXO2)
    nParam <- NXS + NXO1 + NXO2 + 4
                                        # Total # of parameters
    nObs <- length( YS)
    NO1 <- length( YS[YS==0])
    NO2 <- length( YS[YS==1])
    YO <- ifelse(YS==0, YO1, YO2)
    ## indices in for the parameter vector
    iBetaS <- 1:NXS
    iBetaO1 <- seq(tail(iBetaS, 1)+1, length=NXO1)
    iSigma1 <- tail(iBetaO1, 1) + 1
    iRho1 <- tail(iSigma1, 1) + 1
    iBetaO2 <- seq(tail(iRho1, 1) + 1, length=NXO2)
    iSigma2 <- tail(iBetaO2, 1) + 1
    iRho2 <- tail(iSigma2, 1) + 1
    ## split the data by selection
    XS0 <- XS[YS==0,,drop=FALSE]
    XS1 <- XS[YS==1,,drop=FALSE]
    YO1 <- YO[YS==0]
    YO2 <- YO[YS==1]
    XO1 <- XO1[YS==0,,drop=FALSE]
    XO2 <- XO2[YS==1,,drop=FALSE]
    ##
    if(print.level > 0) {
        cat( "Choice 2:", NO1, "times; choice 3:", NO2, "times\n")
    }
    if( print.level > 0) {
        cat( "Initial values:\n")
        cat("beta1\n")
        print(start[iBetaO1])
        cat("sigma1\n")
        print(start[iSigma1])
        cat("rho1\n")
        print(start[iRho1])
        cat("beta2\n")
        print(start[iBetaO2])
        cat("sigma2\n")
        print(start[iSigma2])
        cat("rho2\n")
        print(start[iRho2])
    }
    result <- maxLik(loglik, grad=gradlik,
                     hess=hesslik, start=start,
                     print.level=print.level,
                     method=maxMethod,
                     ...)
#   compareDerivatives(gradlik, hesslik, t0=start)
   result$tobitType <- 5
   result$method <- "ml"
   class( result ) <- c( "selection", class( result ) )
   return( result )
}
