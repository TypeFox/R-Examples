tobit2Bfit <- function(YS, XS, YO, XO, start,
                      weights = NULL, print.level=0,
                      maxMethod="BHHH",
                      ...) {
   ## Fit 2-dimensional sample selection models where the outcome
   ## variable is binary.
   ## The model is as follows (following Amemiya 1985):
   ## 
   ## The latent variables are:
   ## YS* = XS'g + u1
   ## YO* = X'b2 + u2
   ## 
   ## The observables are:
   ##      / 1  if  YS* > 0
   ## YS = \ 0  if  YS* <= 0
   ##      / 1(YO* > 0)  if  YS = 0
   ## YO = \ 0           if  YS = 1
   ## 
   ##  YS      binary or logical vector, 0 (FALSE) corresponds to
   ##          YO not observed, 1 (TRUE) if observed
   ##  XS      matrix of explanatory variables for selection equation,
   ##          should include exclusion restriction
   ##  YO      binary or logical outcome vector
   ##  XS      matrix of explanatory variables for outcomes
   ##  ...     additional parameters for maxLik
   ##  
   ## Result:
   ## Object of class 'tobit2', derived from 'maxLik'.
   ## Includes all the components of maxLik and additionally
   ## twoStep   Results for Heckman two-step estimator, used for initial values
   ## 
   loglik <- function( beta) {
      betaS <- beta[ibetaS]
      betaO <- beta[ibetaO]
      rho <- beta[iRho]
      if( ( rho < -1) || ( rho > 1)) return(NA)
      Sigma <- matrix(c(1,rho,rho,1), 2, 2)
      XS00.b <- drop(XS[i00,] %*% betaS)
      XS10.b <- drop(XS[i10,] %*% betaS)
      XS11.b <- drop(XS[i11,] %*% betaS)
      XO10.b <- drop(XO[i10,] %*% betaO)
      XO11.b <- drop(XO[i11,] %*% betaO)
      lik <- loglik <- numeric(nObs)
         # lik is the *unweighted* likelihood function needed for the gradient
      ## YS == 0, YO == 0
      loglik[i00] <- w00 * pnorm(-XS00.b, log.p=TRUE)
      lik[i00] <- exp(loglik[i00]/w00)
      ## YS == 1, YO == 0
      lowermat <- -cbind(XS10.b, XO10.b)
      f2 <- apply(lowermat, 1,
                  function(x) pmvnorm(lower=x, corr=Sigma))
      lik[i10] <- pnorm(-XS10.b, lower.tail=FALSE) - f2
      loglik[i10] <- w10 * log(lik[i10])
      ## YS == 1, YO == 1
      lowermat <- -cbind(XS11.b, XO11.b)
      f2 <- apply(lowermat, 1,
                  function(x) pmvnorm(lower=x, corr=Sigma))
      lik[i11] <- f2
      loglik[i11] <- w11 * log(f2)
      ## --- gradient ---
      ## The d loglik/d rho was taken from Henry Nyberg "A Bivariate Autoregressive Probit Model:
      ##   Predicting U.S. Business Cycle and Growth Rate Cycle Recessions",
      ##   Helsinki University working paper 272 (2009)
      grad <- matrix(0, nObs, nParam)
      r <- sqrt(1 - rho^2)
      ## YS == 0, YO == 0
      grad[i00,ibetaS] <- w00 * XS[i00,] * (-dnorm(-XS00.b)/lik[i00])
      ## YS == 1, YO == 0
      A <- dnorm(XS10.b)
      B <- A*pnorm((XO10.b - rho*XS10.b)/r, lower.tail=FALSE)
      grad[i10,ibetaS] <- w10 * XS[i10,]*B/lik[i10]
      A <- dnorm(XO10.b)
      B <- A*pnorm((XS10.b - rho*XO10.b)/r)
      grad[i10,ibetaO] <- -w10 * XO[i10,]*B/lik[i10]
      locmat <- -cbind(XS10.b, XO10.b)
      pdf <- apply(locmat, 1,
                   function(x) dmvnorm(x, c(0,0), Sigma))
      grad[i10,iRho] <- -w10 * pdf/lik[i10]
      ## YS == 1, YO == 1
      A <- dnorm(XS11.b)
      B <- A*pnorm((XO11.b - rho*XS11.b)/r)
      grad[i11,ibetaS] <- w11 * XS[i11,]*B/lik[i11]
      A <- dnorm(XO11.b)
      B <- A*pnorm((XS11.b - rho*XO11.b)/r)
      grad[i11,ibetaO] <- w11 * XO[i11,]*B/lik[i11]
      locmat <- -cbind(XS11.b, XO11.b)
      pdf <- apply(locmat, 1,
                   function(x) dmvnorm(x, c(0,0), Sigma))
      grad[i11,iRho] <- w11 * pdf/lik[i11]
      ## loglik <- sum(loglik)
      ## grad <- colSums(grad)
      attr(loglik, "gradient") <- grad
      return(loglik)
   }
   gradlik <- function(x) {
      l <- loglik(x)
      return(attr(l, "gradient"))
   }
    ## ---------------
    NXS <- ncol( XS)
    if(is.null(colnames(XS)))
        colnames(XS) <- rep("XS", NXS)
    NXO <- ncol( XO)
    if(is.null(colnames(XO)))
        colnames(XO) <- rep("XO", NXO)
   nObs <- length( YS)
    NO <- length( YS[YS > 0])
   ## selection indices
   i00 <- YS == 0
   i10 <- (YS == 1) & (YO == 0)
   i11 <- (YS == 1) & (YO == 1)
   ## parameter indices
   ibetaS <- 1:NXS
   ibetaO <- seq(tail(ibetaS, 1)+1, length=NXO)
   iRho <- tail(ibetaO, 1) + 1
   nParam <- iRho
   
   # weights
   if( !is.null( weights ) ) {
      w  <- weights
      w00 <- weights[ i00 ]
      w10 <- weights[ i10 ]
      w11 <- weights[ i11 ]
   } else {
      w  <- rep( 1, nObs )
      w00 <- rep( 1, sum( i00 ) )
      w10 <- rep( 1, sum( i10 ) )
      w11 <- rep( 1, sum( i11 ) )
   }
   
   ## output, if asked for it
   if( print.level > 0) {
      cat("YO observed:", NO, "times; not observed:", nObs - NO,
          "times:\n")
      print(table(YS, YO, exclude=NULL))
      cat( "Initial values:\n")
      print(start)
   }
   
   # browser()
   # compareDerivatives(loglik, gradlik, t0=start )
   # range(numericGradient(loglik, t0=start)-gradlik(start))
   
   ## estimate
   result <- maxLik(loglik, 
                    start=start,
                    method=maxMethod,
                    print.level=print.level, ...)
   result$tobitType <- 2
   result$method <- "ml"
   class( result ) <- c( "selection", class( result ) )
   return( result )
}
