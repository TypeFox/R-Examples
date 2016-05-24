probit <- function( formula, weights = NULL, ...) {
   ## Probit binary choice model.  Essentially a wrapper for "binaryCoice"
   ## formula: model formula, response must be either a logical or numeric vector containing only 0-s and
   ##          1-s
   ## start:      initial value of the parameters
   ## data     dataframe where the variables are defined
   ## x        whether to return model matrix
   ## y                                response vector
   ## model                            frame
   ## method   method for evaluation:
   ##          ML            maximum likelihood
   ##          model.frame   don't evaluate, only return the frame
   ## ...      further arguments for the maxLik algorithm
   ##
   ## return: a list with following components.
   ##  $results: maximisation results
   ##  $LRT:     list with two components:
   ##            LRT: likelihood ration test H0 - none of the variables significant
   ##            df:  corresponding degrees of freedom
   ##  x         model matrix (only if requested)
   ##  call      call
   ##  terms     terms
   ##
   loglik <- function( beta, weights ) {
      ## probit loglik, using Demidenko (2001) robust method
      x0 <- get("x0", envir=probitFrame + 1)
      x1 <- get("x1", envir=probitFrame + 1)
      Y <- get("Y", envir=probitFrame + 1)
      if( !is.null( weights ) ) {
         if( length( Y ) != length( weights ) ) {
            stop( "number of observations (", length( Y ),
               ") is not equal to number of weights (", length( weights ), ")")
         }
         w0 <- weights[ Y == 0 ]
         w1 <- weights[ Y == 1 ]
      } else {
         w0 <- 1
         w1 <- 1
      }
      xb0 <- drop(x0 %*% beta)
      xb1 <- drop(x1 %*% beta)
      loglik <- numeric(length(Y))
      loglik[Y == 0] <- w0 * pnorm(xb0, lower.tail=FALSE, log.p=TRUE)
      loglik[Y == 1] <- w1 * pnorm(xb1, lower.tail=TRUE, log.p=TRUE)
      ##
      f0 <- dnorm(xb0, log = TRUE)
      f1 <- dnorm(xb1, log = TRUE)
      F0 <- pnorm(xb0, lower.tail=FALSE, log.p = TRUE)
      F1 <- pnorm(xb1, lower.tail=TRUE, log.p = TRUE)
      gradlik <- matrix(0, length(Y), length(beta))
      theta3 <- exp( f1 - F1 )
      theta4 <- - exp( f0 - F0 )
      gradlik[Y == 1,] <- w1 * theta3*x1
      gradlik[Y == 0,] <- w0 * theta4*x0
      ##
      theta5 <- - xb1 * exp( f1 - F1 ) - exp(2*(f1 - F1 ))
      theta6 <- - exp(2*(f0 - F0) ) + xb0 * exp( f0 - F0 )
      hesslik <- t( x1) %*% ( w1 * x1 * theta5) + t( x0) %*% ( w0 * x0 * theta6)
                    # note that df/db' = -f (x'b) x'
      ## The following code can be used to compute Fisher Scoring approximation for the Hessian
      ## Note, it may be invertible in case of huge outliers where the Hessian is not.  However,
      ## the value of those standard errors are negligible anyway, so we do not use it by default.
      ## 
      ## theta7 <- -ifelse(xb1 < -nCutoff, -xb1*f1,
      ##                   ifelse(xb1 > nCutoff, xb1*f1,
      ##                          f1^2/F1/pnorm(xb1, lower.tail=FALSE)))
      ## theta8 <- -ifelse(xb0 < -nCutoff, -xb0*f0,
      ##                   ifelse(xb0 > nCutoff, xb0*f0,
      ##                          f0^2/F0/pnorm(xb0, lower.tail=TRUE)))
      ## FScore <- t( x1) %*% ( x1 * theta7) + t( x0) %*% ( x0 * theta8)
      attr(loglik, "gradient") <- gradlik
      attr(loglik, "hessian") <- hesslik
      loglik
   }
   # nCutoff <- 5
   probitFrame <- sys.nframe()
                           # the binaryChoice would create a few necessary variables in the next
                           # frame
   result <- binaryChoice(formula, ..., userLogLik = loglik, weights = weights )
   result$weights <- weights
   
   cl <- class(result)
   result <- c(result,
               family=list(binomial(link="probit"))
                           # NA action and the removed cases
               )
   class(result) <- c("probit", cl)
   result
}
