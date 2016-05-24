
### Binary choice models
binaryChoice <- function(formula, subset, na.action,
                   start=NULL,
                   data=sys.frame(sys.parent()),
                   x=FALSE, y=FALSE, model=FALSE,
                   method="ML",
                         userLogLik=NULL,
                         cdfLower, cdfUpper=function(x) 1 - cdfLower(x),
                         logCdfLower=NULL, logCdfUpper=NULL,
                         pdf,logPdf=NULL,
                         gradPdf,
                         maxMethod="Newton-Raphson",
                   ...) {
   ## formula: model formula, response must be either a logical or numeric vector containing only 0-s and
   ##          1-s
   ## start:      initial value of the parameters
   ## data     dataframe where the variables are defined
   ## x        whether to return model matrix
   ## y                                response vector
   ## model                            frame
   ## method   method for evaluation:
   ##          ML            maximum likelihood
   ## userLogLik   user-supplied loglik function.  If supplied, that one will be used.  Note you might
   ##          want the userLogLik to supply the gradient and Hessian as attributes (see maxNR).
   ##          If not supplied, cdfLower, pdf and gradPdf are used instead.
   ## model.frame   don't evaluate, only return the frame
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
  ##  na.action  na.action used
  ##  cdfLower   lower tail of the cdf of the disturbance terms
  ##  cdfUpper   upper
   ##  logCdfLower, logCdfUpper  corresponding log functions, if available.  Increase precision
   ##              in extreme tails
  ##  pdf        pdf
  ##  gradPdf    gradient of the pdf
   genericLoglik <- function(beta) {
      ## A generic linear index parametric binary choice log-likelihood:
      ## F: disturbance terms cdf
      ## f:                   pdf
      ## 
      ## l = sum(log(1 - F(x0'beta))) + sum(log(F(x1'beta)))
      ## dl/dbeta = sum(-f(x0'beta)/(1 - F(x0'beta)) * x0) + sum(f(x1'beta)/F(x1'beta) * x1)
      ## d2l/dbeta2 =
      ##      sum((-df(x0'beta)/dbeta*(1 - F(x0'beta)) - f^2(x0'beta))/(1 - F(x0'beta))^2 *x0'x0 +
      ##      sum((df(x1'beta)/dbeta*F(x1'beta) - f^2(x1'beta))/F(x1'beta)^2 *x1'x1
      ##
      ## We need to supply 4 specific functions:
      ## F          cdfLower    (pnorm for probit)
      ## 1 - F      cdfUpper    (pnorm(.., lower.tail=FALSE))
      ## f          pdf         (dnorm)
      ## df/dbeta   gradPdf     (-dnorm(x)*x)
      xb0 <- drop(x0 %*% beta)
      xb1 <- drop(x1 %*% beta)
      #
      F0 <- cdfUpper(xb0)
      F1 <- cdfLower(xb1)
      if(!is.null(logCdfUpper)) {
         logF0 <- logCdfUpper(xb0)
         logF1 <- logCdfLower(xb1)
      }
      else {
         logF0 <- log(F0)
         logF1 <- log(F1)
      }
      loglik <- numeric(length(Y))
      loglik[Y == 0] <- logF0
      loglik[Y == 1] <- logF1
      ##
      f0 <- pdf(xb0)
      f1 <- pdf(xb1)
      gradlik <- matrix(0, length(Y), length(beta))
      if(!is.null(logPdf)) {
         r0 <- -exp(logPdf(xb0) - logF0)
         r1 <- exp(logPdf(xb1) - logF1)
      }
      else {
         r0 <- -pdf(xb0)/F0
         r1 <- pdf(xb1)/F1
      }         
      gradlik[Y == 0,] <- r0*x0
      gradlik[Y == 1,] <- r1*x1
      ##
      gradf0 <- gradPdf(xb0)
      gradf1 <- gradPdf(xb1)
      hesslik <-
          t( x0) %*% ( x0 * ( -gradf0*F0 - f0*f0)/F0/F0) +
              t( x1) %*% ( x1 * ( gradf1*F1 - f1*f1)/F1/F1)
      attr(loglik, "gradient") <- gradlik
      attr(loglik, "hessian") <- hesslik
      loglik
   }
   ## Binary choice specific stuff
   cl <- match.call(call=sys.call(sys.parent()))
                           # documentation claims 'sys.call(sys.parent())' is the default
                           # value for 'call'.  But in this way it works -- oterwise we just see the
                           # values for the last (passed through) arguments, not their original values
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf$drop.unused.levels <- TRUE
   mf[[1]] <- as.name("model.frame")
   eval(data)
                                        # we need to eval() data here, otherwise the evaluation of the
                                        # model frame will be wrong if called from inside a function
                                        # inside a function (sorry, I can't understand it either :-(
   mf <- eval(mf, envir=parent.frame())
   if( method == "model.frame" ) {
      if( !is.null( list(...)$weights ) ) {
         mf[[ "(weights)" ]] <- list(...)$weights
      }
       return(mf)
   } else if (method != "ML")
       warning("method = ", method, " is not supported. Using \"ML\"")
   mt <- attr(mf, "terms")
   Y <- model.response( mf )
   YLevels <- levels( as.factor( Y ) )
   if( length( YLevels ) != 2 ) {
      stop( "the left hand side of the 'formula' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   Y <- as.integer(Y == YLevels[ 2 ])
                                        # selection will be kept as integer internally
   X <- model.matrix(mt, mf, contrasts)
   nParam <- ncol( X)
   nObs <- length( Y)
   N1 <- sum(Y == 1)
   N0 <- nObs - N1
   if(N0 == 0 | N1 == 0) {
      stop("No variance in the response variable")
   }
   x0 <- X[Y==0,,drop=FALSE]
   x1 <- X[Y==1,,drop=FALSE]
   if(is.null(start)) {
      start <- rep( 0, nParam)
   }
   if(is.null(names(start))) {
      names(start) <- dimnames(X)[[2]]
   }
   ## Use the following for testing analytic derivatives
   ## 
   ## compareDerivatives(#f=function(theta) sum(userLogLik(theta)),
   ##                    f=function(theta) colSums(extractGrad(theta, userLogLik)),
   ##                    grad=function(theta) extractHess(theta, userLogLik),
   ##                    t0=start)
   ## stop()
   ##
   ## Main estimation
   if(is.null(userLogLik)) {
      estimation <- maxLik(genericLoglik, start=start,
                           method=maxMethod, ...)
   }
   else {
      estimation <- maxLik(userLogLik, start=start,
                           method=maxMethod, ...)
   }
   ## compare.derivatives(gradlik, hesslik, t0=start)
                                        #
   ## Likelihood ratio test: H0 -- all the coefficients, except intercept
   ## are zeros.  
   ll.bar <- N0*log(N0) + N1*log(N1) - (N0 + N1)*log(N0 + N1)
                                        # log-likelihood of the H0
   LRT <- 2*(logLik(estimation) - ll.bar)
                                        # note: this assumes that the model includes the constant
   result <- c(estimation,
               LRT=list(list(LRT=LRT, df=nParam-1)),
                                        # there are df-1 constraints
               param=list(list(nParam=nParam,nObs=nObs, N1=N1, N0=N0,
                               levels=YLevels)),
                           # We estimate probability for Y == levels[2]
                           # as opposite of levels[1]
               df.residual = nObs - sum( activePar( estimation ) ),
               call=cl,
               terms=mt,
               x=switch(x, "1"=list(X), "0"=NULL),
               y=switch(y, "1"=list(Y), "0"=NULL),
               model=switch(model, "1"=list(mf), "0"=NULL),
               na.action=list(attr(mf, "na.action"))
                           # NA action and the removed cases
               )
   class(result) <- c("binaryChoice", class(estimation))
   result
}
