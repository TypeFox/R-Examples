heckitTfit <- function(selection, outcome,
                       data=sys.frame(sys.parent()),
                       ys=FALSE, yo=FALSE,
                       xs=FALSE, xo=FALSE,
                       mfs=FALSE, mfo=FALSE,
                       print.level=0,
                       maxMethod="Newton-Raphson", ... ) {
   ## 2-step all-normal treatment effect estimator
   ## not public API
   ##
   ## maxMethod:   probit method
   ##
   ## Do a few sanity checks...
   if( class( selection ) != "formula" ) {
      stop( "argument 'selection' must be a formula" )
   }
   if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   }
   thisCall <- match.call()
   ## extract selection frame
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("selection", "data", "subset"), names(mf), 0)
   mfS <- mf[c(1, m)]
   mfS$drop.unused.levels <- TRUE
   mfS$na.action <- na.pass
   mfS[[1]] <- as.name("model.frame")
   names(mfS)[2] <- "formula"
                           # model.frame requires the parameter to
                           # be 'formula'
   mfS <- eval(mfS, parent.frame())
   mtS <- attr(mfS, "terms")
   XS <- model.matrix(mtS, mfS)
   YS <- model.response( mfS )
   YSLevels <- levels( as.factor( YS ) )
   if( length( YSLevels ) != 2 ) {
      stop( "the dependent variable of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   ysNames <- names( YS )
   YS <- as.integer(YS == YSLevels[ 2 ])
                           # selection kept as integer internally
   names( YS ) <- ysNames

   ## check for NA-s.  Because we have to find NA-s in several frames, we cannot use the standard 'na.'
   ## functions here.  Find bad rows and remove them later.
   badRow <- !complete.cases(YS, XS)
   badRow <- badRow | is.infinite(YS)
   badRow <- badRow | apply(XS, 1, function(v) any(is.infinite(v)))
   if("formula" %in% class( outcome)) {
      if( length( outcome ) != 3 ) {
         stop( "argument 'outcome1' must be a 2-sided formula" )
      }
      m <- match(c("outcome", "data", "subset"), names(mf), 0)
      mfO <- mf[c(1, m)]
      mfO$drop.unused.levels <- TRUE
      mfO$na.action <- na.pass
      mfO[[1]] <- as.name("model.frame")
      names(mfO)[2] <- "formula"
      mfO <- eval(mfO, parent.frame())
      mtO <- attr(mfO, "terms")
      XO <- model.matrix(mtO, mfO)
      YO <- model.response(mfO, "numeric")
      badRow <- badRow | !complete.cases(YO, XO)
      badRow <- badRow | is.infinite(YO)
      badRow <- badRow | apply(XO, 1, function(v) any(is.infinite(v)))
   }
   else
       stop("argument 'outcome' must be a formula")
   NXS <- ncol(XS)
   NXO <- ncol(XO)
   # Remove rows w/NA-s
   XS <- XS[!badRow,,drop=FALSE]
   YS <- YS[!badRow]
   XO <- XO[!badRow,,drop=FALSE]
   YO <- YO[!badRow]
   nObs <- length(YS)
   # few pre-calculations: split according to selection
   i0 <- YS == 0
   i1 <- YS == 1
   N0 <- sum(i0)
   N1 <- sum(i1)
   ## and run the model: selection
   probitResult <- probit(YS ~ XS - 1, maxMethod = maxMethod )
   if( print.level > 1) {
      cat("The probit part of the model:\n")
      print(summary(probitResult))
   }
   gamma <- coef(probitResult)
   ##
   z <- XS %*% gamma
   ## outcome
   invMillsRatio0 <- lambda(-z[i0])
   invMillsRatio1 <- lambda(z[i1])
   XO <- cbind(XO, .invMillsRatio=0)
   XO[i0,".invMillsRatio"] <- -invMillsRatio0
   XO[i1,".invMillsRatio"] <- invMillsRatio1
   ## if(checkIMRcollinearity(XO)) {
   ##    warning("Inverse Mills Ratio is virtually multicollinear to the rest of explanatory variables in the outcome equation")
   ## }
   olm <- lm(YO ~ -1 + XO)
                           # XO includes the constant (probably)
   if(print.level > 1) {
      cat("Raw outcome equation\n")
      print(summary(olm))
   }
   intercept <- any(apply(model.matrix(olm), 2,
                          function(v) (v[1] > 0) & (all(v == v[1]))))
                           # we have determine whether the outcome model has intercept.
                           # This is necessary later for calculating
                           # R^2
   delta0 <- mean( invMillsRatio0^2 - z[i0]*invMillsRatio0)
   delta1 <- mean( invMillsRatio1^2 + z[i1]*invMillsRatio1)
   betaL <- coef(olm)["XO.invMillsRatio"]
   sigma0.2 <- mean((residuals(olm)[i0])^2)*nObs/(nObs - NXS)
   sigma1.2 <- mean((residuals(olm)[i1])^2)*nObs/(nObs - NXS)
                           # residual variance: differs for
                           # treated/non-treated
   if(print.level > 2) {
      s2 <- 1
      rho <- 0.8
      th0.2 <- s2 + rho^2*s2*mean(z[i0]*invMillsRatio0) -
          rho^2*s2*mean(invMillsRatio0^2)
      th1.2 <- s2 - rho^2*s2*mean(z[i1]*invMillsRatio1) -
          rho^2*s2*mean(invMillsRatio1^2)
      a <- rbind(sd=c("non-participants"=sqrt(sigma0.2),
                 "participants"=sqrt(sigma1.2)),
                 th=c(sqrt(th0.2), sqrt(th1.2))
                 )
      cat("variances:\n")
      print(a)
   }
   sigma.02 <- sigma0.2 - betaL^2*mean(z[i0]*invMillsRatio0) +
       betaL^2*mean(invMillsRatio0^2)
   sigma.12 <- sigma1.2 + betaL^2*mean(z[i1]*invMillsRatio1) +
       betaL^2*mean(invMillsRatio1^2)
   ## take a weighted average over participants/non-participants
   sigma.2 <-
      (sum(i0)*sigma.02 + sum(i1)*sigma.12)/length(i1)
   sigma <- sqrt(sigma.2)
   rho <- betaL/sigma
   ## Now pack the results into a parameter vector
   ## indices in for the parameter vector
   iBetaS <- seq(length=NXS)
   iBetaO <- seq(tail(iBetaS, 1)+1, length=NXO)
   iMills <- tail(iBetaO, 1) + 1
                           # invMillsRatios are counted as parameter
   iSigma <- iMills + 1
   iRho <- tail(iSigma, 1) + 1
   nParam <- iRho
   ## Varcovar matrix.  Fill only a few parts, rest will remain NA
   coefficients <- numeric(nParam)
   coefficients[iBetaS] <- coef(probitResult)
   names(coefficients)[iBetaS] <- gsub("^XS", "",
                                       names(coef(probitResult)))
   coefficients[iBetaO] <- coef(olm)[names(coef(olm)) != "XO.invMillsRatio"]
   names(coefficients)[iBetaO] <- gsub("^XO", "",
                                       names(coef(olm))[names(coef(olm)) != "XO.invMillsRatio"])
   coefficients[c(iMills, iSigma, iRho)] <-
      c(coef(olm)["XO.invMillsRatio"], sigma, rho)
   names(coefficients)[c(iMills, iSigma, iRho)] <-
       c("invMillsRatio", "sigma", "rho")
   vc <- matrix(0, nParam, nParam)
   colnames(vc) <- row.names(vc) <- names(coefficients)
   vc[] <- NA
   if(!is.null(vcov(probitResult)))
       vc[iBetaS,iBetaS] <- vcov(probitResult)
   ## the 'param' component is intended to all kind of technical info
   param <- list(index=list(betaS=iBetaS,
                 betaO=iBetaO,
                 Mills=iMills, sigma=iSigma, rho=iRho,
                 errTerms = c(iMills, iSigma, iRho),
                 outcome = c(iBetaO, iMills) ),
                                        # The location of results in the coef vector
                 oIntercept1=intercept,
                 nObs=nObs, nParam=nParam, df=nObs-nParam + 2,
                 NXS=NXS, NXO=NXO, N0=N0, N1=N1,
                 levels=YSLevels
                           # levels[1]: selection 1; levels[2]: selection 2
                 )
   #
   result <- list(probit=probitResult,
                  lm=olm,
                  rho=rho,
                  sigma=sigma,
                  call = thisCall,
                  termsS=mtS,
                  termsO=mtO,
                  ys=switch(as.character(ys), "TRUE"=YS, "FALSE"=NULL),
                  xs=switch(as.character(xs), "TRUE"=XS, "FALSE"=NULL),
                  yo=switch(as.character(yo), "TRUE"=YO, "FALSE"=NULL),
                  xo=switch(as.character(xo), "TRUE"=XO, "FALSE"=NULL),
                  mfs=switch(as.character(mfs), "TRUE"=list(mfS), "FALSE"=NULL),
                  mfo=switch(as.character(mfs), "TRUE"=mfO, "FALSE"=NULL),
                  param=param,
                  coefficients=coefficients,
                  vcov=vc
                  )
   result$tobitType <- "treatment"
   result$method <- "2step"
   class( result ) <- c( "selection", class(result))
   return( result )
}
