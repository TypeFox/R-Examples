heckit5fit <- function(selection, outcome1, outcome2,
                    data=sys.frame(sys.parent()),
                    ys=FALSE, yo=FALSE,
                    xs=FALSE, xo=FALSE,
                    mfs=FALSE, mfo=FALSE,
                    print.level=0, maxMethod="Newton-Raphson", ... )
{
   ## internal function: 2-step estimator for tobit-5 model
   ## 
   ## maxMethod     maximization method for probit regression
   ## 
   checkIMRcollinearity <- function(X, tol=1e6) {
      ## This small utility checks whether inverse Mills ratio is (virtually) collinear to the other explanatory
      ## variables.  IMR is in the last columns.
      ## In case of collinearity it returns TRUE, otherwise FALSE
      X <- X[!apply(X, 1, function(row) any(is.na(row))),]
      if(kappa(X) < tol)
          return(FALSE)
      if(kappa(X[,-ncol(X)]) > tol)
          return(FALSE)
                                        # it is multicollinear, but not (just) due to IMR
      return(TRUE)
   }
   ## Do a few sanity checks...
   if( class( selection ) != "formula" ) {
      stop( "argument 'selection' must be a formula" )
   }
   if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   }
   thisCall <- match.call()
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
   badRow <- is.na(YS)
   badRow <- badRow | apply(XS, 1, function(v) any(is.na(v)))
   if("formula" %in% class( outcome1 )) {
      if( length( outcome1 ) != 3 ) {
         stop( "argument 'outcome1' must be a 2-sided formula" )
      }
      m <- match(c("outcome1", "data", "subset"), names(mf), 0)
      mf1 <- mf[c(1, m)]
      mf1$drop.unused.levels <- TRUE
      mf1$na.action <- na.pass
      mf1[[1]] <- as.name("model.frame")
      names(mf1)[2] <- "formula"
      mf1 <- eval(mf1, parent.frame())
      mt1 <- attr(mf1, "terms")
      XO1 <- model.matrix(mt1, mf1)
      YO1 <- model.response(mf1, "numeric")
      badRow <- badRow | (is.na(YO1) & (!is.na(YS) & YS == 0))
      badRow <- badRow | (apply(XO1, 1, function(v) any(is.na(v))) & (!is.na(YS) & YS == 0))
      if("formula" %in% class( outcome2 )) {
         if( length( outcome2 ) != 3 ) {
            stop( "argument 'outcome2' must be a 2-sided formula" )
         }
         m <- match(c("outcome2", "data", "subset"), names(mf), 0)
         mf2 <- mf[c(1, m)]
         mf2$drop.unused.levels <- TRUE
         mf2$na.action <- na.pass
         mf2[[1]] <- as.name("model.frame")
         names(mf2)[2] <- "formula"
         mf2 <- eval(mf2, parent.frame())
         mt2 <- attr(mf2, "terms")
         XO2 <- model.matrix(mt2, mf2)
         YO2 <- model.response(mf2, "numeric")
         badRow <- badRow | (is.na(YO2) & (!is.na(YS) & YS == 1))
         badRow <- badRow | (apply(XO2, 1, function(v) any(is.na(v))) & (!is.na(YS) & YS == 1))
      }
      else
          stop("argument 'outcome2' must be a formula")
   }
   else if("list" %in% class(outcome1)) {
      if(length(outcome1) != 2) {
         stop("argument 'outcome1' must be either a formula or a list of two formulas")
      }
      if("formula" %in% class(outcome1[[1]])) {
         if( length( outcome1[[1]] ) != 3 ) {
            stop( "argument 'outcome1[[1]]' must be a 2-sided formula" )
         }
      }
      else
          stop( "argument 'outcome1[[1]]' must be a formula" )
      if("formula" %in% class(outcome1[[2]])) {
         if( length( outcome1[[2]] ) != 3 ) {
            stop( "argument 'outcome[[2]]' must be a 2-sided formula" )
         }
         formula1 <- outcome1[[1]]
         formula2 <- outcome1[[2]]
                                        # Now we have extracted both formulas
         m <- match(c("outcome1", "data", "subset",
                      "offset"), names(mf), 0)
         ## replace the outcome list by the first equation and evaluate it
         oArg <- match("outcome1", names(mf), 0)
                                        # find the outcome argument
         mf[[oArg]] <- formula1
         mf1 <- mf[c(1, m)]
         mf1$drop.unused.levels <- TRUE
         mf1$na.action = na.pass
         mf1[[1]] <- as.name("model.frame")
                                        # eval it as model frame
         names(mf1)[2] <- "formula"
         mf1 <- eval(mf1, parent.frame())
         mt1 <- attr(mf1, "terms")
         XO1 <- model.matrix(mt1, mf1)
         YO1 <- model.response(mf1, "numeric")
         badRow <- badRow | (is.na(YO1) & (!is.na(YS) & YS == 0))
         badRow <- badRow | (apply(XO1, 1, function(v) any(is.na(v))) & (!is.na(YS) & YS == 0))
         ## repeat all the stuff with second equation
         mf[[oArg]] <- formula2
         mf2 <- mf[c(1, m)]
         mf2$drop.unused.levels <- TRUE
         mf2$na.action <- na.pass
         mf2[[1]] <- as.name("model.frame")
                                        # eval it as model frame
         names(mf2)[2] <- "formula"
         mf2 <- eval(mf2, parent.frame())
         mt2 <- attr(mf2, "terms")
         XO2 <- model.matrix(mt2, mf2)
         YO2 <- model.response(mf2, "numeric")
         badRow <- badRow | (is.na(YO2) & (!is.na(YS) & YS == 1))
         badRow <- badRow | (apply(XO2, 1, function(v) any(is.na(v))) & (!is.na(YS) & YS == 1))
      }
      else
          stop( "argument 'outcome[[2]]' must be a formula" )
   }
   else
       stop("argument 'outcome1' must be a formula or a list of two formulas")
   NXS <- ncol(XS)
   NXO1 <- ncol(XO1)
   NXO2 <- ncol(XO2)
   # Remove rows w/NA-s
   XS <- XS[!badRow,,drop=FALSE]
   YS <- YS[!badRow]
   XO1 <- XO1[!badRow,,drop=FALSE]
   YO1 <- YO1[!badRow]
   XO2 <- XO2[!badRow,,drop=FALSE]
   YO2 <- YO2[!badRow]
   nObs <- length(YS)
   # few pre-calculations: split according to selection
   i1 <- YS == 0
   i2 <- YS == 1
   XS1 <- XS[i1,,drop=FALSE]
   XS2 <- XS[i2,,drop=FALSE]
   XO1 <- XO1[i1,,drop=FALSE]
   XO2 <- XO2[i2,,drop=FALSE]
   YO1 <- YO1[i1]
   YO2 <- YO2[i2]
   N1 <- length(YO1)
   N2 <- length(YO2)
   ## and run the model: selection
   probitResult <- probit(YS ~ XS - 1, maxMethod = maxMethod )
   if( print.level > 0) {
      cat("The probit part of the model:\n")
      print(summary(probitResult))
   }
   gamma <- coef(probitResult)
   ## outcome
   invMillsRatio1 <- dnorm( -XS1%*%gamma)/pnorm( -XS1%*%gamma)
   invMillsRatio2 <- dnorm( XS2%*%gamma)/pnorm( XS2%*%gamma)
   colnames(invMillsRatio1) <- colnames(invMillsRatio2) <- "invMillsRatio"
                                        # lambdas are inverse Mills ratios
   XO1 <- cbind(XO1, invMillsRatio1)
   XO2 <- cbind(XO2, invMillsRatio2)
                                        # lambda1 is a matrix -- we need to remove the dim in order to
   if(checkIMRcollinearity(XO1)) {
      warning("Inverse Mills Ratio is virtually multicollinear to the rest of explanatory variables in the outcome equation 1")
   }
   if(checkIMRcollinearity(XO2)) {
      warning("Inverse Mills Ratio is virtually multicollinear to the rest of explanatory variables in the outcome equation 2")
   }
   lm1 <- lm(YO1 ~ -1 + XO1)
   lm2 <- lm(YO2 ~ -1 + XO2)
                                        # XO includes the constant
   intercept1 <- any(apply(model.matrix(lm1), 2,
                           function(v) (v[1] > 0) & (all(v == v[1]))))
   intercept2 <- any(apply(model.matrix(lm2), 2,
                           function(v) (v[1] > 0) & (all(v == v[1]))))
                                        # we have determine whether the outcome model has intercept.
                                        # This is necessary later for calculating R^2
   se1 <- summary(lm1)$sigma
   se2 <- summary(lm2)$sigma
                                        # residual variance
   delta1 <- mean( invMillsRatio1^2 - XS1%*%gamma *invMillsRatio1)
   delta2 <- mean( invMillsRatio2^2 + XS2%*%gamma *invMillsRatio2)
   betaL1 <- coef(lm1)["XO1invMillsRatio"]
   betaL2 <- coef(lm2)["XO2invMillsRatio"]
   sigma1 <- sqrt( se1^2 + ( betaL1*delta1)^2)
   sigma2 <- sqrt( se2^2 + ( betaL2*delta2)^2)
   rho1 <- -betaL1/sigma1
   rho2 <- betaL2/sigma2
   if( rho1 <= -1) rho1 <- -0.99
   if( rho2 <= -1) rho2 <- -0.99
   if( rho1 >= 1) rho1 <- 0.99
   if( rho2 >= 1) rho2 <- 0.99
   ## Now pack the results into a parameter vector
   ## indices in for the parameter vector
   iBetaS <- 1:NXS
   iBetaO1 <- seq(tail(iBetaS, 1)+1, length=NXO1)
   iMills1 <- tail(iBetaO1, 1) + 1
   iSigma1 <- iMills1 + 1
   iRho1 <- tail(iSigma1, 1) + 1
   iBetaO2 <- seq(tail(iRho1, 1) + 1, length=NXO2)
   iMills2 <- tail(iBetaO2, 1) + 1
   iSigma2 <- iMills2 + 1
   iRho2 <- tail(iSigma2, 1) + 1
   nParam <- iRho2
                                        # invMillsRatios are counted as parameter
   #
      ## Varcovar matrix.  Fill only a few parts, rest will remain NA
   coefficients <- numeric(nParam)
   coefficients[iBetaS] <- coef(probitResult)
   names(coefficients)[iBetaS] <- gsub("^XS", "", names(coef(probitResult)))
   coefficients[iBetaO1] <- coef(lm1)[names(coef(lm1)) != "XO1invMillsRatio"]
   names(coefficients)[iBetaO1] <- gsub("^XO1", "",
                                        names(coef(lm1))[names(coef(lm1)) != "XO1invMillsRatio"])
   coefficients[iBetaO2] <- coef(lm2)[names(coef(lm2)) != "XO2invMillsRatio"]
   names(coefficients)[iBetaO2] <- gsub("^XO2", "",
                                        names(coef(lm2))[names(coef(lm2)) != "XO2invMillsRatio"])
   coefficients[c(iMills1, iSigma1, iRho1, iMills2, iSigma2, iRho2)] <-
       c(coef(lm1)["XO1invMillsRatio"], sigma1, rho1,
         coef(lm2)["XO2invMillsRatio"], sigma2, rho2)
   names(coefficients)[c(iMills1, iSigma1, iRho1, iMills2, iSigma2, iRho2)] <-
       c("invMillsRatio1", "sigma1", "rho1", "invMillsRatio2", "sigma2", "rho2")
   vc <- matrix(0, nParam, nParam)
   colnames(vc) <- row.names(vc) <- names(coefficients)
   vc[] <- NA
   if(!is.null(vcov(probitResult)))
       vc[iBetaS,iBetaS] <- vcov(probitResult)
   ## the 'param' component is intended to all kind of technical info
   param <- list(index=list(betaS=iBetaS,
                 betaO1=iBetaO1, betaO2=iBetaO2,
                 Mills1=iMills1, sigma1=iSigma1, rho1=iRho1,
                 Mills2=iMills2, sigma2=iSigma2, rho2=iRho2,
                 errTerms = c( iMills1, iMills2, iSigma1, iSigma2, iRho1, iRho2 ),
                 outcome = c( iBetaO1, iMills1, iBetaO2, iMills2 ) ),
                                        # The location of results in the coef vector
                 oIntercept1=intercept1, oIntercept2=intercept2,
                 nObs=nObs, nParam=nParam, df=nObs-nParam + 2,
                 NXS=NXS, NXO1=NXO1, NXO2=NXO2, N1=N1, N2=N2,
                 levels=YSLevels
                           # levels[1]: selection 1; levels[2]: selection 2
                 )
   #
   result <- list(probit=probitResult,
                  lm1=lm1,
                  rho1=rho1,
                  sigma1=sigma1,
                  lm2=lm2,
                  rho2=rho2,
                  sigma2=sigma2,
                  call = thisCall,
                  termsS=mtS,
                  termsO=list(mt1, mt2),
                  ys=switch(as.character(ys), "TRUE"=YS, "FALSE"=NULL),
                  xs=switch(as.character(xs), "TRUE"=XS, "FALSE"=NULL),
                  yo=switch(as.character(yo), "TRUE"=list(YO1, YO2), "FALSE"=NULL),
                  xo=switch(as.character(xo), "TRUE"=list(XO1, XO2), "FALSE"=NULL),
                  mfs=switch(as.character(mfs), "TRUE"=list(mfS), "FALSE"=NULL),
                  mfo=switch(as.character(mfs), "TRUE"=list(mf1, mf2), "FALSE"=NULL),
                  param=param,
                  coefficients=coefficients,
                  vcov=vc
                  )
   result$tobitType <- 5
   result$method <- "2step"
   class( result ) <- c( "selection", class(result))
   return( result )
}
