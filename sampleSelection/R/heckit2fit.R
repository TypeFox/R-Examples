heckit2fit <- function( selection, outcome,
                   data=sys.frame(sys.parent()),
                   weights = NULL, inst = NULL,
                   print.level = 0, maxMethod="Newton-Raphson" ) {
   ## selection     formula, selection equation
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
   ## What is the role of na.action here?  We cannot use na.omit -- we must not omit the observation
   # where outcome is not observed.  na-s cannot be passed either.
   # However, we can (and should?) omit the na-s in explanatory and probit outcomes.  This needs
   # a bit of refinement.
   if( class( outcome ) != "formula" ) {
      stop( "argument 'outcome' must be a formula" )
   } else if( length( outcome ) != 3 ) {
      stop( "argument 'outcome' must be a 2-sided formula" )
   } else if( "invMillsRatio" %in% all.vars( outcome ) ) {
      stop( "argument 'outcome' may not include a variable name",
         " 'invMillsRatio'" )
   } else if( class( selection ) != "formula" ) {
      stop( "argument 'selection' must be a formula" )
   } else if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   } else if( "probit" %in% substr( all.vars( selection ), 1, 6 ) ) {
      stop( "argument 'selection' may not include a variable",
         " names starting with 'probit'" )
   } else if( !is.null( inst ) ) {
      if( class ( inst ) != "formula" || length( inst ) != 2 ) {
         stop( "argument 'inst' must be a 1-sided formula" )
      }
   }
   result <- list()
   result$call <- match.call()
   ## Now extract model frames etc.
   ## Selection equation
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("selection", "data", "subset",
                "offset"), names(mf), 0)
   mfS <- mf[c(1, m)]
   mfS$na.action <- na.pass
   mfS$drop.unused.levels <- TRUE
   mfS[[1]] <- as.name("model.frame")
   names(mfS)[2] <- "formula"
                                        # model.frame requires the parameter to
                                        # be 'formula'
   mfS <- eval(mfS, parent.frame())
   mtS <- attr(mfS, "terms")
   XS <- model.matrix(mtS, mfS)
   NXS <- ncol(XS)
   YS <- model.response(mfS)
   badRow <- is.na(YS)
   badRow <- badRow | apply(XS, 1, function(v) any(is.na(v)))
                                        # check for NA-s.  Because we have to find NA-s in several
                                        # frames, we cannot use the standard na.* functions here.
                                        # Find bad rows and remove them later.
   probitLevels <- levels( as.factor( YS ) )
   if( length( probitLevels ) != 2 ) {
      stop( "the dependent variable of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   ysNames <- names( YS )
   YS <- as.integer(YS == probitLevels[ 2 ])
                                        # selection kept as integer internally
   names( YS ) <- ysNames

   ## Outcome equation
   m <- match(c("outcome", "data", "subset",
                "offset"), names(mf), 0)
   mfO <- mf[c(1, m)]
   mfO$na.action <- na.pass
   mfO$drop.unused.levels <- TRUE
   mfO$na.action <- na.pass
                                        # Here we have to keep NA-s: unobserved outcome variables may
                                        # be marked as NA-s.
   mfO[[1]] <- as.name("model.frame")
                                        # eval it as model frame
   names(mfO)[2] <- "formula"
   mfO <- eval(mfO, parent.frame())
   mtO <- attr(mfO, "terms")
   XO <- model.matrix(mtO, mfO)
                                        # the explanatory variables in matrix form
   NXO <- ncol(XO)
   YO <- model.response(mfO)
   if(is.factor(YO))
       YO <- as.numeric(as.character(YO))
   ## Remove NA observations
   badRow <- badRow | (is.na(YO) & (!is.na(YS) & YS == 1))
   badRow <- badRow | (apply(XO, 1, function(v) any(is.na(v))) & (!is.na(YS) & YS == 1))
                                        # rows in outcome, which contain NA and are observable -> bad too
   
   if( !is.null( weights ) ) {
      if( length( weights ) != length( badRow ) ) {
         stop( "number of weights (", length( weights ), ") is not equal",
            " to the number of observations (", length( badRow ), ")" )
      }
      badRow <- badRow | is.na( weights )
   }   
   
   if(print.level > 0) {
      cat(sum(badRow), "invalid observations\n")
   }
   XS <- XS[!badRow,,drop=FALSE]
   YS <- YS[!badRow]
   XO <- XO[!badRow,,drop=FALSE]
   YO <- YO[!badRow]
   weightsNoNA <- weights[ !badRow ]
   ## Now indices for packing the separate outcomes into full outcome vectors.  Note we treat
   ## invMillsRatio as a separate parameter
   iBetaS <- seq(length=NXS)
   iBetaO <- NXS + seq(length=NXO)
   iMills <- NXS + NXO + 1
   iSigma <- iMills + 1
   iRho <- iSigma + 1
   ##
   nObs <- length(YS)
   nParam <- iRho
   N0 <- sum(YS == 0)
   N1 <- nObs - N0
                                        # sigma, rho
   if( print.level > 0 ) {
      cat ( "\nEstimating 1st step Probit model . . ." )
   }
   result$probit <- probit(YS ~ XS - 1, x=TRUE, 
      weights = weightsNoNA, print.level=print.level - 1, iterlim=30,
      maxMethod = maxMethod )
                                        # a large iterlim may help with weakly identified models
   if( print.level > 0 ) {
       cat( " OK\n" )
   }
   imrData <- invMillsRatio( result$probit )
   result$imrDelta <- drop(imrData$delta1)
   ## ---- Outcome estimation -----
   if( is.null( inst ) ) {
      if( print.level > 0 ) {
         cat ( "Estimating 2nd step (outcome) OLS model . . ." )
      }
      if(checkIMRcollinearity(cbind(XO, imrData$IMR1))) {
         warning("Inverse Mills Ratio is (virtually) collinear to the rest of the explanatory variables")
      }
      outcomeMod <- lm(YO ~ -1 + XO + imrData$IMR1, weights = weightsNoNA,
                      subset = YS == 1 )
      intercept <- any(apply(model.matrix(outcomeMod), 2,
                             function(v) (v[1] > 0) & (all(v == v[1]))))
                                        # we have determine whether the outcome model has intercept.
                                        # This is necessary later for calculating R^2
      resid <- residuals( outcomeMod )
      step2coef <- coef( outcomeMod )
      names(step2coef) <- c(colnames(XO), "invMillsRatio")
      if(print.level > 0)
          cat( " OK\n" )
   }
   else {
      if( !is.null( weights ) ) {
         warning( "weights are ignored in instrumental variable estimations" )
      }
      data$invMillsRatio <- imrData$IMR1
      if( print.level > 0 ) {
         cat ( "Estimating 2nd step 2SLS/IV model . . ." )
      }
      step2formula <- as.formula( paste( outcome[ 2 ], "~", outcome[ 3 ],
                                        "+ invMillsRatio" ) )
      formulaList <- list( step2formula )
      instImr <- as.formula( paste( "~", inst[ 2 ], "+ invMillsRatio" ) )
      outcomeMod <- systemfit( formulaList, method = "2SLS", inst = instImr,
                             data = data[ YS == 1, ] )
      intercept = FALSE
                                        # we calculate R^2 differently here (hopefully)
      resid <- residuals( outcomeMod )[ , 1 ]
      step2coef <- coefficients( outcomeMod$eq[[ 1 ]] )
      if( print.level > 0 ) cat( " OK\n" )
   }
   result$sigma <- drop( sqrt( crossprod( resid ) /
                              sum( YS ) +
                              mean(result$imrDelta[ YS == 1] ) *
                              step2coef[ "invMillsRatio" ]^2 ) )
   result$rho <-  step2coef[ "invMillsRatio" ] / result$sigma
   names(result$rho) <- NULL
   result$invMillsRatio <- drop( imrData$IMR1 )
   ## Stack all final coefficients to 'coefficients'
   coefficients <- c(coef(result$probit),
                     step2coef[names(step2coef) != "invMillsRatio"],
                     step2coef["invMillsRatio"],
                     sigma=result$sigma, rho=result$rho)
   names(coefficients)[iBetaS] <- gsub("^XS", "", names(coefficients)[iBetaS])
   if( print.level > 0 ) {
      cat ( "Calculating coefficient covariance matrix . . ." )
   }
                                        # the following variables are named according to Greene (2003), p. 785
   xMat <- model.matrix( outcomeMod )
   ## Varcovar matrix.  Fill only a few parts, rest will remain NA
   vc <- matrix(0, nParam, nParam)
   colnames(vc) <- row.names(vc) <- names(coefficients)
   vc[] <- NA
   if(!is.null(vcov(result$probit)))
       vc[iBetaS,iBetaS] <- vcov(result$probit)
   vc[c(iBetaO, iMills), c(iBetaO, iMills)] <- heckitVcov( xMat,
                                                          model.matrix( result$probit )[ YS == 1, ],
                                                          vcov( result$probit ),
                                                          result$rho,
                                                          result$imrDelta[ YS == 1],
                                                          result$sigma )
                                        # here we drop invMillsRatio part
   result$vcov <- vc
   ##
   if( print.level > 0 )
       cat( " OK\n" )
   result$coefficients <- coefficients
                                        # for coef() etc. methods
   ## the 'param' component is intended to all kind of technical info
   result$param <- list(index=list(betaS=iBetaS, betaO=iBetaO, 
                        Mills=iMills, sigma=iSigma, rho=iRho,
                        errTerms = c( iMills, iSigma, iRho ),
                        outcome = c( iBetaO, iMills ) ),
                                        # The location of results in the coef vector
                        oIntercept=intercept,
                        N0=N0, N1=N1,
                        nParam=nParam, nObs=nObs, df=nObs-nParam+1)
   result$lm <- outcomeMod
   result$tobitType <- 2
   result$method <- "2step"
   result$weights <- weights
   class( result ) <- c( "selection", class(result))
   return( result )
}
