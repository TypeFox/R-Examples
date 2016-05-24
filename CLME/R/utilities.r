

##
## Various support functions to assist the main functions
##

#' Create model matrices for \code{clme}
#'
#' @description
#' Parses formulas to creates model matrices for \code{clme}. 
#'
#' @param formula a formula defining a linear fixed or mixed effects model. The constrained effect(s) must come before any unconstrained covariates on the right-hand side of the expression. The first \code{ncon} terms will be assumed to be constrained. 
#' @param data data frame containing the variables in the model.
#' @param ncon the number of variables in \code{formula} that are constrained.
#' 
#' 
#' @note
#' The first term on the right-hand side of the formula should be the fixed effect 
#' with constrained coefficients. Random effects are represented with a vertical bar, 
#' so for example the random effect \code{U} would be included by 
#' \code{Y ~ X1 + (1|U)}.
#' 
#' The intercept is removed automatically. This is done to ensure that parameter 
#' estimates are of the means of interest, rather than being expressed as a mean 
#' with offsets.
#' 
#' @return
#' A list with the elements:
#' \tabular{rl}{
#'   Y       \tab response variable \cr
#'   X1      \tab design matrix for constrained effect \cr
#'   X2      \tab design matrix for covariates \cr
#'   P1      \tab number of constrained coefficients \cr
#'   U       \tab matrix of random effects \cr
#'   formula \tab the final formula call (automatically removes intercept) \cr
#'   dframe  \tab the dataframe containing the variables in the model \cr
#'   REidx   \tab an element to define random effect variance components \cr
#'   REnames \tab an element to define random effect variance components \cr
#' }
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' data( rat.blood )
#' model_terms_clme( mcv ~ time + temp + sex + (1|id) , data = rat.blood )
#' 
#' @importFrom lme4 lFormula
#' @export
#' 
model_terms_clme <- function( formula, data, ncon=1 ){
  
  formula2 <- update.formula( formula , . ~ . - 1 )
  
  nn    <- nrow(data)
  Unull <- c( rep("RemoveAAA",round(nn/2)), rep("RemoveBBB", nn-round(nn/2)) )
  
  data  <- cbind( data , Unull )
  
  formula3 <- update.formula( formula , . ~ . + (1|Unull) - 1)
  
  suppressMessages( mterms   <- lme4::lFormula( formula3, data=data ) )
  
  Y  <- mterms$fr[,1]
  X  <- mterms$X
  P1 <- sum( attr(X, "assign") <= ncon )
  X1 <- X[,     1:P1     , drop=FALSE]
  
  if( ncol(X) > P1 ){
    X2 <- X[,(P1+1):ncol(X), drop=FALSE]  
  } else{
    X2 <- NULL
  }
  
  U            <- t( as.matrix(mterms$reTrms$Zt) )
  dframe       <- mterms$fr
  dframe$Unull <- NULL
  
  if( ncol(U)==2 ){
    mmat <- list( Y=Y, X1=X1, X2=X2, P1=P1, U=NULL, formula=formula2, dframe=dframe, 
                  REidx=NULL, REnames=NULL )
  } else{
    drop.col <- which( colnames(U) %in% c("RemoveAAA","RemoveBBB") )
    drop.nam <- which( colnames(mterms$reTrms$flist)=="Unull" )
    
    U <- U[, -drop.col]
    
    REnames <- colnames(mterms$reTrms$flist)[-drop.nam]
    REidx   <- mterms$reTrms$Lind[-drop.col]
    
    mmat <- list( Y=Y, X1=X1, X2=X2, P1=P1, U=U, formula=formula2, dframe=dframe, 
                  REidx=REidx, REnames=REnames )
  }
  
  return( mmat )
  
}


##
## Some methods for class CLME
##

#' Constructor method for objects S3 class clme
#' 
#' @rdname as.clme
#' @export
#' 
is.clme <- function(x) inherits(x, "clme")


#' Constructor method for objects S3 class clme
#'
#' @description
#' Test if an object is of class \code{clme} or coerce an object to be such.
#' 
#' @rdname as.clme
#' 
#' @param x list with the elements corresponding to the output of \code{\link{clme}}.
#' @param ... space for additional arguments.
#' 
#' @return
#' Returns an object of the class \code{clme}.
#' 
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' data( rat.blood )
#' 
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#' 
#' is.clme( clme.out )
#' as.clme( clme.out )
#' @export
#' 
as.clme <- function( x , ... ){
  
  if( is.clme(x) ){
    return(x)
  } else{
    
    err.flag  <- 0
    flagTheta <- flagSsq <- flagTsq <- flagCov <- flagW1 <- flagW2 <- flagP1 <- flagP2 <- flagConst <- ""
    
    if( !is.numeric(x$theta) ){
      err.flag  <- 1
      flagTheta <- " theta must be numeric \n"
      x$theta   <- numeric(0)
    }
    
    if( !is.numeric(x$ssq) ){
      err.flag <- 1
      flagSsq  <- " ssq must be numeric \n"
      x$ssq    <- numeric(0)
    } 
    
    if( !is.null(x$tsq) & !is.numeric(x$tsq) ){
      err.flag <- 1
      flagTsq  <- " if present, tau must be numeric \n"
      x$tsq    <- NULL
    }
    
    if( !is.matrix(x$cov.theta) || !is.numeric(x$cov.theta) ||
          nrow(x$cov.theta) != ncol(x$cov.theta) ||
          nrow(x$cov.theta) != length(x$theta)   ||
          sum(sum(abs(x$cov.theta - t(x$cov.theta)))) > sqrt(.Machine$double.eps) ){
      err.flag    <- 1
      flagCov     <- " cov.theta must be square, symmetric, numeric matrix with dimensions equal to length of theta\n"
      x$cov.theta <- matrix( numeric(0) , nrow=length(x$theta) , ncol=length(x$theta) )
    }
    
    if( !is.numeric(x$ts.glb) ){
      err.flag <- 1
      flagW1   <- " ts.glb must be numeric \n"
      x$ts.glb  <- numeric(0)
    } 
    
    if( !is.numeric(x$ts.ind) ){
      err.flag <- 1
      flagW2   <- " ts.ind must be numeric \n"
      x$ts.ind  <- numeric(0)
    }
    
    if( !is.numeric(x$p.value) || length(x$p.value) != length(x$ts.glb) ){
      err.flag   <- 1
      flagP1     <- " p.value must be numeric and of same length as ts.glb \n"
      x$p.value  <- numeric(0)
    } 
    
    if( !is.numeric(x$p.value.ind) || length(x$p.value.ind) != length(x$ts.ind) ){
      err.flag       <- 1
      flagP2         <- " p.value.ind must be numeric and of same length as ts.ind \n"
      x$p.value.ind  <- numeric(0)
    } 
    
    if( !is.list(x$constraints) ){
      err.flag        <- 1
      flagConst       <- " constraints must be list \n"
      x$constraints   <- list( A=matrix( numeric(0) ) )
    } else{
      cnames <- names(x$constraints)
      if( sum(cnames=="A") != 1 ){
        err.flag        <- 1
        flagConst       <- " constraints must contain element A\n"
        x$constraints$A <- matrix( numeric(0) , nrow=length(x$ts.ind ))
      }
    }
    
    
    if( err.flag==1 ){
      err.mssg <- paste( "coercing 'x' to class 'clme' produced errors: \n", 
                         flagTheta, flagSsq, flagTsq, flagCov, flagW1,
                         flagW2, flagP1, flagP2, flagConst, "output may not be valid." , sep = "")
      # warning(warn, sys.call(-1))
      warning( err.mssg )
    }
    
    class(x) <- "clme"
    return(x)
    
  }
}


################################################################################

#' Akaike information criterion
#'
#' @description
#' Calculates the Akaike and Bayesian information criterion for objects of class \code{clme}. 
#' 
#' @param object object of class \code{\link{clme}}.
#' @param ... space for additional arguments.
#' @param k value multiplied by number of coefficients
#' 
#' @details
#' The log-likelihood is assumed to be the Normal distribution. The model uses residual bootstrap methodology, and Normality is neither required nor assumed. Therefore the log-likelihood and these information criterion may not be useful measures for comparing models.
#' For \code{k=2}, the function computes the AIC. To obtain BIC, set \eqn{k = log( n/(2*pi) )}; which the method \code{BIC.clme} does.
#' 
#' 
#' @return
#' Returns the information criterion (numeric).
#' 
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' 
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#' 
#' AIC( clme.out )
#' AIC( clme.out, k=log( nobs(clme.out)/(2*pi) ) )
#' 
#' 
#' @method AIC clme
#' @export
#' 
AIC.clme <- function( object, ..., k=2 ){
  ## For BIC, set k = ln( n/(2*pi) )
  logl <- logLik.clme( object, ...)[1]
  kk   <- ncol(model.matrix.clme(object)) + length(object$tsq) + length(object$ssq) 
  # aic  <- 2*(kk - logl)
  aic  <- k*kk - 2*logl
  return(aic)
}

#' Akaike information criterion
#'
#' @description
#' Calculates the Akaike and Bayesian information criterion for objects of class \code{clme}. 
#' 
#' @rdname AIC.clme
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @method AIC summary.clme
#' @export
#' 
AIC.summary.clme <- function( object, ..., k=2 ){
  class(object) <- "clme"
  AIC( object, ..., k=k )
}



#' Bayesian information criterion
#'
#' @description
#' Calculates the Bayesian information criterion for objects of class \code{clme}. 
#' 
#' @param object object of class \code{\link{clme}}.
#' @param ... space for additional arguments.
#' @param k value multiplied by number of coefficients
#' 
#' @details
#' The log-likelihood is assumed to be the Normal distribution. The model uses residual bootstrap methodology, and Normality is neither required nor assumed. Therefore the log-likelihood and these information criterion may not be useful measures for comparing models.
#' For \code{k=2}, the function computes the AIC. To obtain BIC, set \eqn{k = log( n/(2*pi) )}; which the method \code{BIC.clme} does.
#' 
#' 
#' @return
#' Returns the Bayesian information criterion (numeric).
#' 
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' 
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#' 
#' BIC( clme.out )
#' BIC( clme.out, k=log( nobs(clme.out)/(2*pi) ) )
#' 
#' 
#' @method BIC clme
#' @export
#' 
BIC.clme <- function( object, ..., k=log(nobs(object)/(2*pi)) ){
  ## For BIC, set k = ln( n/(2*pi) )
  logl <- logLik( object, ...)[1]
  bic <- AIC( object, k=k )
  return(bic)
}

#' Bayesian information criterion
#'
#' @description
#' Calculates the Akaike and Bayesian information criterion for objects of class \code{clme}. 
#' 
#' @rdname BIC.clme
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @method BIC summary.clme
#' @export
#' 
BIC.summary.clme <- function( object, ..., k=log(nobs(object)/(2*pi)) ){
  class(object) <- "clme"
  BIC( object, ..., k=k )
}


#' Individual confidence intervals
#'
#' @description
#' Calculates confidence intervals for fixed effects parameter estimates in objects of class \code{clme}.
#' @rdname confint
#' 
#' @param object object of class \code{\link{clme}}.
#' @param parm parameter for which confidence intervals are computed (not used).
#' @param level nominal confidence level.
#' @param ... space for additional arguments.
#' 
#' @details
#' Confidence intervals are computed using Standard Normal critical values. 
#' Standard errors are taken from the covariance matrix of the unconstrained parameter estimates.
#' 
#' 
#' @return
#' Returns a matrix with two columns named lcl and ucl (lower and upper confidence limit).
#' 
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#' 
#' confint( clme.out )
#' 
#' 
#' @method confint clme
#' @export
#' 
confint.clme <- function(object, parm, level=0.95, ...){
  ## More types of confidence intervals (e.g., bootstrap) may be added in the future.
  ## If so, default confidence interval will be the current methods
  
  ## Confidence intervals based on unconstrained variance-covariance matrix
  ## Actual covarage >= Nominal coverage
  
  cc     <- match.call()
  digits <- cc$digits
  if( is.null(digits) ){     digits <- 3 }
  if( !is.numeric(digits) ){ digits <- 3 }
  if( digits < 0  ){         digits <- 3 }
  
  
  alpha <- 1 - level
  theta <- fixef(object)
  cv    <- qnorm(1-alpha/2)
  varco <- vcov( object )
  lcl   <- as.numeric( format( round(theta - cv*sqrt(diag(varco)) , digits=digits)) )
  ucl   <- as.numeric( format( round(theta + cv*sqrt(diag(varco)) , digits=digits)) )
  
  ints <- cbind( lcl, ucl)
  
  ## Return intervals
  return( ints )
}

#' Individual confidence intervals
#'
#' @description
#' Calculates confidence intervals for fixed effects parameter estimates in objects of class \code{clme}.
#' @rdname confint
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @method confint summary.clme
#' @export
#' 
confint.summary.clme <- function(object, parm, level=0.95, ...){
  class(object) <- "clme"
  confint( object, parm, level, ... )
}


#' Extract fixed effects
#' 
#' @importFrom lme4 fixef
#' @method fixef clme
#' @export 
#' 
fixef.clme <- function( object, ...){ UseMethod("fixef") }


#' Extract fixed effects
#' 
#' @rdname fixef.clme
#' @importFrom lme4 fixef
#' @method fixef summary.clme
#' @export 
#' 
fixef.summary.clme <- function( object, ...){ 
  class(object) <- "clme"
  fixef(object, ...)
}


#' Extract fixed effects
#'
#' @description
#' Extracts the fixed effects estimates from objects of class \code{clme}. 
#' 
#' @rdname fixef.clme
#' 
#' @param object object of class \code{\link{clme}}.
#' @param ... space for additional arguments
#' 
#' @return
#' Returns a numeric vector.
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#' 
#' fixef( clme.out )
#' 
#' @importFrom lme4 fixef
#' @method fixef clme
#' @export
#' 
fixef.clme <- function( object, ... ){
  ## Print out the fixed effects
  if( is.clme(object) ){
    return( object$theta )
  } else{
    stop("'object' is not of class clme")
  }
}

#' Extract fixed effects
#' 
#' @rdname fixef.clme
#' @importFrom nlme fixed.effects
#' @export
#' 
fixed.effects <- function( object, ...){ UseMethod("fixed.effects") }

#' Extract fixed effects
#' 
#' @rdname fixef.clme
#' @importFrom nlme fixed.effects
#' @export
#' 
fixed.effects.summary.clme <- function( object, ...){
  class(object) <- "clme"
  fixef(object, ...)
}


#' 
#' @rdname fixef.clme
#' 
#' @importFrom nlme fixed.effects
#' @method fixed.effects clme
#' @export
#' 
fixed.effects.clme <- function( object , ... ){
  fixef.clme( object, ... )
}

#' @rdname fixef.clme
#' @method coefficients clme
#' @export
#' 
fixed.effects.clme <- function( object , ... ){
  fixef.clme( object, ... )
}

#' @rdname fixef.clme
#' @method coefficients clme
#' @export
#' 
coefficients.clme <- function( object, ... ){
  fixef.clme( object, ... )
}

#' @rdname fixef.clme
#' @method coef clme
#' @export
#' 
coef.clme <- function( object, ... ){
  fixef.clme( object, ... )
}

#' @rdname fixef.clme
#' @method coefficients summary.clme
#' @export
#' 
coefficients.summary.clme <- function( object, ... ){
  class(object) <- "clme"
  fixef.clme( object, ... )
}

#' @rdname fixef.clme
#' @method coef summary.clme
#' @export
#' 
coef.summary.clme <- function( object, ... ){
  class(object) <- "clme"
  fixef.clme( object, ... )
}


#' Extract formula
#'
#' @description
#' Extracts the formula from objects of class \code{clme}. 
#' 
#' @param x object of class \code{\link{clme}}.
#' @param ... space for additional arguments
#' 
#' @details
#' The package \pkg{CLME} parametrizes the model with no intercept term. 
#' If an intercept was included, it will be removed automatically.
#' 
#' @return
#' Returns a formula object
#' 
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#' 
#' formula( clme.out )
#' 
#' @method formula clme
#' @export
#' 
formula.clme <- function(x, ...){
  return( x$formula )
}


#' Log-likelihood
#'
#' @description
#' Computes the log-likelihood of the fitted model for objects of class \code{clme}. 
#' 
#' @param object object of class \code{\link{clme}}.
#' @param ... space for additional arguments
#' 
#' @details
#' The log-likelihood is computed using the Normal distribution. The model uses residual bootstrap 
#' methodology, and Normality is neither required nor assumed. Therefore the log-likelihood may 
#' not be a useful measure in the context of \pkg{CLME}.
#' 
#' @return
#' Numeric.
#' 
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#' 
#' logLik( clme.out )
#' 
#' @method logLik clme
#' @export
#' 
logLik.clme <- function( object, ...){
  ## Residuals
  YY <- object$dframe[,1]
  XX <- model.matrix( object )
  TT <- fixef( object )
  RR <- YY - apply( XX , 1 , FUN=function(xx,tht){ sum(xx*tht) }, tht=TT )
  nn <- nobs(object)
  
  ## Covariance matrix (piecewise)
  ssq <- object$ssq
  Nks <- object$gfix
  ssqvec <- rep( ssq, Nks )
  RSiR <- c( t(RR) %*% (RR/ssqvec) )
  detS <- sum( log(ssqvec) )
  
  if( is.null(object$tsq) ){
    ## Fixed effects only
    RPiR   <- RSiR
    detPhi <- detS
  } else{
    ## Mixed Effects    
    UU     <- model.matrix( object , type="ranef" )
    tsq    <- object$tsq
    Qs     <- object$gran
    tsqvec <- rep( tsq, Qs )
    RSiU   <- matrix( apply( UU, 2, FUN=function(uu,rr){ sum(uu*rr) }, rr=(RR/ssqvec) ), nrow=1 )
    
    U1         <- apply( UU, 2, FUN=function(uu,sq){uu/sq}, sq=sqrt(ssqvec) )
    tusu       <- t(U1) %*% U1
    diag(tusu) <- diag(tusu) + 1/tsqvec
    tusui      <- solve( tusu )
    
    RPiR   <- RSiR - c( RSiU%*%(tusui%*%t(RSiU)) )
    detPhi <- log(det( tusu )) + sum( log(tsqvec) ) + detS
  }
  
  logL <- -0.5*( nn*log(2*pi) + detPhi + RPiR )
  
  return( logL )
  
}

#' Log-likelihood
#' 
#' @rdname logLik.clme
#' 
#' @seealso
#' \code{\link{logLik.clme}} 
#' 
#' @method logLik summary.clme
#' @export
#' 
logLik.summary.clme <- function( object, ...){
  class(object) <- "clme"
  logLik(object, ...)
}

#' Extracts the model frame
#'
#' @description
#' Extracts the model frame from objects of class \code{clme}. 
#' 
#' @param formula a formula expression.
#' @param ... space for additional arguments
#' 
#' 
#' @return
#' Returns a data frame with the variables in the model.
#' 
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' \dontrun{
#' data( rat.blood )
#' model.frame.clme( mcv ~ time + temp + sex + (1|id), data = rat.blood )
#' }
#' @method model.frame clme
#' @export
#' 
model.frame.clme <- function( formula , ...){
  ## Return the data frame
  mmat     <- model_terms_clme( formula, ... )  
  return( mmat$dframe )
}

#' Extracts the model frame
#' 
#' @rdname model.frame.clme
#' 
#' @seealso
#' \code{\link{model.frame.clme}} 
#' 
#' @method model.frame summary.clme
#' @export
#' 
model.frame.summary.clme <- function( formula, ...){
  class(object) <- "clme"
  model.frame(formula, ...)
}



#' Extract the model design matrix.
#'
#' @description
#' Extracts the fixed-effects design matrix from objects of class \code{clme}.
#' 
#' @param object an object of class \code{clme}.
#' @param type specify whether to return the fixed-effects or random-effects matrix.
#' @param ... space for additional arguments
#' 
#' 
#' @return
#' Returns a matrix.
#' 
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' \dontrun{
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#' 
#' model.matrix( clme.out )
#' }
#' @method model.matrix clme
#' @export
#' 
model.matrix.clme <- function( object, type="fixef", ...){
  mmat <- model_terms_clme( object$formula, object$dframe )  
  if( type=="fixef" ){
    ## Return the fixed-effects matrix
    X1   <- mmat$X1
    X2   <- mmat$X2
    return( cbind(X1, X2) )  
  } else if( type=="ranef" ){
    ## Return the random-effects matrix
    return(mmat$U)
  }
}

#' Extract the model design matrix.
#' 
#' @rdname model.matrix.clme
#' 
#' @seealso
#' \code{\link{model.matrix.clme}} 
#' 
#' @method model.matrix summary.clme
#' @export
#' 
model.matrix.summary.clme <- function( object, ...){
  class(object) <- "clme"
  model.matrix(object, ...)
}

#' Number of observations
#'
#' @description
#' Obtains the number of observations used to fit an model for objects of class \code{clme}.
#' 
#' @param object an object of class \code{clme}.
#' @param ... space for additional arguments
#' 
#' 
#' @return
#' Numeric.
#' 
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#' 
#' nobs( clme.out )
#' 
#' @method nobs clme
#' @export
#' 
nobs.clme <- function(object, ...){
  nrow( model.matrix.clme(object) )
}

#' Number of observations
#' 
#' @rdname nobs.clme
#' 
#' @seealso
#' \code{\link{nobs.clme}} 
#' 
#' @method nobs summary.clme
#' @export
#' 
nobs.summary.clme <- function( object, ...){
  class(object) <- "clme"
  nobs(object, ...)
}



#' Printout of fitted object.
#'
#' @description
#' Prints basic information on a fitted object of class \code{clme}.
#' 
#' @param x an object of class \code{clme}.
#' @param ... space for additional arguments
#' 
#' @return
#' Text printed to console.
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' \dontrun{
#' data( rat.blood )
#' set.seed( 42 )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 10)
#' 
#' print( clme.out )
#' }
#' 
#' @method print clme
#' @export
#' 
print.clme <- function(x, ...){
  #cc     <- match.call()
  #digits <- cc$digits
  
  object <- x
  ## Print out residuals of specified type
  if( !is.clme(object) ){
    stop("'object' is not of class clme")
  }
  
  cat( "Linear mixed model subject to order restrictions\n")
  cat( "Formula: ")
  print( object$formula )  
  crit <- c(logLik.clme(object),
            AIC.clme(object),
            AIC( object, k=log(nobs.clme(object)/(2*pi)) ) )  
  
  critc <- format( crit , digits=5)
  
  cat( "\nlog-likelihood:", critc[1] )
  cat( "\nAIC:           ", critc[2] )
  cat( "\nBIC:           ", critc[3] )
  cat( "\n(log-likelihood, AIC, BIC computed under normality)")  
  
  cat( "\n\nFixed effect coefficients (theta): \n")
  print( fixef.clme(object)  )
  
  cat( "\nVariance components: \n")
  print( VarCorr.clme(object) )
  #cat( "\n\nModel based on", object$nsim, "bootstrap samples." )
}


#' Extract random effects
#' 
#' @param object object of class clme.
#' @param ... space for additional arguments
#' 
#' @rdname ranef.clme
#' @importFrom nlme ranef
#' @export
#' 
ranef.clme <- function( object, ...){ UseMethod("ranef") }

#' Extract random effects
#' 
#' @rdname ranef.clme
#' 
#' @importFrom nlme ranef
#' @export
#' 
ranef.summary.clme <- function( object, ...){
  class(object) <- "clme"
  ranef(object, ...)
}


#' Extract random effects
#'
#' @description
#' Extracts the random effects estimates from objects of class \code{clme}.
#' 
#' @rdname ranef.clme
#' 
#' @return
#' Returns a numeric vector.
#' 
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#' 
#' ranef( clme.out )
#' 
#' 
#' @importFrom nlme ranef
#' @method ranef clme
#' @export
#' 
ranef.clme <- function( object, ... ){
  ## Print out the random effects
  if( is.clme(object) ){
    return( object$random.effects )
  } else{
    stop("'object' is not of class clme")
  }
}

#' Extract random effects
#' 
#' @param object object of class clme.
#' @param ... space for additional arguments
#'  
#' @rdname ranef
#' @importFrom nlme random.effects
#' @export
#' 
random.effects <- function( object, ... ){ UseMethod("random.effects") }

#' Extract random effects
#' 
#' @rdname ranef
#' @importFrom nlme random.effects
#' @export
#' 
random.effects.summary.clme <- function( object, ...){
  class(object) <- "clme"
  ranef(object, ...)
}


#' 
#' @rdname ranef.clme
#' @importFrom nlme random.effects
#' @method random.effects clme
#' @export
#' 
random.effects.clme <- function( object , ... ){
  ranef.clme( object, ... )
}



#' Various types of residuals 
#'
#' @description
#' Computes several types of residuals for objects of class \code{clme}. 
#' 
#' @param object object of class \code{\link{clme}}.
#' @param type type of residual (for mixed-effects models only).
#' @param ... space for additional arguments
#' 
#' @details
#' For fixed-effects models \eqn{Y = X\beta + \epsilon}{Y = X*b + e}, residuals are given as \deqn{\hat{e} = Y - X\hat{\beta}}{ ehat = Y - X*betahat}.
#' For mixed-effects models \eqn{Y = X\beta + + U\xi + \epsilon}{Y = X*b + U*xi + e}, three types of residuals are available.
#' \eqn{PA = Y - X\hat{\beta}}{ PA = Y - X*betahat}\\
#' \eqn{SS = U\hat{\xi}}{ SS = U*xihat}\\
#' \eqn{FM = Y - X\hat{\beta} - U\hat{\xi}}{ FM = Y - X*betahat - U*xihat}
#' 
#' @return
#' Returns a numeric matrix.
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' \dontrun{
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#'                  
#' residuals( clme.out, type='PA' )
#' }
#' @method residuals clme
#' @export
#' 
residuals.clme <- function( object, type="FM", ... ){
  ## Print out residuals of specified type
  if( is.clme(object) ){
    ridx <- which( c("PA", "SS", "FM")==type )
    if( ncol(object$residuals)<ridx ) ridx <- 1
    return( object$residuals[,ridx] )
  } else{
    stop("'object' is not of class clme")
  }
}

#' Various types of residuals 
#'
#' @rdname residuals.clme
#' 
#' @method residuals summary.clme
#' @export
#' 
residuals.summary.clme <- function( object, type="FM", ... ){
  class(object) <- "clme"
  residuals( object, type, ...)
}



#' Residual variance components
#' 
#' @param object object of class \code{\link{clme}}.
#' @param ... space for additional arguments
#' 
#' @description
#' Extract residual variance components for objects of class \code{clme}.
#'
#' @export
#' 
sigma <- function( object, ...){
  UseMethod("sigma")
}



#' Residual variance components
#'
#' @description
#' Extract residual variance components for objects of class \code{clme}.
#' 
#' @param object object of class \code{\link{clme}}.
#' @param ... space for additional arguments
#' 
#' 
#' @return
#' Numeric.
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#'                  
#' sigma( clme.out )
#' 
#' @method sigma clme
#' @export
#' 
sigma.clme <- function( object, ...){
  return( object$ssq )
}

#' Residual variance components
#'
#' @description
#' Extract residual variance components for objects of class \code{clme}.
#' 
#' @param object object of class \code{\link{clme}}.
#' @param ... space for additional arguments
#' 
#' 
#' @return
#' Numeric.
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#'                  
#' sigma( clme.out )
#' 
#' @method sigma summary.clme
#' @export
#' 
sigma.summary.clme <- function( object, ...){
  return( object$ssq )
}



#' Variance components
#' 
#' @param object object from which to extract variance components
#' @param ... space for additional arguments.
#' 
#' @rdname VarCorr
#' @export
#' 
VarCorr <- function( object, ...){ UseMethod("VarCorr") }


#' Variance components
#' 
#' @rdname VarCorr
#' @export
#' 
VarCorr.summary.clme <- function( object, ...){
  class(object) <- "clme"
  VarCorr(object, ...)
}


#' Variance components.
#'
#' @description
#' Extracts variance components for objects of class \code{clme}. 
#' 
#' @rdname VarCorr
#' 
#' @return
#' Numeric.
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#'                  
#' VarCorr( clme.out )
#' 
#' 
#' @importFrom lme4 VarCorr
#' @method VarCorr clme
#' @export
#' 
VarCorr.clme <- function(object, ...){
  ## Print out variances or SDs
  ## Defines tiny class "varcorr_clme" to handle printing
  ## using the method: print.varcorr_clme
  if( !is.clme(object) ){
    stop("'object' is not of class clme")
  } else{
    varcomps <- matrix( c(object$tsq, object$ssq ), ncol=1 )
    rnames   <- c( "Source", names(object$tsq), names(object$ssq) )
    rownames(varcomps) <- rnames[-1]
    colnames(varcomps) <- "Variance"
    class(varcomps)    <- "varcorr_clme"
    return( varcomps )
  }
}

## Leave this method out of alphabetical order so that
## is it right next to the VarCorr.clme method

#' Printout for variance components
#'
#' @description
#' Prints variance components of an objects of \code{clme}. 
#' 
#' @param object object of class \code{\link{clme}}.
#' @param rdig number of digits to round to.
#' @param ... space for additional arguments.
#' 
#' @return
#' Text printed to console.
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' \dontrun{
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#'                  
#' print.varcorr_clme( clme.out )
#' }
#' @importFrom stringr str_pad
#' @exportMethod print varcorr_clme
#' 
print.varcorr_clme <- function(object, rdig=5, ...){

  rnames   <- c( "Source", rownames( object ) )
  rnames   <- str_pad(rnames, width=max(nchar(rnames)), side = "right", pad = " ")
  vars     <- format( object , digits=rdig )
    
  cat( rnames[1], "\t" , "Variance" )
  for( ii in 1:length(vars) ){
    cat( "\n", rnames[ii+1], "\t" , vars[ii] )
  }
  
}



#' Variance-covariance matrix
#'
#' @description
#' Extracts variance-covariance matrix for objects of class \code{clme}.
#' 
#' @param object object of class \code{\link{clme}}.
#' @param ... space for additional arguments
#' 
#' 
#' @return
#' Numeric matrix.
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' 
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#'                  
#' vcov( clme.out )
#' 
#' @method vcov clme
#' @export
#' 
vcov.clme <- function(object, ...){
  ## Print out covariance matrix of theta
  if( is.clme(object) ){
    return( object$cov.theta )
  } else{
    stop("'object' is not of class clme")
  }  
}

#' Variance-covariance matrix
#'
#' @rdname vcov.clme
#' 
#' @method vcov summary.clme
#' @export
#' 
vcov.summary.clme <- function(object, ...){
  class(object) <- "clme"
  vcov(object, ...)
}


##
## Hidden functions to format characters / decimals
##



## Align the length of header with length of values
.align_table.clme <- function( tble, digits=4, ... ){
  cnames <- colnames( tble )
  for( ii in 1:length(cnames) ){
    ntitle <- nchar( cnames[ii] )
    maxc   <- max( nchar(tble[,ii]) )    
    if( ntitle > maxc ){
      tble[,ii] <- str_pad( tble[,ii], width=ntitle, side = "left", pad = " ")    
    }
    if( ntitle < maxc ){
      cnames[ii] <- str_pad( cnames[ii], width=(maxc+1), side = "right", pad = " ")    
    }
  }
  colnames(tble) <- cnames
  return( tble )
}





