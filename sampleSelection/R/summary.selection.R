summary.selection <- function(object, ...) {
   ## object      object of class "selection"

   if( object$method == "ml" ) {
      s <- NextMethod( "summary", object, ...)
   } else if( object$method == "2step" )  {
      s <- list()  # list for results that will be returned
      RSq <- function(model, intercept) {
         ## Calculate r-squared.  Note that the way lm() finds R2 is
         ## a bit naive -- it checks for intercept
         ## in the formula, but not whether the intercept is present
         ## in any of the data vectors (or matrices)
         if(class(model) == "lm") {
            y <- model.response(model.frame(model))
            if(intercept) {
               R2 <- 1 - sum(residuals(model)^2)/sum((y - mean(y))^2)
               R2adj <- 1 - (1 - R2)*(nObs(model) - 1)/(nObs(model) - nParam(model))
            }
            else {
               R2 <- 1 - sum(residuals(model)^2)/sum(y^2)
               R2adj <- 1 - (1 - R2)*(nObs(model))/(nObs(model) - nParam(model))
            }
         }
         else {
            R2 <- summary( model$eq[[ 1 ]] )$r.squared
            R2adj <- summary( model$eq[[ 1 ]] )$adj.r.squared
         }
         c(R2, R2adj)
      }
      if( object$tobitType == 2 ) {
         r <- RSq( object$lm, object$param$oIntercept )
         R2 <- r[1]
         R2adj <- r[2]
         s$rSquared <- list(R2=R2, R2adj=R2adj)
      } else {
         r <- RSq(object$lm1, object$param$oIntercept1)
         R21 <- r[1]
         R2adj1 <- r[2]
         r <- RSq(object$lm2, object$param$oIntercept2)
         R22 <- r[1]
         R2adj2 <- r[2]
         iBetaS <- object$param$index$betaS
         s$rSquared <- list(R21=R21, R2adj1=R2adj1, R22=R22, R2adj2=R2adj2)
      }
      stdd <- sqrt(diag(vcov(object, part="full")))
      estcoef <- coef( object, part="full" )
      names( estcoef ) <- sub( "^[SO][12]?:", "", names( estcoef ) )
      s$estimate <- coefTable( estcoef, stdd, object$param$df)
   }
   s$param    <- object$param
   s$tobitType <- object$tobitType
   s$method <- object$method
   s$activePar <- activePar(object)
   class( s ) <- c( "summary.selection", class( s ) )
   return( s )
}

print.summary.selection <- function(x,
                                 digits=max(3, getOption("digits") - 3),
                                 part="full",
                                 ... ) {

   cat("--------------------------------------------\n")
   cat("Tobit", x$tobitType, "model" )
   if( x$tobitType == 2 ) {
      cat( " (sample selection model)\n" )
   } else {
      cat( " (switching regression model)\n" )
   }
   if( x$method == "ml" ) {
      cat( "Maximum Likelihood estimation\n" )
      cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
      cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
      if(!is.null(x$estimate)) {
         cat("Log-Likelihood:", logLik(x), "\n")
      }
   } else if( x$method == "2step" ) {
      cat( "2-step Heckman / heckit estimation\n" )
   }
   if(!is.null(x$estimate)) {
      cat( x$param$nObs, "observations" )
      if( x$tobitType == 2 ) {
         cat( " (", x$param$N0, " censored and ", x$param$N1, " observed)\n",
             sep = "" )
      }
      else if(x$tobitType == 5) {
         cat( ": ", x$param$N1, " selection 1 (", x$param$levels[1], ") and ",
             x$param$N2, " selection 2 (", x$param$levels[2], ")\n", sep = "" )
      }
      else if(x$tobitType == "treatment") {
         cat( ": ", x$param$N0, " non-participants (selection ",
             x$param$levels[1], ") and ",
             x$param$N1, " participants (selection ",
             x$param$levels[2], ")\n", sep = "", fill=TRUE)
      }
      else
         stop("Tobit type must be either '2', '5', or 'treatment'")
      cat(sum(x$activePar), "free parameters" )
      cat( " (df = ", x$param$df, ")\n", sep="")
      if(part == "full") {
         cat("Probit selection equation:\n")
         printCoefmat( x$estimate[ x$param$index$betaS,,drop=FALSE],
            signif.legend = FALSE, digits = digits )
      }
      if( x$tobitType == 2 ) {
         cat("Outcome equation:\n")
         printCoefmat( x$estimate[ x$param$index$betaO,,drop=FALSE],
            signif.legend = ( part != "full" ), digits = digits )
         if( x$method == "2step" ) {
            cat("Multiple R-Squared:", round(x$rSquared$R2, digits),
               ",\tAdjusted R-Squared:", round(x$rSquared$R2adj, digits),
               "\n", sep="")
         }
      } else if( x$tobitType == 5 ) {
         cat("Outcome equation 1:\n")
         printCoefmat( x$estimate[ x$param$index$betaO1,,drop=FALSE],
            signif.legend = FALSE, digits = digits )
         if( x$method == "2step" ) {
            cat("Multiple R-Squared:", round(x$rSquared$R21, digits),
               ",\tAdjusted R-Squared:", round(x$rSquared$R2adj1, digits),
               "\n", sep="")
         }
         cat("Outcome equation 2:\n")
         printCoefmat( x$estimate[ x$param$index$betaO2,,drop=FALSE],
            signif.legend = ( part != "full" ), digits = digits )
         if( x$method == "2step" ) {
            cat("Multiple R-Squared:", round(x$rSquared$R22, digits),
               ",\tAdjusted R-Squared:", round(x$rSquared$R2adj2, digits),
            "\n", sep="")
         }
      }
      else if(x$tobitType == "treatment") {
         cat("Outcome equation:\n")
         printCoefmat( x$estimate[ x$param$index$betaO,,drop=FALSE],
            signif.legend = ( part != "full" ), digits = digits )
         if( x$method == "2step" ) {
            cat("Multiple R-Squared:", round(x$rSquared$R2, digits),
               ",\tAdjusted R-Squared:", round(x$rSquared$R2adj, digits),
               "\n", sep="")
         }
      }
      if(part=="full") {
         cat("   Error terms:\n")
         printCoefmat( x$estimate[ x$param$index$errTerms,,drop=FALSE],
            digits = digits )
      }
   }
   cat("--------------------------------------------\n")
   invisible( x )
}
