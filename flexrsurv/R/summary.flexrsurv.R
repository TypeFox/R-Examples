
# version du 13.06.13

summary.flexrsurv <- function (object, 
                               correlation = FALSE,
                               symbolic.cor = FALSE, ...) {  
  
  # definition of fbis (used if alternative iterative algorithm)
  # ----------------------------------------------------------------------------
  f <- as.formula(object$formula)
  termsf <- terms(f, specials=c("NLL","NPH","NPHNLL"))  
  n <- length(attr(termsf,"order"))
  fbis <- vector("character",n)
  
  if (terms(f)[[3]]==1) {
    n <- 1 
    fbis <- 1
  }
  else {
    for (i in 1:n) {
      fbis[i] <- deparse(terms(f)[i][[3]],width.cutoff=500L)
    }
  }
      
  model <- object$call[["model"]]
  if (is.null(object$call[["model"]])) {model <- eval(formals(flexrsurv)$model)[1]}
  
  method <- object$call[["method"]]
  if (is.null(object$call[["method"]])) {method <- eval(formals(flexrsurv)$method)[1]}


  aliased <- is.na(coef(object))  # used in print method


  
  # 1- if method=="glm" with NPHNLL --> none variance matrix
  # ----------------------------------------------------------------------------

      if (inherits(object, "flexrsurv.glmiterative")){
        warning("Variance matrix is not available with NPHNLL effect and GLM method.")
      }
  
      if ( inherits(object, "glm") ){
       # modified call to summary.glm
       # in order to prevent computation of residuals
      if(object$df.residual>0){
        object$df.residual <-  -object$df.residual
      }  
      ans <- summary.glm(object, dispersion = 1,
                         correlation = correlation,
                         symbolic.cor = symbolic.cor, ...)
      ans$cov <- ans$cov.unscaled
      ans$loglik <- object$loglik
      ans$family <- NULL
      ans$null.deviance <- NULL
      ans$dispersion <- NULL
      ans$cov.scaled <- NULL
      ans$cov.unscaled <- NULL

      attr(ans, "fitclass") <- class(object)
            
      class(ans) <- "summary.flexrsurv"
      return(ans)

    }
    else {
           # maximum likelyhood estimation or glmiterative
      # buid output matrix of coef


      aliased <- is.na(coef(object))  # used in print method
      coef <- object$coefficients
      if( is.null(dim(object$var))){
        coef.table <- cbind(coef, NaN, NaN, NaN)
        covmat <- object$var
        dimnames(coef.table) <- list(names(coef),
                                     c("Estimate", "Std. Error", "z value","Pr(>|z|)"))
      }
      else {
        covmat <- object$var
        dimnames(covmat) <- list(names(coef),names(coef))
        var.coef <- diag(covmat)
        std.err <- sqrt(var.coef)
        tvalue <- coef/std.err
        pvalue <- 2*pnorm(-abs(tvalue))
        coef.table <- cbind(coef, std.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef),
                                     c("Estimate", "Std. Error", "z value","Pr(>|z|)"))
      }
      keep <- match(c("call","terms","loglik", 
                      "contrasts", "df.residual",
                      "iter", "na.action"), names(object), 0L)
      
      ans <- c(object[keep],
               list(coefficients = coef.table,
                    cov = covmat,
                    aliased = aliased
                    ))
      if(correlation > 0 & !is.null(object$var)) {
        dd <- sqrt(diag(covmat))
        ans$correlation <-  covmat/outer(dd,dd)
        symbolic.cor <- symbolic.cor
      }
    }
  
  ans$optim.control <- object$optim.control
  attr(ans, "fitclass") <- class(object)
  class(ans) <- "summary.flexrsurv"
  return(ans)
}

      
