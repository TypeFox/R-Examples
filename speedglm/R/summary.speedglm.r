summary.speedglm <-
  function(object,correlation=FALSE,...)
  {
    
    if (!inherits(object,"speedglm")) stop("object is not of class speedglm")
    z <- object
    var_res <- as.numeric(z$RSS/z$df)
    dispersion <- if (z$family$family %in% c("poisson","binomial")) 1 else var_res
    if(z$method=='qr') {
      z$XTX <- z$XTX[z$ok,z$ok]
    }
    inv <- solve(z$XTX,tol=z$tol.solve)
    covmat <- diag(inv)
    se_coef <- rep(NA,length(z$coefficients)) 
    se_coef[z$ok] <- sqrt(dispersion*covmat)   
    
    if (z$family$family %in% c("binomial","poisson"))  { 
      z1      <- z$coefficients/se_coef
      p       <- 2*pnorm(abs(z1),lower.tail=FALSE)    
    } else    {
      t1      <- z$coefficients/se_coef
      p       <- 2*pt(abs(t1),df=z$df,lower.tail=FALSE)       
    }
    ip <- !is.na(p)  
    p[ip] <- as.numeric(format(p[ip],digits=3))
    dn <- c("Estimate", "Std. Error")
    if (z$family$family %in% c("binomial","poisson")) 
    {
      format.coef <- if (any(na.omit(abs(z$coef))<0.0001)) format(z$coefficients,
                                                                  scientific=TRUE,digits=4) else round(z$coefficients,digits=7)
      format.se <- if (any(na.omit(se_coef)<0.0001)) format(se_coef,scientific=TRUE,
                                                            digits=4) else round(se_coef,digits=7)
      format.pv <- if (any(na.omit(p)<0.0001)) format(p,scientific=TRUE,digits=4) else 
        round(p,digits=4)                      
      param <- data.frame(format.coef,format.se, round(z1,digits=4),format.pv)
      
      dimnames(param) <- list(names(z$coefficients), c(dn,"z value", "Pr(>|z|)"))                            
      
    } else 
    {
      format.coef <- if (any(abs(na.omit(z$coefficients))<0.0001)) format(z$coefficients,
                                                                          scientific=TRUE,digits=4) else round(z$coefficients,digits=7)
      format.se <- if (any(na.omit(se_coef)<0.0001)) format(se_coef,
                                                            scientific=TRUE,digits=4) else round(se_coef,digits=7)
      format.pv <- if (any(na.omit(p)<0.0001)) format(p,scientific=TRUE,digits=4) else 
        round(p,digits=4)                      
      
      param <- data.frame(format.coef,format.se, round(t1,digits=4),format.pv)
      
      dimnames(param) <- list(names(z$coefficients), c(dn,"t value", "Pr(>|t|)"))                            
    }                
    eps <- 10 * .Machine$double.eps
    if (z$family$family == "binomial") {
      if (any(z$mu > 1 - eps) || any(z$mu < eps))
        warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (z$family$family == "poisson") {
      if (any(z$mu < eps))  warning("fitted rates numerically 0 occurred")
    }  
    
    keep <- match(c("call", "terms", "family", "deviance", "aic", 
                    "df", "nulldev", "nulldf", "iter","tol","n","convergence",
                    "ngoodobs","logLik","RSS","rank"), names(object), 0)
    
    ans <- c(object[keep], list(coefficients = param, dispersion = dispersion, 
                                correlation = correlation, cov.unscaled = inv, cov.scaled = inv * var_res ))
    
    if (correlation) {         
      ans$correl <- ( inv * var_res)/outer(na.omit(se_coef), na.omit(se_coef))
    }          
    
    class(ans) <- "summary.speedglm"
    return(ans)
  }

print.summary.speedglm <- function(x,digits = max(3, getOption("digits") - 3),...)
{
  cat("Generalized Linear Model of class 'speedglm':\n")
  if (!is.null(x$call)) cat("\nCall: ", deparse(x$call), "\n\n")
  if (length(x$coef)) {
    cat("Coefficients:\n")
  cat(" ------------------------------------------------------------------", "\n")
  sig <- function(z) {
    if (z < 0.001) 
      "***"
    else if (z < 0.01) 
      "** "
    else if (z < 0.05) 
      "*  "
    else if (z < 0.1) 
      ".  "
    else "   "
  }
  sig.1 <- sapply(as.numeric(as.character(x$coefficients$"Pr(>|t|)")), sig)
  est.1 <- cbind(format(x$coefficients, digits = digits), sig.1)
  colnames(est.1)[ncol(est.1)] <- ""
  print(est.1)
  cat("\n")
  cat("-------------------------------------------------------------------", 
      "\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", 
      "\n")
  cat("\n")  
  } else cat("No coefficients\n")
  cat("---\n")
  cat("null df: ",x$nulldf,"; null deviance: ",round(x$nulldev,digits=2),";\n",
      "residuals df: ",x$df,"; residuals deviance: ",round(x$deviance,
                                                           digits=2),";\n",
      "# obs.: ",x$n,"; # non-zero weighted obs.: ",x$ngoodobs,";\n",
      "AIC: ",x$aic,"; log Likelihood: ",x$logLik,";\n",
      "RSS: ",round(x$RSS,digits=1),"; dispersion: ",x$dispersion,"; iterations: ",
      x$iter,";\n","rank: ",round(x$rank,digits=1),"; max tolerance: ",
      format(x$tol,scientific=TRUE,digits=3),"; convergence: ",x$convergence,".\n",sep="")
  invisible(x)
  if (x$correlation) {
    cat("---\n")
    cat("Correlation of Coefficients:\n")
    x$correl[upper.tri(x$correl,diag=TRUE)]<-NA
    print(x$correl[-1,-nrow(x$correl)],na.print="",digits=2)
    
  }
}

print.speedglm <- function(x,digits = max(3, getOption("digits") - 3),...)
{
  cat("Generalized Linear Model of class 'speedglm':\n")
  if (!is.null(x$call)) cat("\nCall: ", deparse(x$call), "\n\n")
  if (length(x$coef)) {
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2,
                  quote = FALSE)
  } else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

coef.speedglm <- function(object,...) object$coefficients

vcov.speedglm <- function(object,...) object$dispersion*solve(object$XTX)

logLik.speedglm <- function(object,...) object$logLik

AIC.speedglm <- function(object,...) {
  if (!(length(list(...)))) object$aic else { 
    aic<-function(x) x$aic
    object <- list(object,...)
    val <- sapply(object, aic)
    val
  } 
}

extractAIC.speedglm<-function(fit, scale=0, k=2,...) 
{
  n <- fit$n
  edf <- n - fit$df
  aic <- fit$aic
  c(edf, aic +  (k - 2) * edf)
}



