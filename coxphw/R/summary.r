summary.coxphw <- function
(
  object,              # object of class coxphf
  print = TRUE,
  ...                  
)
  ### DD: 2014-01
  ### MP and GH: 2007-07
{
  if (print) { print(object$call) }
  if (object$template != "none")  { 
    if (object$template == "AHR")           { method <- "weighted estimation (AHR template)" } else
    if (object$template == "ARE")           { method <- "weighted estimation (ARE template)" } else
    if (object$template == "PH")            { method <- "unweighted estimation (PH template)" }
    if (print) { cat("\nModel fitted by", method, "\n\n") }
  }  
  
  se<-diag(object$var)^0.5
  out <- cbind(object$coefficients, se, exp(object$coefficients), object$ci.lower, object$ci.upper, object$coefficients/se, object$prob)
  dimnames(out) <- list(names(object$coefficients), c("coef", "se(coef)", "exp(coef)", paste(c("lower", "upper"), 1 - object$alpha), "z", "p"))
    
  if (!is.null(object$betafix)) { out[!is.na(object$betafix), -c(1,3)] <- NA }   
    
  if (print) { print(out) }
  
  if("loglik" %in% names(object)) {
    LL <- 2 * diff(object$loglik)
    if (print) { 
      cat("\nLikelihood ratio test=", LL, " on ", object$df, " df, p=", 1 - pchisq(LL, object$df), ", n=", object$n, "\n", sep = "")
      cat("\nScore test=", object$score, " on ", object$df, " df, p=", 1 - pchisq(object$score, object$df), ", n=", object$n, "\n", sep = "")
    }
  }
  
  if (print) { 
    if ( is.null(object$betafix)) { 
      cat("\nWald Chi-square =", object$Wald, "on", object$df, " df, p =", 1 - pchisq(object$Wald, object$df))
      cat("\n\nCovariance-Matrix:\n")
      print(object$var)
    } else
    if (!is.null(object$betafix)) { 
      w <- wald(coeff=object$coefficients[is.na(object$betafix)], cov=object$var[is.na(object$betafix), is.na(object$betafix)])
      cat("\nWald Chi-square =", w[1], "on", w[2], " df, p =", w[3], " (based on:", names(object$coefficients[is.na(object$betafix)]), ")")
      cat("\n\nCovariance-Matrix:\n")
      print(object$var[is.na(object$betafix), is.na(object$betafix)])
    }
    
    cat("\nGeneralized concordance probability:")
    if (max(object$template %in% c("PH", "ARE")) == 1) { cat("   Estimates may be biased!\n") } else { cat("\n") }
    cc <- concord(object)
    if (!is.null(object$betafix)) { cc[!is.na(object$betafix), -1] <- NA }
    print(cc)
  }

  invisible(object)
}
