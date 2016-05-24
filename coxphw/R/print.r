print.coxphw <- function
(
  x,             # object of class coxphw
  ...            # dummy
)
  ### DD 2014-01
  ### MP and GH 2007-07
  ### overall test is score test
{
  print(x$call)
  if (x$template != "none")  { 
    if (x$template == "AHR")           { method <- "weighted estimation (AHR template)" } else
    if (x$template == "ARE")           { method <- "weighted estimation (ARE template)" } else
    if (x$template == "PH")            { method <- "unweighted estimation (PH template)" }
    cat("\nModel fitted by", method, "\n\n")  
  }  

  se<- diag(x$var)^0.5
  out <- cbind(x$coefficients, se, exp(x$coefficients), x$ci.lower, x$ci.upper, x$coefficients/se, x$prob)
  dimnames(out) <- list(names(x$coefficients), c("coef", "se(coef)", "exp(coef)", paste(c("lower", "upper"), 1 - x$alpha), "z", "p"))
    
  if (!is.null(x$betafix)) { out[!is.na(x$betafix), -c(1,3)] <- NA }   
    
  print(out)
    
  if ( is.null(x$betafix)) { cat("\nWald Chi-square=", x$Wald, " on ", x$df, "df, p=", 1 - pchisq(x$Wald, x$df), ", n=", x$n, "\n\n", sep = "") } else
  if (!is.null(x$betafix)) { 
    w <- wald(coeff=x$coefficients[is.na(x$betafix)], cov=x$var[is.na(x$betafix), is.na(x$betafix)])
    cat("\nWald Chi-square =", w[1], "on", w[2], " df, p =", w[3], " (based on:", names(x$coefficients[is.na(x$betafix)]), ")")
  }
  
  
  invisible(x)
}
