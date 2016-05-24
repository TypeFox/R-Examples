### This files contains functions summary() and print.summary() for commonly
### used objects created in other files.
### Written: Wei-Chen Chen on 2008/10/14.

summary.emret <- function(object, ...){
  ret <- object
  ret$mean <- t(ret$Mu)
  ret$variance <- LTSigma2variance(ret$LTSigma)

  ret$tp <- length(ret$pi) - 1 + length(ret$Mu) + length(ret$LTSigma)
  if(! is.null(ret$llhdval)){
    ret$logL <- ret$llhdval
    ret$AIC <- -2 * ret$logL + 2 * ret$tp
    ret$BIC <- -2 * ret$logL + ret$tp * log(ret$n)
  }

  class(ret) <- "summary.emret"
  ret
}

print.summary.emret <- function(x,
    digits = max(4, getOption("digits") - 3), ...){
  emret <- x
  my.cat("Method: ", emret$method,
         "\n")
  my.cat(" n = ", emret$n,
         ", p = ", emret$p,
         ", nclass = ", emret$nclass,
         ", flag = ", emret$flag,
         ", total parameters = ", emret$tp,
         ",\n")
  if(!is.null(emret$logL)){
    my.cat(" logL = ", my.format(emret$logL, digits = digits),
           ", AIC = ", my.format(emret$AIC, digits = digits),
           ", BIC = ", my.format(emret$BIC, digits = digits),
           ".\n")
  }
  my.cat("nc: \n")
  my.print(emret$nc, digits = digits)
  my.cat("pi: \n")
  my.print(emret$pi, digits = digits)
}

print.emret <- function(x,
    digits = max(4, getOption("digits") - 3), ...){
  emret <- x
  my.cat("Method: ", emret$method,
         "\n")
  my.cat(" n = ", emret$n,
         ", p = ", emret$p,
         ", nclass = ", emret$nclass,
         ", flag = ", emret$flag,
         ", logL = ", my.format(emret$llhdval, digits = digits),
         ".\n")
  my.cat("nc: \n")
  my.print(emret$nc, digits = digits)
  my.cat("pi: \n")
  my.print(emret$pi, digits = digits)
}



### This function is only for advance users.
summary.emret.wt <- function(object, ...){
  ret <- wt2wot(object)
  ret$mean <- t(ret$Mu)
  ret$variance <- LTSigma2variance(ret$LTSigma)

  ret$tp <- length(ret$pi) - 1 + length(ret$Mu) + length(ret$LTSigma)
  if(! is.null(ret$llhdval)){
    ret$logL <- ret$llhdval
    ret$AIC <- -2 * ret$logL + 2 * ret$tp
    ret$BIC <- -2 * ret$logL + ret$tp * log(ret$n)
  }

  class(ret) <- "summary.emret"
  ret
}
