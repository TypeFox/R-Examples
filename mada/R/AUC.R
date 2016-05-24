AUC <- function(x, ...) UseMethod("AUC")

AUC.default <- function(x, fpr = 1:99/100, ...){
  sroc <- x
  stopifnot(is.vector(fpr), all(fpr <= 1), all(fpr >= 0), is.function(sroc))
  if(! length(formals(sroc)) == 1){stop("expected a function with exactly one argument")}
  n <- length(fpr)
  if(n < 10){stop("specify at least 10 FPR values!")}
  s <- numeric(n)
  for(i in 1: n){
    temp <- try(sroc(fpr[i]))
    if(class(temp) == "try-error"){stop(paste("calculation of sroc failed for value of FPR ", fpr[i]))}
    if(temp < 0 | temp > 1){stop("expected values of sroc to be >= 0 and <= 1, but this is not the case for FPR value ", fpr[i])}      
    s[i] <- temp
  } # end of loop 
  AUC <-  (s[1]/2 + sum(s[2:(n-1)]) + s[n]/2)/n
  ret <- list(AUC = AUC)
  class(ret) <- "AUC"
  return(ret)
}

AUC.phm <- function(x, level = 0.95, ...){
  theta <- coef(x)[1]
  ci <- confint(x, level = level)["theta",]
  AUC <- 1/(theta+1)
  AUCci <- 1/(ci+1)
  obsfprrange <- range(fpr(x$data))
  pAUC <- function(theta,a,b){1/(theta+1)*(b^{theta+1}-a^{theta+1})/(b-a)}
  ret <- list(AUC = AUC, ci = AUCci, 
              pAUC = pAUC(theta,obsfprrange[1],obsfprrange[2]), 
              pci = c(pAUC(ci[1],obsfprrange[1],obsfprrange[2]),
                      pAUC(ci[2],obsfprrange[1],obsfprrange[2])))
  class(ret) <- "AUC"
  return(ret)
}

## AUC.reitsma in file reitsma.R
