
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   equal.kappa.test function                               #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: November, 06, 2013                                #
#   Version: 0.4                                            #
#                                                           #
#   Copyright (C) 2013 Claudio Agostinelli                  #
#                                                           #
#############################################################

equal.kappa.test <- function(x, group) {
    # Handling missing values
    ok <- complete.cases(x, group)
    x <- x[ok]
    group <- group[ok, drop = TRUE]
    if (length(x)==0 | length(table(group)) < 2) {
        warning("No observations or no groups (at least after removing missing values)")
        return(NULL)
    }
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(x, "class") <- attr(x, "circularp") <- NULL
    result <- EqualKappaTestRad(x, group)
    result$call <- match.call()
    class(result) <- "equal.kappa.test"
    return(result)
}

EqualKappaTestRad <- function(x, group) {
   x <- x%%(2*pi)
   ns      <- tapply(x, group, FUN=length)
   r.bars  <- tapply(x, group, FUN=RhoCircularRad)
   rs      <- r.bars*ns
   kappas  <- tapply(x, group, FUN=function(x) MlevonmisesRad(x)[4])
   grps    <- length(r.bars)
   n       <- length(group)
   r.bar.all  <- RhoCircularRad(x)
   kappa.all  <- MlevonmisesRad(x)[4]
   warn1 <- 0
    
   if (r.bar.all < 0.45){
      g1 <- function(x){asin(sqrt(3/8)*x)}
      ws <- 4*(ns-4)/3
      g1s <- g1(2*r.bars)
      U  <- sum(ws*g1s^2) - sum(ws*g1s)^2/sum(ws)
      if (any(is.na(g1s))) {
         warn1 <- 1
         warning("An argument outside of [-1,1] was passed to asin function in calculation of approximate chi-squared test statistic. Bartlett's test of homogeneity was used instead of the approximation using asin.")
      }
   }
   if (r.bar.all >= 0.45 & r.bar.all <= 0.70){
      g2 <- function(x){
          c1 <- 1.089
          c2 <- 0.258
          asinh((x-c1)/c2)
      }
      ws <- (ns-3)/0.798
      g2s <- g2(r.bars)
      U  <- sum(ws*g2s^2) - sum(ws*g2s)^2/sum(ws)
   }    
   if (r.bar.all > 0.70 | warn1==1){
      vs <- ns-1
      v  <- n-grps
      d  <- 1/(3*(grps-1))*(sum(1/vs)-1/v)
      U  <- 1/(1+d)*(v*log((n-sum(rs))/v) - sum(vs*log((ns-rs)/vs)))
   }
   p.value <- 1-pchisq(U, grps-1)
   result <- list(kappa=kappas, kappa.all=kappa.all, rho=r.bars, rho.all=r.bar.all, df=grps-1, statistic=U, p.value=p.value)
   return(result)
}

#############################################################
#                                                           #
#   print.equal.kappa.test functio                          #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: April, 13, 2005                                   #
#   Version: 0.1-1                                          #
#                                                           #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.equal.kappa.test <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")    
    cat("\n", "Test for Homogeneity of Concentration Parameters", "\n \n")
    cat(" df:     ", format(x$df, digits=digits), "\n ChiSq:  ", format(x$statistic, digits=digits), "\n p.value:",  format(x$p.value, digits=digits), "\n \n") 
    invisible(x)
}

