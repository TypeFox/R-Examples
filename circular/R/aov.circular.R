
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   aov.circular function                                   #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: August, 10, 2006                                  #
#   Version: 0.2-1                                          #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

aov.circular <- function(x, group, kappa=NULL, method=c("F.test", "LRT"), F.mod=TRUE, control.circular=list()) {

    method <- match.arg(method)  
    # Handling missing values
    ok <- complete.cases(x, group)
    x <- x[ok]
    group <- group[ok]
    if (length(x)==0 | length(table(group)) < 2) {
        warning("No observations or no groups (at least after removing missing values)")
        return(NULL)
    }
    if (is.circular(x)) {
       datacircularp <- circularp(x)     
    } else {
       datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
    }

    dc <- control.circular
    if (is.null(dc$type))
       dc$type <- datacircularp$type
    if (is.null(dc$units))
       dc$units <- datacircularp$units
    if (is.null(dc$template))
       dc$template <- datacircularp$template
    if (is.null(dc$modulo))
       dc$modulo <- datacircularp$modulo
    if (is.null(dc$zero))
       dc$zero <- datacircularp$zero
    if (is.null(dc$rotation))
       dc$rotation <- datacircularp$rotation
    
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(x, "class") <- attr(x, "circularp") <- NULL

    result <- AovCircularRad(x, group, kappa=NULL, method, F.mod)
    result$mu <- conversion.circular(circular(result$mu), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    result$mu.all <- conversion.circular(circular(result$mu.all), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation) 
    result$call <- match.call()
    class(result) <- "aov.circular"
    return(result)
}

AovCircularRad <- function(x, group, kappa=NULL, method, F.mod) {
### x must be in radians, modulo 2pi
   ns        <- tapply(x, group, FUN=length)
   resultant <- tapply(x, group, FUN=function(x) RhoCircularRad(x)*length(x))
   mean.dirs <- tapply(x, group, FUN=MeanCircularRad)
   kappas    <- tapply(x, group, FUN=function(x) MlevonmisesRad(x)[4])
   grps <- length(resultant)
   n <- length(group)
   res.all <- RhoCircularRad(x)*n
   mean.dir.all <- MeanCircularRad(x)
   kappa.all <- MlevonmisesRad(x)[4]

   if (method=="F.test"){
      if (!is.null(kappa))
         warning("Specified value of kappas is not used in the F-test")
      sum.res <- sum(resultant)
      df <- c(grps-1, n-grps, n-1)
      SS <- c(sum.res - res.all, n-sum.res, n-res.all) 
      MS <- SS/df
      if (F.mod==TRUE) {
         stat <- (1+3/(8*kappa.all))*MS[1]/MS[2]
      } else {
         stat <- MS[1]/MS[2]
      }
      p.value <- 1-pf(stat, grps-1,n-grps)
   } else {
      SS <- NA
      MS <- NA
      if (is.null(kappa))
         kappa <- kappa.all
      stat1 <- 1-1/(4*kappa)*A1(kappa)*(sum(1/ns)-1/n)
      stat2 <- 2*kappa*sum(resultant*(1-cos(mean.dirs-mean.dir.all)))
      stat <- stat1*stat2
      df <- grps-1
      p.value <- 1-pchisq(stat, df)
   }
   result <- list()
   result$mu <- mean.dirs 
   result$mu.all <- mean.dir.all
   result$kappa <- kappas 
   result$kappa.all <- kappa.all
   result$rho <- resultant
   result$rho.all <- res.all
   result$method <- method
   result$df <- df 
   result$SS <- SS
   result$MS <- MS
   result$statistic <- stat
   result$p.value <- p.value
   return(result)
}

#############################################################
#                                                           #
#   print.aov.circular function                             #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: April, 30, 2005                                   #
#   Version: 0.2-1                                          #
#                                                           #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.aov.circular <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="") 
 
    if (x$method=="F.test") {
        result.matrix <- cbind(x$df, x$SS, x$MS, c(x$statistic,NA,NA), c(x$p.value,NA,NA))
        dimnames(result.matrix) <- list(c("Between","Within","Total"),c("df", "SS", "MS", "F", "p"))
        cat("\n", "Circular Analysis of Variance: High Concentration F-Test", "\n", "\n")
        print(result.matrix, digits=digits)
        cat("\n \n")
    } else {
        cat("\n", "Circular Analysis of Variance: Likelihood Ratio Test", "\n", "\n")
        cat(" df:     ", format(x$df, digits=digits), "\n ChiSq:  ", format(x$statistic, digits=digits), "\n p.value:",  format(x$p.value, digits=digits), "\n \n")
    }
    invisible(x)
}

