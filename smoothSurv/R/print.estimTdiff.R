#####################################################
#### AUTHOR:    Arnost Komarek                   ####
####            25/02/2004                       ####
###             03/05/2004                       ####
###             28/06/2014  minor change         ####
####                                             ####
#### FILE:      print.estimTdiff.R               ####
####                                             ####
#### FUNCTIONS: print.estimTdiff                 ####
#####################################################

### ====================================================================
### print.estimTdiff: Print objects of class 'estimTdiff'
### ====================================================================
## x .......... object of class 'estimTdiff'
## digits ..... # of printed digits
## ... ........ other arguments passed to 'print' function
print.estimTdiff <- function(x, digits = min(options()$digits, 4), ...)
{

    if(is.null(digits))
        digits <- min(options()$digits, 4)

    if (is.null(attr(x, "cov1"))){
      #cat("   Only intercept was in the model.\n")
    }
    else{
      cat("\nCovariate Values Compared:\n")      
      cat("   Covariate values for T1:\n")
      print(attr(x, "cov1"), digits = digits, ...)
      cat("\n")
      cat("   Covariate values for T2:\n")
      print(attr(x, "cov2"), digits = digits, ...)            
    }
    cat("\n")

    if (!is.null(attr(x, "logscale.cov1"))){
      cat("\nLog-Scale Covariate Values Compared:\n")
      cat("   Log-scale covariate values for T1:\n")
      print(attr(x, "logscale.cov1"), digits = digits, ...)
      cat("\n")
      cat("   Log-scale covariate values for T2:\n")
      print(attr(x, "logscale.cov2"), digits = digits, ...)
      cat("\n")
    }         
    
    if (is.null(attr(x, "cov1"))){
      cat("\nEstimate of Expectation\n")
      cat("and Its ", attr(x, "conf.level")*100, "% Confidence Intervals\n\n", sep="")
      
      Z1 <- x$ET1 / x$sd.ET1
      p1 <- 2 * pnorm(-abs(Z1))
      show1 <- c(x$ET1, x$sd.ET1, x$ET1.lower, x$ET1.upper, Z1, p1)
      names(show1) <- c("ET", "Std.Error", "Lower", "Upper", "Z", "p")
      print(show1, digits = digits, ...); cat("\n")            
    }else{
      cat("\nEstimates of Expectations\n")
      cat("and Their ", attr(x, "conf.level")*100, "% Confidence Intervals\n", sep="")
      
      Z1 <- x$ET1 / x$sd.ET1
      Z2 <- x$ET2 / x$sd.ET2
      Zdiff <- x$diffT / x$sd.diffT

      p1 <- 2 * pnorm(-abs(Z1))
      p2 <- 2 * pnorm(-abs(Z2))
      pdiff <- 2 * pnorm(-abs(Zdiff))

      show1 <- data.frame(x$ET1, x$sd.ET1, x$ET1.lower, x$ET1.upper, Z1, p1)
      show2 <- data.frame(x$ET2, x$sd.ET2, x$ET2.lower, x$ET2.upper, Z2, p2)
      show3 <- data.frame(x$diffT, x$sd.diffT, x$diffT.lower, x$diffT.upper, Zdiff, pdiff)        
      
      colnames(show1) <- c("ET1", "Std.Error", "Lower", "Upper", "Z", "p")
      colnames(show2) <- c("ET2", "Std.Error", "Lower", "Upper", "Z", "p")
      colnames(show3) <- c("E(T1 - T2)", "Std.Error", "Lower", "Upper", "Z", "p")
      if (is.null(attr(x, "cov1"))) rncov1 <- paste("Value ", 1:length(x$ET1), sep = "") else rncov1 <- rownames(attr(x, "cov1"))
      if (is.null(attr(x, "cov2"))) rncov2 <- paste("Value ", 1:length(x$ET2), sep = "") else rncov2 <- rownames(attr(x, "cov2"))
      rownames(show1) <- rncov1
      rownames(show2) <- rncov2     
      rownames(show3) <- paste(rncov1, " - ", rncov2, sep = "")
      print(show1, digits = digits, ...); cat("\n")
      print(show2, digits = digits, ...); cat("\n")
      print(show3, digits = digits, ...); cat("\n")
    }

    return(invisible(x))
}

