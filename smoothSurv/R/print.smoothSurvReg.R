###########################################
#### AUTHOR:    Arnost Komarek         ####
####            01/05/2004             ####
####                                   ####
#### FILE:      print.smoothSurvReg.R  ####
####                                   ####
#### FUNCTIONS: print.smoothSurvReg    ####
###########################################

### ====================================================================
### print.smoothSurvReg: Print objects of class 'smoothSurvReg'
### ====================================================================
## x .......... object of class 'smoothSurvReg'
## spline ..... T/F, do I want to print an information concerning the fitted spline?
## digits ..... # of printed digits
## ... ........ other arguments passed to 'print' function
print.smoothSurvReg <- function(x, spline, digits = min(options()$digits, 4), ...)
{
    if (x$fail >= 99) {
        cat("No summary, smoothSurvReg failed.\n")
        return(invisible(x))
    }

    if(missing(digits))
        digits <- min(options()$digits, 4)

    if(missing(spline)) spline <- (nrow(x$spline) <= 31)

    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }

    est.scale <- x$estimated["Scale"]
    est.c <- x$estimated["ccoef"]

    ## Estimates of regres parameters
    if (!is.null(x$regres)){
       nregres <- dim(x$regres)[1]
       cat("\nEstimated Regression Coefficients:\n")
#       cat("-----------------------------\n")
       scale <- NULL
       if ("Scale" %in% rownames(x$regres)){
	 regres <- x$regres[1:(nregres-1),]
	 scale <- x$regres[nregres, "Value"]
       }
       else{
         regres <- x$regres
       }
       Z.P <- regres$Value/regres[["Std.Error"]]
       pv.P <- 2 * pnorm(-abs(Z.P))
       Z.V <- regres$Value/regres[["Std.Error2"]]
       pv.V <- 2 * pnorm(-abs(Z.V))
       regres[["Z"]] <- Z.P
       regres[["Z2"]] <- Z.V
       regres[["p"]] <- pv.P
       regres[["p2"]] <- pv.V
       print(regres, digits = digits, ...)
       if (!is.null(scale)){
          cat("\nScale =", format(scale, digits=digits), "\n")
       }
    }

    ## Fixed regression components of the model
    regres.fixed <- NULL
    tfr.names <- character(0)
    if (!est.scale){
      tfr.names <- c(tfr.names, "Log(scale)", "Scale")
      regres.fixed <- c(regres.fixed, x$init.regres["Log(scale)", "Value"], x$init.regres["Scale", "Value"])
    }

    if (!is.null(regres.fixed)){
      regres.fixed <- data.frame(Value = regres.fixed)
      rownames(regres.fixed) <- tfr.names
      nregres <- dim(regres.fixed)[1]

      cat("\nFixed Regression Coefficients:\n")
#      cat("-------------------------\n")
      scale <- NULL
      if ("Scale" %in% rownames(regres.fixed)){
        regres <- as.data.frame(regres.fixed[1:(nregres-1),])
	colnames(regres) <- "Value"
	scale <- regres.fixed[nregres, "Value"]
      }
      else{
        regres <- regres.fixed
      }
      print(regres, digits = digits, ...)
      if (!is.null(scale)){
         cat("\nScale:", format(scale, digits=digits), "\n")
      }
    }

    ## Possibly adjusted intercept and scale
#    cat("\nAdjusted Intercept and Scale:\n")
#    cat("-----------------------------\n")
#    print(x$adjust, digits = digits, ...)

    ## Fitted error distribution
#    cat("\n(Fitted) Error Distribution:\n")
#    cat("----------------------------\n")
#    m.err <- x$error.dist$Mean
#    sd.err <- x$error.dist$SD
#    cat("Mean:", format(m.err, digits=digits), "\n")
#    cat("Scale:", format(sd.err, digits=digits), "\n")

    if(spline){
       Z.P <- x$spline[["c coef."]]/x$spline[["Std.Error.c"]]
       pv.P <- 2 * pnorm(-abs(Z.P))
       Z.V <- x$spline[["c coef."]]/x$spline[["Std.Error2.c"]]
       pv.V <- 2 * pnorm(-abs(Z.V))
       spline.print <- x$spline[,1:5]     ## do not print columns with a's
       spline.print[["Z"]] <- Z.P
       spline.print[["Z2"]] <- Z.V
       spline.print[["p"]] <- pv.P
       spline.print[["p2"]] <- pv.V
       cat("\nDetails on (Fitted) Error Distribution:\n")
#       cat("---------------------------------------\n")
       print(spline.print, digits = digits, ...)
    }

    ## Likelihood and iterations
    digits <- digits+3
    cat("\nPenalized Loglikelihood and Its Components:\n")
#    cat("----------------------------------\n")
    pll <- x$loglik[1, "Penalized Log Likelihood"]
    ll <- x$loglik[1, "Log Likelihood"]
    penalty <- x$loglik[1, "Penalty"]
    df <- x$degree.smooth$df
    nparam <- x$degree.smooth[["Number of parameters"]]
    nU <- x$degree.smooth[["Mean param."]]
    nUscale <- x$degree.smooth[["Scale param."]]
    nUc <- x$degree.smooth[["Spline param."]]
    lambda <- x$degree.smooth$Lambda
    loglambda <- x$degree.smooth[, "Log(Lambda)"]
    aic <- x$aic
    n <- dim(x$y)[1]
    nmiss <- length(x$na.action)
    cat("     Log-likelihood:", format(ll, digits=digits), "\n")
    cat("            Penalty:", format(penalty, digits=digits), "\n")

    cat("   Penalized Log-likelihood:", format(pll, digits=digits), "\n\n")

    cat("Degree of smoothing:\n")
    cat("   Number of parameters:", format(nparam, digits=digits), "\n")
    cat("                   Mean parameters:", format(nU, digits=digits), "\n")
    cat("                  Scale parameters:", format(nUscale, digits=digits), "\n")
    cat("                 Spline parameters:", format(nUc, digits=digits), "\n\n")

    cat("                   Lambda:", format(lambda, digits=digits), "\n")
    cat("              Log(Lambda):", format(loglambda, digits=digits), "\n")
    cat("                       df:", format(df, digits=digits), "\n\n")

    cat("AIC (higher is better): ", format(aic, digits=digits), "\n\n")

    cat("Number of Newton-Raphson Iterations: ", x$iter, "\n")
    if (nmiss > 0) cat("n = ", n, " (", nmiss, "observations deleted due to missing)", "\n", sep = "")
    else           cat("n =", n, "\n")
}

