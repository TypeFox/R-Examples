summary.DirichletRegModel <- function(object, ...) {

  if(object$optimization$convergence > 2L) stop("\n", paste(strwrap(paste("\nOptimization did not converge in",object$optimization$bfgs.it,"+",object$optimization$iterations,"iterations and exited with code",object$optimization$convergence), getOption("width")),sep="\n",collapse="\n"))

  res <- structure(list(
    call = object[["call"]],
    terms = terms(object$formula),
    logLik = object$logLik,
    df = object$npar,
    deviance = -2.0 * object$logLik,
    aic = AIC(object),
    bic = BIC(object),
    residuals = residuals(object, type="standardized"),
    coefficients = object$coefficients,
    varnames = object$varnames,
    base = object$base,
    n.vars = object$n.vars,
    npar = object$npar,
    coef.ind = cumsum(object$n.vars),
    nobs = nobs(object),
    parametrization = object$parametrization,
    optimization = object$optimization
  ), class = "summary_DirichletRegModel")

  names(object$coefficients) <- object$coefnames

  resid.mat <- round(t(apply(residuals(object, type = "standardized"), 2, quantile)), 4)
  colnames(resid.mat) <- c("Min", "1Q", "Median", "3Q", "Max")

  res$resid.mat <- resid.mat

  z.values <- object$coefficients / object$se
  p.values <- 2 * pnorm(-abs(z.values))
  coef.mat <- cbind(object$coefficients, object$se, z.values, p.values)
  colnames(coef.mat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  res$coef.mat <- coef.mat

  print(res, ...)
  invisible(res)

}



print.summary_DirichletRegModel <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  .wd <- getOption("width")

  signif_codes <- paste0("Significance codes: ", attr(symnum(0, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")), "legend"))

  if(interactive()) writeLines("")

  writeLines("Call:")
  writeLines(paste(strwrap(deparse(x$call, width.cutoff=500), .wd), sep="\n", collapse="\n"))

  cat("\nStandardized Residuals:\n")
  print(x$resid.mat, print.gap=2)
  cat("\n")

################################################################################
################################################################### COMMON PARAM

  if(x$parametrization == "common"){

    for(i in seq_len(length(x$varnames))){
      writeLines(paste0(rep("-", min(66L, .wd)), collapse=""))
      writeLines(paste0("Beta-Coefficients for variable no. ", i, ": ", x$varnames[i]))

      printCoefmat(x$coef.mat[ ifelse(i==1L, 1L, x$coef.ind[i-1L]+1L):x$coef.ind[i] , , drop = FALSE],
                   digits = digits, cs.ind=1:2, tst.ind=3, P.values = TRUE, signif.legend = FALSE)


    }
    writeLines(paste0(rep("-", min(66L, .wd)), collapse=""))
    writeLines(signif_codes)

  } else {
################################################################################
############################################################## ALTERNATIVE PARAM
    printed.var <- 1L
    set.size    <- x$n.vars[1L]

    cat("MEAN MODELS:\n",sep="",collapse="")

    for(i in seq_len(length(x$varnames))){
      if(i == x$base){
        writeLines(paste0(rep("-", min(66L, .wd)), collapse=""))
        writeLines(paste0("Coefficients for variable no. ",i,": ",x$varnames[i]))
        writeLines("- variable omitted (reference category) -")
      } else {
        writeLines(paste0(rep("-", min(66L, .wd)), collapse=""))
        writeLines(paste0("Coefficients for variable no. ",i,": ",x$varnames[i]))

        printCoefmat(x$coef.mat[printed.var:(printed.var+set.size-1),,drop=FALSE],
                     digits = digits, cs.ind=1:2, tst.ind=3, P.values = TRUE, signif.legend = FALSE)

        printed.var <- printed.var + set.size
      }
    }

    writeLines(paste0(rep("-", min(66L, .wd)), collapse=""))

    writeLines("")

    writeLines("PRECISION MODEL:")
    writeLines(paste0(rep("-", min(66L, .wd)), collapse=""))
    printCoefmat(x$coef.mat[printed.var:length(x$coefficients), , drop = FALSE],
                 digits = digits, cs.ind=1:2, tst.ind=3, P.values = TRUE, signif.legend = FALSE)
    writeLines(paste0(rep("-", min(66L, .wd)), collapse=""))
    writeLines(signif_codes)

  }

################################################################################
############################################################################ FIN

  writeLines("")
  writeLines(paste0("Log-likelihood: ",format(x$logLik,digits=digits)," on ",x$npar," df (", x$optimization$bfgs.it," BFGS + ",x$optimization$iterations," NR Iterations)"))
  writeLines(paste0("AIC: ", format(x$aic, digits=digits),", BIC: ", format(x$bic, digits=digits)))
  writeLines(paste0("Number of Observations: ",x$nobs))
  if(x$parametrization == "common"){
    writeLines(paste0("Link: Log\nParametrization: ", x$parametrization))
  } else {
    writeLines(paste0("Links: Logit (Means) and Log (Precision)\nParametrization: ", x$parametrization))
  }

  if(interactive()) writeLines("")

}
