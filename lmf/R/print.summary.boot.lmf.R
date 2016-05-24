print.summary.boot.lmf <-
function(x,
                                   digits = max(3, getOption("digits") - 3),
                                   signif.stars = getOption("show.signif.stars"),
                                   ...)
{
  #Title
  cat("\nESTIMATING FLUCTUATING SELECTION IN AGE-STRUCTURED POPULATIONS\n",
      sep = "")
  #Display number of bootstraps
  cat(sprintf(ngettext(x$nboot, "%s Bootstrap replicate generated\n",
                       "%s Bootstrap replicates generated\n"), paste(x$nboot)), sep = "")
  #Display call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  #Header
  cat("BOOTSTRAP SUMMARY STATISTICS:\n\n")
  #Print transition matrix and associated components
  if(!is.null(x$lest))
  {
    #Projection matrix
    cat("Projection matrix (l):\n")
    cat("-estimate:\n")
    print.default(format(x$lest, digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("\n")
    cat("-bias (estimate - bootstrap mean):\n")
    print.default(format(x$lbias, digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("\n")
    cat("-bootstrap mean:\n")
    print.default(format(x$lboot.mean, digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("\n")
    cat("-bootstrap st.dev.:\n")
    print.default(format(x$lboot.sd, digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("\n")
    #Lambda, u and v
    cat("Lambda, stable age dist.(u) and reprod. values(v):\n")
    print.default(format(x$luv, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n")
  }
  #Print the remaining parameters (all related to the alpha estimates)
  if(!is.null(x$sigma2.e) & x$call$what != "H0")
  {
    #Print variance components
    cat("Variance components:\n")
    #Environmental
    #Normal
    cat(" Environmental\n")
    print.data.frame(format(x$sigma2.e, digits = digits),
                     print.gap = 2, quote = FALSE, row.names = FALSE)
    #Demographic
    cat(" Demographic\n")
    print.data.frame(format(x$sigma2.dd, digits = digits),
                     print.gap = 2, quote = FALSE, row.names = FALSE)
    cat("\n")
    #Print alphas (aM) - Fluctuating selection
    cat("Temporal mean alpha (a(M)):\n")
    print.default(format(x$aM, digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("\n")
    #Print alphas covariance matrix (M) - Fluctuating selection
    cat("Temporal covariance matrix (M):\n")
    cat("-estimate:\n")
    print.default(format(x$Mest, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("-bias (estimate - bootstrap mean):\n")
    print.default(format(x$Mbias, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("-bootstrap mean:\n")
    print.default(format(x$Mboot.mean, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("-bootstrap st.dev.:\n")
    print.default(format(x$Mboot.sd, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n")
    #Print alphas (a(M=0)) - No fluctuating selection
    cat("Temporal alpha estimates assuming no fluct. selection (a(M=0)):\n")
    print.default(format(x$anf, digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("\n")
    #Print alphas covariance matrix (A) - No fluctuating selection
    cat("Covariance matrix assuming no fluct. selection (A):\n")
    cat("-estimate:\n")
    print.default(format(x$Anfest, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("-bias (estimate - bootstrap mean):\n")
    print.default(format(x$Anfbias, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("-bootstrap mean:\n")
    print.default(format(x$Anfboot.mean, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("-bootstrap st.dev.:\n")
    print.default(format(x$Anfboot.sd, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n")
  }
  #Print tests of significance
  if(length(x$coefficients.H0aMboot[, 1]) |
       length(x$coefficients.H0anfboot[, 1]) |
       length(x$coefficients.H0Mnfboot[, 1]))
  {
    cat("TESTS OF SIGNIFICANCE:\n")
    #Print tests of significance for alpha under fluctuating selection
    #H0: a=0|M
    if(length(x$coefficients.H0aMboot[, 1]))
    {
      cat("Test of directional selection under fluctuating selection (H0: a=[H0exp]|M):\n")
      printCoefmat(x$coefficients.H0aMboot, digits = digits, signif.stars = signif.stars,
                   na.print = "NA", ...)
      cat("\n")
    }
    #Print tests of significance for alpha under no fluctuating selection
    #H0: a=0|M=0
    if(length(x$coefficients.H0anfboot[, 1]))
    {
      cat("Test of directional selection assuming no fluctuating selection (H0: a=[H0exp]|M=0):\n")
      printCoefmat(x$coefficients.H0anfboot, digits = digits, signif.stars = signif.stars,
                   na.print = "NA", ...)
      cat("\n")
    }
    #Print tests of significance for alpha covariance matrix under the
    #assumption of directional selection and no fluctuating selection
    #H0: M=0|a
    if(length(x$coefficients.H0Mnfboot[, 1]))
    {
      cat("Test of fluctuating selection under directional selection (H0: M=[H0exp]|a):\n")
      printCoefmat(x$coefficients.H0Mnfboot, digits = digits, signif.stars = signif.stars,
                   na.print = "NA", ...)
      cat("\n")
    }
  }
  #Print bootstraps in a clear way if desired
  if(!is.null(x$lluvboot) & x$call$what != "H0")
  {
    #Print bootstrapped l, lambda, u and v
    cat("Boostrap replicates of projection matrix (l), lambda, u and v:\n")
    print.data.frame(format(as.data.frame(x$lluvboot), digits = digits), print.gap = 2,
                     quote = FALSE, row.names = FALSE)
    cat("\n")
  }
  if(!is.null(x$deboot) & x$call$what != "H0")
  {
    #Print bootstrapped sigma2.dj, sigma2.d and sigma2.e
    cat("Boostrap replicates of demograpic and environmeltal variances:\n")
    print.data.frame(format(as.data.frame(x$deboot), digits = digits), print.gap = 2,
                     quote = FALSE, row.names = FALSE)
    cat("\n")
    #Print bootstrapped at and At
    cat("Boostrap replicates of yearly alphas (at) and covariance matrices (At):\n")
    print.data.frame(format(as.data.frame(x$atAboot), digits = digits), print.gap = 2,
                     quote = FALSE, row.names = FALSE)
    cat("\n")
    #Print bootstrapped aM and M
    cat("Boostrap replicates of temporal alpha (a(M)) and covariance matrix (M):\n")
    print.data.frame(format(as.data.frame(x$aMMboot), digits = digits), print.gap = 2,
                     quote = FALSE, row.names = FALSE)
    cat("\n")
    #Print bootstrapped atC
    cat("Boostrap replicates of yearly alphas corrected for sampling error (atC):\n")
    print.data.frame(format(as.data.frame(x$atCboot), digits = digits), print.gap = 2,
                     quote = FALSE, row.names = FALSE)
    cat("\n")
    #Print bootstrapped anf and Anf
    cat("Boostrap replicates of temporal alpha (a(M=0)) and covariance matrix (A) assuming no fluct. selection:\n")
    print.data.frame(format(as.data.frame(x$anfAboot), digits = digits), print.gap = 2,
                     quote = FALSE, row.names = FALSE)
    cat("\n")
  }
  #The following prints bootstraps under a null hypothesis
  if(!is.null(x$H0aMboot))
  {
    #Print bootstrapped alphas under the null hypothesis given M != 0
    cat("Boostrap replicates of temporal alpha (a(M)) under the null hypothesis:\n")
    print.data.frame(format(as.data.frame(x$H0aMboot), digits = digits), print.gap = 2,
                     quote = FALSE, row.names = FALSE)
    cat("\n")
  }
  if(!is.null(x$H0anfboot))
  {
    #Print bootstrapped alphas under the null hypothesis given M != 0
    cat("Boostrap replicates of temporal alpha (a(M=0)) under the null hypothesis:\n")
    print.data.frame(format(as.data.frame(x$H0anfboot), digits = digits), print.gap = 2,
                     quote = FALSE, row.names = FALSE)
    cat("\n")
  }
  if(!is.null(x$H0atnfboot))
  {
    #Print bootstrapped yearly alphas under the null hypothesis given M = 0
    cat("Boostrap replicates of yearly alpha (at|a(M=0)) under the null hypothesis:\n")
    print.data.frame(format(as.data.frame(x$H0atnfboot), digits = digits), print.gap = 2,
                     quote = FALSE, row.names = FALSE)
    cat("\n")
    #Print bootstrapped M under the null hypothesis given M = 0
    cat("Boostrap replicates of temporal covariance matrix (M|a) under the null hypothesis:\n")
    print.data.frame(format(as.data.frame(x$H0Mnfboot), digits = digits), print.gap = 2,
                     quote = FALSE, row.names = FALSE)
    cat("\n")
  }
  #End of print
  cat("-End-\n")
  #Return 'x' invisibly
  invisible(x)
}
