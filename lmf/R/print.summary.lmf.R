print.summary.lmf <-
function(x,
                              digits = max(3, getOption("digits") - 3),
                              signif.stars = getOption("show.signif.stars"),
                              ...)
{
  #Print header
  cat("\nESTIMATING FLUCTUATING SELECTION IN AGE-STRUCTURED POPULATIONS\n",
      sep = "")
  #Print call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  #Print transition matrix
  cat("Projection matrix (l):\n")
  print.default(format(x$l, digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
  #Print lambda
  cat("Lambda:\n", paste(" ", format(x$lambda, digits = digits), sep = ""), "\n\n")
  #Print u and v
  cat("Stable age dist.(u) and reprod. values(v):\n")
  print.data.frame(format(data.frame(age = x$uage, u = x$u, v = x$v),
                          digits = digits), print.gap = 2, quote = FALSE, row.names = FALSE)
  cat("\n")
  #Print variance components
  cat("Variance components:\n")
  cat(" Environmental\n", paste(" ", format(x$sigma2.e, digits = digits),
                                sep = ""), "\n")
  cat(" Demographic\n")
  print.data.frame(format(data.frame(age = c(x$uage, "(total)"),
                                     estimate = c(unlist(x$sigma2.dj), x$sigma2.d), SD = c(unlist(x$sigma2.dj.sd),
                                                                                           x$sigma2.d.sd), df = c(unlist(x$sigma2.dj.dof), x$sigma2.d.dof)),
                          digits = digits), print.gap = 2, quote = FALSE,
                   row.names = FALSE)
  cat("\n")
  #Check if estimates by age and year should be printed
  if(x$what.level == "age")
  {
    #Print yearly alpha estimates by age
    cat("Yearly alpha estimates by age (ajt):\n")
    printCoefmat(x$coefficients.ajt, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
    cat("\n")
    #Print yearly covariance matrices by age
    cat("Yearly covariance matrix by age (Ajt):\n")
    sapply(1 : x$nage, function(i) {sapply(1 : x$nyear, function(j)
    {cat(paste("-age ", x$uage[i], ", year ", x$uyear[j], ":", sep = ""), "\n");
     print.default(format(x$Ajt[[i]][[j]], digits = digits), print.gap = 2,
                   quote = FALSE); cat("\n")})})
  }
  #Check if estimates by year should be printed
  if(x$what.level == "year" | x$what.level == "age")
  {
    #Print yearly alpha estimates
    cat("Yearly alpha estimates (at):\n")
    printCoefmat(x$coefficients.at, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
    cat("\n")
    #Print yearly alpha estimates corrected for sampling error
    cat("Yearly alpha estimates corrected for sampling error (atC):\n")
    printCoefmat(x$coefficients.atC, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
    cat("\n")
    #Print yearly covariance matrices (At + M)
    cat("Yearly covariance matrix (At + M):\n")
    sapply(1 : x$nyear, function(i)
    {cat(paste("-year ", x$uyear[i], ":", sep = ""), "\n");
     print.default(format(x$At[[i]], digits = digits), print.gap = 2,
                   quote = FALSE); cat("\n")})
  }
  #Print temporal alpha estimates (aM)
  if (length(x$coefficients.aM[, 1]))
  {
    cat("Temporal alpha estimates (a(M)):\n")
    printCoefmat(x$coefficients.aM, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
  }
  else cat("No temporal alpha estimates (a(M))\n")
  cat("\n")
  #Print temporal covariance matrix (M)
  if (length(x$M) > 1)
  {
    cat("Temporal covariance matrix (M):\n")
    print.default(format(x$M, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  else cat("No temporal covariance matrix (M)\n")
  cat("\n")
  #Print temporal alpha estimates assuming no fluctuating selection (M = 0)
  if (length(x$coefficients.anf[, 1]))
  {
    cat("Temporal alpha estimates assuming no fluct. selection (a(M=0)):\n")
    printCoefmat(x$coefficients.anf, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
  }
  else cat("No temporal alpha estimates assuming no fluct. selection (a(M=0))\n")
  cat("\n")
  #Print temporal covariance matrix assuming no fluctuating selection (M = 0)
  if (length(x$Anf) > 1)
  {
    cat("Temporal covariance matrix assuming no fluct. selection (A):\n")
    print.default(format(x$Anf, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  else cat("No temporal covariance matrix assuming no fluct. selection (A)\n")
  cat("\n")
  #End of print
  cat("-End-\n")
  invisible(x)
}
