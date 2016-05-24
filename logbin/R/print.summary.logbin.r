print.summary.logbin <- function(x, digits = max(3L, getOption("digits") - 3L), 
          signif.stars = getOption("show.signif.stars"), ...) 
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (!is.null(x$knots)) {
    cat("Smooth term knots:\n")
	for(smth in names(x$knots)) {
		cat(smth, ": ", sep = "")
		l <- length(x$knots[[smth]])
		if(l > 10) cat(x$knots[[smth]][1:5], "...", x$knots[[smth]][(l-4):l])
		else cat(x$knots[[smth]])
		cat("\n")
	}
	cat("\n")
  }
  cat("Deviance Residuals: \n")
  if (x$df.residual > 5) {
    x$deviance.resid <- setNames(quantile(x$deviance.resid, 
                                          na.rm = TRUE), c("Min", "1Q", "Median", "3Q", "Max"))
  }
  xx <- zapsmall(x$deviance.resid, digits + 1L)
  print.default(xx, digits = digits, na.print = "", print.gap = 2L)
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    df <- if ("df" %in% names(x)) 
      x[["df"]]
    else NULL
    if (!is.null(df) && (nsingular <- df[3L] - df[1L])) 
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
          sep = "")
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn, 
                                                               colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
  }
  cat("\n", apply(cbind(paste(format(c("Null", 
                                       "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance", 
                                                                                                        "deviance")]), digits = max(5L, digits + 1L)), " on", 
                                       format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"), 
                                       1L, paste, collapse = " "), sep = "")
  if (nzchar(mess <- naprint(x$na.action))) 
    cat("  (", mess, ")\n", sep = "")
  cat("\n", apply(cbind(paste(format(c("AIC:","AIC_c:"), justify = "right"), format(unlist(x[c("aic","aic.c")]), digits = max(4L, digits + 1L)),"\n")), 1L, paste, collapse = " "),  
      "\n", "Number of CEM iterations: ", x$iter, 
      "\n", sep = "")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      correl <- format(round(correl, 2L), nsmall = 2L, 
                       digits = digits)
      correl[!lower.tri(correl)] <- ""
      print(correl[-1, -p, drop = FALSE], quote = FALSE)
    }
  }
  cat("\n")
  invisible(x)
}