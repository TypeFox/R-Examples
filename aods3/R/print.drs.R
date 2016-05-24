print.drs <- function(x, ...) {
	type <- x$type
  tab <- x$tab
  z <- switch(type, d = "Donner", rs = "Rao-Scott")
  cat("\n", "Type: ", z, "\n\n", sep = "")
	## computation of test details
  df <- nrow(tab) - 1   
  P <- pchisq(q = x$X2, df = df, lower.tail = FALSE)
  cat("N = ", sum(tab$N), " clusters, n = ",
      sum(tab$n), " subjects, m = ",
      sum(tab$m), " cases, I = ", nrow(tab), " groups.\n", sep = "")
  if(type == "d")
    cat("\nData and correction factors:\n")
  if(type == "rs")
    cat("\nData and design effects:\n")
  print(format(tab, digits = 4))
  cat("\nAdjusted chi2 test:\n")
  cat("X2 = ", round(x$X2, 1), ",
      df = ", df, ",
      P(> X2) = ", round(P, digits = 4), "\n", sep = "")
  if(type == "d")
    cat("\nIntra-cluster correlation (anova estimate):",
        round(x$rho, digits = 4), "\n")
  }
